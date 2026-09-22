#!/01_software/miniconda3/bin/python3

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

try:
    import numpy as np
    from rapidfuzz import fuzz, process
except ImportError as e:
    sys.stderr.write(
        f"ERROR: missing dependency: {e}\n"
        "Install with:\n"
        "    pip install rapidfuzz numpy\n"
    )
    sys.exit(1)


DEFAULT_DATABASE = (
    "all.gene.csv"
)


def default_threads():
    """
    优先读取集群任务分配的 CPU 数。
    """
    for key in (
        "NSLOTS",
        "SLURM_CPUS_PER_TASK",
        "PBS_NP"
    ):
        value = os.environ.get(key)

        if value:
            try:
                return max(1, int(value))
            except ValueError:
                pass

    return min(4, os.cpu_count() or 1)


def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Fast fuzzy matching of a two-column query TSV "
            "against a two-column gene database."
        )
    )

    parser.add_argument(
        "query",
        help="Query file: two-column TSV"
    )

    parser.add_argument(
        "ratio",
        type=float,
        help="Minimum similarity threshold (0-100)"
    )

    parser.add_argument(
        "-d",
        "--database",
        default=DEFAULT_DATABASE,
        help=f"Database file [default: {DEFAULT_DATABASE}]"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="Output TSV [default: stdout]"
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="CPU threads"
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=1,
        help="Top matches retained for each query [default: 1]"
    )

    parser.add_argument(
        "--all-matches",
        action="store_true",
        help=(
            "Output every match >= threshold. "
            "This can generate very large output."
        )
    )

    parser.add_argument(
        "--scorer",
        choices=[
            "ratio",
            "token_sort",
            "token_set",
            "wratio"
        ],
        default="wratio",
        help=(
            "Similarity algorithm [default: wratio]. "
            "Use ratio for maximum speed."
        )
    )

    parser.add_argument(
        "--min-delta",
        type=float,
        default=5.0,
        help=(
            "If best score - second score < this value, "
            "mark as AMBIGUOUS [default: 5]"
        )
    )

    parser.add_argument(
        "--query-header",
        action="store_true",
        help="Query TSV has a header"
    )

    parser.add_argument(
        "--database-header",
        action="store_true",
        help="Database file has a header"
    )

    parser.add_argument(
        "--database-sep",
        default=",",
        help=(
            "Database delimiter [default: comma]. "
            "Use '\\t' for TSV."
        )
    )

    parser.add_argument(
        "--batch-size",
        type=int,
        default=500,
        help=(
            "Maximum number of query records processed "
            "per batch [default: 500]"
        )
    )

    parser.add_argument(
        "--max-matrix-mb",
        type=float,
        default=256.0,
        help=(
            "Maximum cdist score matrix memory "
            "per batch [default: 256 MB]"
        )
    )

    parser.add_argument(
        "--no-normalize",
        action="store_true",
        help="Disable text normalization"
    )

    return parser.parse_args()


def decode_sep(sep):

    if sep == r"\t":
        return "\t"

    return sep


def normalize_text(text):
    """
    标准化 annotation。

    例如：
        Wilms' tumor protein-1
        wilms tumor protein 1

    会转成相近形式。
    """

    text = text.strip().lower()

    text = text.replace(
        "_",
        " "
    )

    text = text.replace(
        "-",
        " "
    )

    text = re.sub(
        r"[^\w\s]",
        " ",
        text
    )

    text = re.sub(
        r"\s+",
        " ",
        text
    )

    return text.strip()


def identity(text):
    return text.strip()


def read_two_column_file(
    path,
    delimiter,
    has_header=False,
    label="input"
):

    records = []

    try:
        fh = open(
            path,
            "r",
            encoding="utf-8-sig",
            newline=""
        )

    except OSError as e:

        sys.stderr.write(
            f"ERROR: cannot open {label} file:\n"
            f"{path}\n"
            f"{e}\n"
        )

        sys.exit(1)

    with fh:

        reader = csv.reader(
            fh,
            delimiter=delimiter
        )

        start_line = (
            2 if has_header else 1
        )

        if has_header:
            next(
                reader,
                None
            )

        for line_no, row in enumerate(
            reader,
            start=start_line
        ):

            if not row:
                continue

            if row[0].lstrip().startswith("#"):
                continue

            if len(row) < 2:

                sys.stderr.write(
                    f"WARNING: skip {label} "
                    f"line {line_no}: "
                    f"fewer than 2 columns\n"
                )

                continue

            gene_id = row[0].strip()
            text = row[1].strip()

            if not gene_id or not text:

                sys.stderr.write(
                    f"WARNING: skip {label} "
                    f"line {line_no}: "
                    f"empty ID/text\n"
                )

                continue

            records.append(
                (
                    gene_id,
                    text
                )
            )

    return records


def get_scorer(name):

    scorers = {
        "ratio":
            fuzz.ratio,

        "token_sort":
            fuzz.token_sort_ratio,

        "token_set":
            fuzz.token_set_ratio,

        "wratio":
            fuzz.WRatio
    }

    return scorers[name]


def auto_batch_size(
    requested,
    max_mb,
    n_database
):
    """
    cdist 使用 float32：
    每个 similarity score = 4 bytes。

    自动控制矩阵大小，避免一次占用过多内存。
    """

    if n_database <= 0:
        return 1

    max_rows = int(
        (
            max_mb
            * 1024
            * 1024
        )
        /
        (
            n_database
            * 4
        )
    )

    return max(
        1,
        min(
            requested,
            max_rows
        )
    )


def sort_hit_indices(
    row,
    cutoff,
    top_n=None
):
    """
    从一行 similarity matrix 中快速取得 top hits。
    """

    if top_n is None:

        indices = np.flatnonzero(
            row >= cutoff
        )

    else:

        k = min(
            max(
                top_n,
                1
            ),
            row.size
        )

        if k == row.size:

            indices = np.arange(
                row.size
            )

        else:

            indices = np.argpartition(
                row,
                -k
            )[-k:]

        indices = indices[
            row[indices] >= cutoff
        ]

    if indices.size == 0:

        return indices

    # 第一排序条件：score 从高到低
    # 第二排序条件：database index 从小到大
    order = np.lexsort(
        (
            indices,
            -row[indices]
        )
    )

    return indices[order]


def main():

    args = parse_args()

    if not 0 <= args.ratio <= 100:

        sys.exit(
            "ERROR: ratio must be between 0 and 100"
        )

    if args.top_n < 1:

        sys.exit(
            "ERROR: --top-n must be >= 1"
        )

    if args.threads < 1:

        sys.exit(
            "ERROR: --threads must be >= 1"
        )

    if args.batch_size < 1:

        sys.exit(
            "ERROR: --batch-size must be >= 1"
        )

    if args.max_matrix_mb <= 0:

        sys.exit(
            "ERROR: --max-matrix-mb must be > 0"
        )

    if args.no_normalize:

        normalizer = identity

    else:

        normalizer = normalize_text

    scorer = get_scorer(
        args.scorer
    )

    #
    # 读取 query
    #

    query_records = read_two_column_file(
        args.query,
        delimiter="\t",
        has_header=args.query_header,
        label="query"
    )

    #
    # 读取 database
    #

    database_records = read_two_column_file(
        args.database,
        delimiter=decode_sep(
            args.database_sep
        ),
        has_header=args.database_header,
        label="database"
    )

    if not query_records:

        sys.exit(
            "ERROR: no valid query records"
        )

    if not database_records:

        sys.exit(
            "ERROR: no valid database records"
        )

    #
    # ==========================================
    # 核心优化 1
    #
    # database annotation 只 normalize 一次。
    #
    # 旧方法：
    #
    # query1 -> normalize database N 次
    # query2 -> normalize database N 次
    # ...
    #
    # 新方法：
    #
    # database -> normalize 一次
    #
    # ==========================================
    #

    database_normalized = [
        normalizer(text)
        for _, text
        in database_records
    ]

    #
    # ==========================================
    # 核心优化 2
    #
    # 建立 exact-match hash。
    #
    # annotation 完全一致时：
    #
    # 不需要进入 fuzzy matching。
    #
    # ==========================================
    #

    exact_index = defaultdict(
        list
    )

    for index, text in enumerate(
        database_normalized
    ):

        if text:

            exact_index[
                text
            ].append(
                index
            )

    #
    # 自动调整 batch
    #

    batch_size = auto_batch_size(
        args.batch_size,
        args.max_matrix_mb,
        len(
            database_records
        )
    )

    sys.stderr.write(
        "\n[INFO]\n"
        f"  queries      : {len(query_records)}\n"
        f"  database     : {len(database_records)}\n"
        f"  threads      : {args.threads}\n"
        f"  scorer       : {args.scorer}\n"
        f"  cutoff       : {args.ratio}\n"
        f"  batch size   : {batch_size}\n"
        f"  matrix limit : {args.max_matrix_mb:g} MB\n\n"
    )

    #
    # 输出文件
    #

    if args.output == "-":

        out_handle = sys.stdout
        close_output = False

    else:

        out_handle = open(
            args.output,
            "w",
            encoding="utf-8",
            newline=""
        )

        close_output = True

    writer = csv.writer(
        out_handle,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "query_id",
        "database_id",
        "query_text",
        "database_text",
        "similarity",
        "rank",
        "delta_to_second",
        "match_type",
        "status"
    ])

    #
    # 统计量
    #

    n_matched = 0
    n_unmatched = 0
    n_ambiguous = 0

    n_exact = 0
    n_fuzzy = 0

    n_output_hits = 0

    #
    # ==========================================
    # 分 batch 处理 query
    # ==========================================
    #

    for start in range(
        0,
        len(query_records),
        batch_size
    ):

        chunk = query_records[
            start:
            start + batch_size
        ]

        #
        # 保存每个 query 的候选
        #

        chunk_results = [
            None
        ] * len(chunk)

        #
        # 只把没有 exact match 的 query
        # 放入 fuzzy matching。
        #

        fuzzy_positions = []
        fuzzy_queries = []

        for position, (
            query_id,
            query_text
        ) in enumerate(
            chunk
        ):

            query_normalized = normalizer(
                query_text
            )

            exact_hits = exact_index.get(
                query_normalized
            )

            if exact_hits:

                if args.all_matches:

                    keep = exact_hits

                else:

                    #
                    # 内部至少保留两个 exact hit，
                    # 用于判断 ambiguity。
                    #

                    keep = exact_hits[
                        :max(
                            args.top_n,
                            2
                        )
                    ]

                chunk_results[position] = [
                    (
                        index,
                        100.0
                    )
                    for index in keep
                ]

            else:

                fuzzy_positions.append(
                    position
                )

                fuzzy_queries.append(
                    query_normalized
                )

        #
        # ======================================
        # 核心优化 3
        #
        # RapidFuzz C++ 批量 similarity matrix
        #
        # ======================================
        #

        if fuzzy_queries:

            score_matrix = process.cdist(
                fuzzy_queries,
                database_normalized,

                scorer=scorer,

                #
                # 已经提前 normalize，
                # 不再调用 Python processor。
                #
                processor=None,

                score_cutoff=args.ratio,

                dtype=np.float32,

                workers=args.threads
            )

            if args.all_matches:

                required_hits = None

            else:

                #
                # 即便最终只输出 top1，
                # 也内部计算 top2，
                # 用于 delta 判断。
                #

                required_hits = max(
                    args.top_n,
                    2
                )

            for matrix_row, position in enumerate(
                fuzzy_positions
            ):

                row = score_matrix[
                    matrix_row
                ]

                indices = sort_hit_indices(
                    row,
                    args.ratio,
                    required_hits
                )

                chunk_results[position] = [
                    (
                        int(index),
                        float(
                            row[index]
                        )
                    )
                    for index in indices
                ]

            #
            # 及时释放矩阵
            #

            del score_matrix

        #
        # ======================================
        # 按 query 原始顺序输出
        # ======================================
        #

        for position, (
            query_id,
            query_text
        ) in enumerate(
            chunk
        ):

            hits_all = (
                chunk_results[position]
                or []
            )

            if not hits_all:

                writer.writerow([
                    query_id,
                    ".",
                    query_text,
                    ".",
                    ".",
                    ".",
                    ".",
                    "NONE",
                    "UNMATCHED"
                ])

                n_unmatched += 1

                continue

            #
            # 最终输出多少个
            #

            if args.all_matches:

                hits = hits_all

            else:

                hits = hits_all[
                    :args.top_n
                ]

            #
            # 判断 exact / fuzzy
            #

            query_normalized = normalizer(
                query_text
            )

            exact_hits = exact_index.get(
                query_normalized
            )

            if exact_hits:

                match_type = "EXACT"

                n_exact += 1

            else:

                match_type = "FUZZY"

                n_fuzzy += 1

            n_matched += 1

            #
            # best - second best
            #

            if len(hits_all) >= 2:

                delta = (
                    hits_all[0][1]
                    -
                    hits_all[1][1]
                )

            else:

                delta = None

            ambiguous = (
                delta is not None
                and
                delta < args.min_delta
            )

            if ambiguous:

                n_ambiguous += 1

            #
            # 输出
            #

            for rank, (
                database_index,
                score
            ) in enumerate(
                hits,
                start=1
            ):

                database_id, database_text = (
                    database_records[
                        database_index
                    ]
                )

                if rank == 1:

                    if delta is None:

                        delta_text = "."

                    else:

                        delta_text = (
                            f"{delta:.2f}"
                        )

                    if ambiguous:

                        status = "AMBIGUOUS"

                    else:

                        status = "MATCHED"

                else:

                    delta_text = "."

                    status = "CANDIDATE"

                writer.writerow([
                    query_id,
                    database_id,
                    query_text,
                    database_text,
                    f"{score:.2f}",
                    rank,
                    delta_text,
                    match_type,
                    status
                ])

                n_output_hits += 1

    if close_output:

        out_handle.close()

    #
    # Summary
    #

    sys.stderr.write(
        "\n[SUMMARY]\n"
        f"  query genes       : {len(query_records)}\n"
        f"  matched            : {n_matched}\n"
        f"  exact matches      : {n_exact}\n"
        f"  fuzzy matches      : {n_fuzzy}\n"
        f"  unmatched          : {n_unmatched}\n"
        f"  ambiguous          : {n_ambiguous}\n"
        f"  output hits        : {n_output_hits}\n"
    )

    if args.output != "-":

        sys.stderr.write(
            f"  output             : {args.output}\n"
        )


if __name__ == "__main__":
    main()
