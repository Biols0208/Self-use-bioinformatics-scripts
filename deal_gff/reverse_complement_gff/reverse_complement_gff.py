#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
from pathlib import Path

def open_maybe_gz(path, mode='rt'):
    path = str(path)
    if path.endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode, encoding='utf-8')

def read_fasta_lengths(fasta_path):
    """
    读取FASTA并返回 {seqid: length}
    """
    lengths = {}
    seqid = None
    seqlen = 0

    with open_maybe_gz(fasta_path, 'rt') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if seqid is not None:
                    lengths[seqid] = seqlen
                # header 取第一个空白前作为 seqid
                seqid = line[1:].split()[0]
                seqlen = 0
            else:
                seqlen += len(line.strip())
        if seqid is not None:
            lengths[seqid] = seqlen

    return lengths

def transform_coords(start, end, chrom_len):
    """
    1-based闭区间坐标镜像:
    new_start = L - end + 1
    new_end   = L - start + 1
    """
    new_start = chrom_len - end + 1
    new_end = chrom_len - start + 1
    return new_start, new_end

def flip_strand(strand):
    if strand == '+':
        return '-'
    elif strand == '-':
        return '+'
    return strand  # '.', '?' 等保持不变

def process_gff(gff_in, gff_out, fasta_lengths, rc_chroms):
    """
    rc_chroms: 需要反向互补并修正GFF的染色体集合 set
    """
    with open_maybe_gz(gff_in, 'rt') as fin, open_maybe_gz(gff_out, 'wt') as fout:
        for line in fin:
            if line.startswith('#') or not line.strip():
                fout.write(line)
                continue

            cols = line.rstrip('\n').split('\t')
            if len(cols) != 9:
                # 非标准行原样输出
                fout.write(line)
                continue

            seqid, source, ftype, start, end, score, strand, phase, attributes = cols

            if seqid in rc_chroms:
                if seqid not in fasta_lengths:
                    raise ValueError(f"GFF中的染色体 {seqid} 在FASTA里找不到长度。")
                chrom_len = fasta_lengths[seqid]

                s = int(start)
                e = int(end)
                if s > e:
                    raise ValueError(f"发现start>end: {line.strip()}")

                new_s, new_e = transform_coords(s, e, chrom_len)
                new_strand = flip_strand(strand)

                cols[3] = str(new_s)
                cols[4] = str(new_e)
                cols[6] = new_strand

            fout.write('\t'.join(cols) + '\n')

def main():
    parser = argparse.ArgumentParser(
        description="对指定染色体进行GFF坐标镜像和链方向翻转（用于染色体反向互补后注释修正）"
    )
    parser.add_argument('--fasta', required=True, help='基因组FASTA文件（可.gz）')
    parser.add_argument('--gff-in', required=True, help='输入GFF/GTF文件（可.gz）')
    parser.add_argument('--gff-out', required=True, help='输出GFF/GTF文件（可.gz）')
    parser.add_argument(
        '--rc-chroms', required=True,
        help='需要反向互补的染色体名，逗号分隔，如 chr1,chr3,scaffold_8'
    )

    args = parser.parse_args()

    rc_chroms = {x.strip() for x in args.rc_chroms.split(',') if x.strip()}
    if not rc_chroms:
        raise ValueError("--rc-chroms 不能为空")

    fasta_lengths = read_fasta_lengths(args.fasta)
    process_gff(args.gff_in, args.gff_out, fasta_lengths, rc_chroms)

    print("Done.")
    print(f"反向互补修正染色体: {', '.join(sorted(rc_chroms))}")
    print(f"输出文件: {args.gff_out}")

if __name__ == '__main__':
    main()
