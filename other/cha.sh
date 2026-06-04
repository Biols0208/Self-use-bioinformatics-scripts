#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
NCBI Datasets -> TSV. Uses:
  datasets summary ... --as-json-lines | dataformat tsv ... --fields ... > out.tsv

Examples:
  # Gene by symbol (default mode), restrict to species:
  ### Adding `gene-taxon` as a qualifier is highly recommended; otherwise, the output file will be very large.

  cha.sh gene --gene BRCA1,TP53 --gene-taxon human --out gene.tsv --fields go-bp-id,go-bp-name,go-mf-id,go-mf-name,gene-id,gene-type

  # Gene by gene-id:
  cha.sh gene --gene-mode gene-id --gene 672,2597 --out gene.tsv

  # Genome by taxon:
  cha.sh genome --taxon "mus musculus" --out genome.tsv

  # Genome by accession (assembly or BioProject):
  cha.sh genome --accession GCF_000001405.40 --out genome.tsv

Common options:
  --out <file.tsv>            Output TSV path (required)
  --limit <n>                 datasets summary --limit (default: 20)
  --fields <csv>              Override default fields for this entity
  --datasets <path>           Path to `datasets` executable (optional; default: from PATH)
  --dataformat <path>         Path to `dataformat` executable (optional; default: from PATH)
  --api-key <key>             NCBI Datasets API key (optional; default: env NCBI_API_KEY; if none, omit)

gene subcommand options:
  --gene <csv>                Query values (required). e.g. BRCA1,TP53 or 672,2597
  --gene-mode <mode>          symbol|gene-id|accession|taxon|locus-tag (default: symbol)
  --gene-taxon <taxon>        Optional restriction for symbol/locus-tag (recommended)

genome subcommand options (choose one):
  --taxon <name|taxid>        genome taxon query
  --accession <acc>           genome accession query

Default fields (can be overridden with --fields):
  gene:     annotation-assembly-accession,annotation-genomic-range-accession,annotation-genomic-range-range-start,annotation-genomic-range-range-stop,annotation-release-name,chromosomes,common-name,ensembl-geneids,gene-id,gene-type,go-bp-id,go-bp-name,go-bp-reference-pmid,go-cc-id,go-cc-name,go-mf-id,go-mf-name,group-id,summary-description,swissprot-accessions,symbol,synonyms,tax-id,tax-name,transcript-count
  genome:   accession,annotinfo-busco-complete,annotinfo-busco-lineage,annotinfo-featcount-gene-non-coding,annotinfo-featcount-gene-other,annotinfo-featcount-gene-protein-coding,annotinfo-featcount-gene-pseudogene,annotinfo-featcount-gene-total,annotinfo-release-date,annotinfo-report-url,assminfo-assembly-method,assminfo-atypicalwarnings,assminfo-bioproject-lineage-accession,assminfo-biosample-accession,assminfo-biosample-description-organism-common-name,assminfo-biosample-description-organism-tax-id,assminfo-name,assminfo-paired-assm-accession,assminfo-paired-assm-status,assminfo-refseq-category,assminfo-type,assmstats-atgc-count,assmstats-contig-n50,assmstats-gc-percent,assmstats-genome-coverage,assmstats-scaffold-n50,assmstats-total-number-of-chromosomes,assmstats-total-sequence-len,assmstats-total-ungapped-len,current-accession,organelle-total-seq-length,organism-common-name,organism-name,organism-tax-id

Notes:
- If dataformat fails due to field name mismatch, the script will exit non-zero.
  (Because TSV is the only output; you can override --fields to fix.)
EOF
}

die(){ echo "Error: $*" >&2; exit 1; }

resolve_bin() {
  # Usage: resolve_bin <flag-name> <value-or-empty> <command-name>
  local flag="$1"
  local val="${2:-}"
  local cmd="$3"

  if [[ -n "$val" ]]; then
    [[ -f "$val" ]] || die "$flag not found: $val"
    [[ -x "$val" ]] || die "$flag is not executable: $val"
    printf '%s\n' "$val"
  else
    local found=""
    found="$(command -v "$cmd" 2>/dev/null || true)"
    [[ -n "$found" ]] || die "Missing dependency in PATH: $cmd (or pass $flag /path/to/$cmd)"
    printf '%s\n' "$found"
  fi
}

split_csv_to_args() {
  local s="${1:-}"
  [[ -z "$s" ]] && return 0
  awk -v s="$s" 'BEGIN{n=split(s,a,","); for(i=1;i<=n;i++){gsub(/^[ \t]+|[ \t]+$/,"",a[i]); if(length(a[i])) print a[i];}}'
}

[[ $# -lt 1 ]] && { usage; exit 1; }
SUB="$1"; shift

LIMIT=20
OUT=""

# entity-specific defaults
FIELDS_GENE_DEFAULT="annotation-assembly-accession,annotation-genomic-range-accession,annotation-genomic-range-range-start,annotation-genomic-range-range-stop,annotation-release-name,chromosomes,common-name,ensembl-geneids,gene-id,gene-type,go-bp-id,go-bp-name,go-bp-reference-pmid,go-cc-id,go-cc-name,go-mf-id,go-mf-name,group-id,summary-description,swissprot-accessions,symbol,synonyms,tax-id,tax-name,transcript-count"
FIELDS_GENOME_DEFAULT="accession,annotinfo-busco-complete,annotinfo-busco-lineage,annotinfo-featcount-gene-non-coding,annotinfo-featcount-gene-other,annotinfo-featcount-gene-protein-coding,annotinfo-featcount-gene-pseudogene,annotinfo-featcount-gene-total,annotinfo-release-date,annotinfo-report-url,assminfo-assembly-method,assminfo-atypicalwarnings,assminfo-bioproject-lineage-accession,assminfo-biosample-accession,assminfo-biosample-description-organism-common-name,assminfo-biosample-description-organism-tax-id,assminfo-name,assminfo-paired-assm-accession,assminfo-paired-assm-status,assminfo-refseq-category,assminfo-type,assmstats-atgc-count,assmstats-contig-n50,assmstats-gc-percent,assmstats-genome-coverage,assmstats-scaffold-n50,assmstats-total-number-of-chromosomes,assmstats-total-sequence-len,assmstats-total-ungapped-len,current-accession,organelle-total-seq-length,organism-common-name,organism-name,organism-tax-id"
FIELDS_TAX_DEFAULT="tax-id,scientific-name,rank"

# shared inputs
TAXON=""
ACCESSION=""
FIELDS_OVERRIDE=""

# gene inputs
GENE_MODE="symbol"
GENE_QUERY=""
GENE_TAXON=""

# explicit bins (optional)
DATASETS_PATH="/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/datasets"
DATAFORMAT_PATH="/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/dataformat"

# api key (optional)
API_KEY="a99aaa85c59305e417219ef841f5fff89408"

# parse options (common + subcommand-specific)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) OUT="${2:-}"; shift 2;;
    --limit) LIMIT="${2:-}"; shift 2;;
    --fields) FIELDS_OVERRIDE="${2:-}"; shift 2;;

    --datasets) DATASETS_PATH="${2:-}"; shift 2;;
    --dataformat) DATAFORMAT_PATH="${2:-}"; shift 2;;
    --api-key) API_KEY="${2:-}"; shift 2;;

    --taxon) TAXON="${2:-}"; shift 2;;
    --accession) ACCESSION="${2:-}"; shift 2;;

    --gene-mode) GENE_MODE="${2:-}"; shift 2;;
    --gene) GENE_QUERY="${2:-}"; shift 2;;
    --gene-taxon) GENE_TAXON="${2:-}"; shift 2;;

    --help|-h) usage; exit 0;;
    *) die "Unknown option: $1 (use --help)";;
  esac
done

[[ -z "$OUT" ]] && die "--out <file.tsv> is required"

# Resolve binaries
DATASETS_BIN="$(resolve_bin --datasets "$DATASETS_PATH" datasets)"
DATAFORMAT_BIN="$(resolve_bin --dataformat "$DATAFORMAT_PATH" dataformat)"

# Resolve API key: flag overrides env; if still empty, omit.
if [[ -z "$API_KEY" ]]; then
  API_KEY="${NCBI_API_KEY:-}"
fi

# Build optional api-key args (as an array to preserve quoting safely)
API_KEY_ARGS=()
if [[ -n "$API_KEY" ]]; then
  API_KEY_ARGS=(--api-key "$API_KEY")
fi

run_pipe() {
  # $1: entity for dataformat (gene|genome|taxonomy)
  # remaining args: full datasets summary command args AFTER 'datasets summary'
  local entity="$1"; shift
  local fields="$1"; shift

  "$DATASETS_BIN" summary "$@" --limit "$LIMIT" --as-json-lines "${API_KEY_ARGS[@]}" \
    | "$DATAFORMAT_BIN" tsv "$entity" --fields "$FIELDS" \
    | awk '!a[$0]++' > "$OUT"
}

case "$SUB" in
  gene)
    [[ -z "$GENE_QUERY" ]] && die "gene requires --gene <csv>"
    mapfile -t qlist < <(split_csv_to_args "$GENE_QUERY")
    [[ ${#qlist[@]} -eq 0 ]] && die "--gene provided but empty"

    FIELDS="${FIELDS_OVERRIDE:-$FIELDS_GENE_DEFAULT}"

    case "$GENE_MODE" in
      gene-id|accession|taxon)
        run_pipe gene "$FIELDS" gene "$GENE_MODE" "${qlist[@]}"
        ;;
      symbol|locus-tag)
        if [[ -n "$GENE_TAXON" ]]; then
          "$DATASETS_BIN" summary gene "$GENE_MODE" "${qlist[@]}" --taxon "$GENE_TAXON" --limit "$LIMIT" --as-json-lines "${API_KEY_ARGS[@]}" \
            | "$DATAFORMAT_BIN" tsv gene --fields "$FIELDS" \
            | awk '!a[$0]++' > "$OUT"
        else
          run_pipe gene "$FIELDS" gene "$GENE_MODE" "${qlist[@]}"
        fi
        ;;
      *)
        die "Invalid --gene-mode: $GENE_MODE (symbol|gene-id|accession|taxon|locus-tag)"
        ;;
    esac
    ;;
  genome)
    FIELDS="${FIELDS_OVERRIDE:-$FIELDS_GENOME_DEFAULT}"

    if [[ -n "$TAXON" && -n "$ACCESSION" ]]; then
      die "genome: choose only one of --taxon OR --accession"
    elif [[ -n "$TAXON" ]]; then
      run_pipe genome "$FIELDS" genome taxon "$TAXON"
    elif [[ -n "$ACCESSION" ]]; then
      run_pipe genome "$FIELDS" genome accession "$ACCESSION"
    else
      die "genome requires --taxon <...> OR --accession <...>"
    fi
    ;;
  taxonomy)
    [[ -z "$TAXON" ]] && die "taxonomy requires --taxon <...>"
    FIELDS="${FIELDS_OVERRIDE:-$FIELDS_TAX_DEFAULT}"
    run_pipe taxonomy "$FIELDS" taxonomy taxon "$TAXON"
    ;;
  *)
    die "Unknown subcommand: $SUB (use: gene|genome|taxonomy)"
    ;;
esac

echo "Wrote TSV: $OUT"
