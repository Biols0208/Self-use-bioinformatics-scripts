#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# NCBI best genome downloader
#
# Purpose:
#   1) Query all annotated assemblies under a taxon.
#   2) Map assembly organism TaxIDs to species TaxIDs.
#   3) Optionally keep one "best" assembly per species.
#   4) Download user-selected NCBI Datasets file types.
#
# Requirements:
#   datasets, dataformat, python3, unzip, awk, sort, find
#
# NCBI Datasets include types:
#   genome,protein,cds,rna,gff3,gtf,gbff,seq-report,all,none
#
# Default best-assembly ranking:
#   RefSeq category (reference > representative > none)
#   > assembly level (complete > chromosome > scaffold > contig)
#   > RefSeq/GCF over GenBank/GCA
#   > scaffold N50
#   > contig N50
#   > annotation release date
#   > assembly release date
# ============================================================

VERSION="1.0.0"

TAXON=""
INCLUDE="protein,rna,cds,gff3"
OUTDIR="NCBI_best_genomes"
THREADS=20
ASSEMBLY_SOURCE="all"
ASSEMBLY_LEVEL=""
BEST_PER_SPECIES="yes"
EXCLUDE_ATYPICAL="yes"
EXCLUDE_MULTI_ISOLATE="yes"
METADATA_ONLY="no"
FORCE_METADATA="no"
FORCE_PACKAGE="no"
NO_PROGRESS="no"
API_KEY="${NCBI_API_KEY:-}"

usage() {
    cat <<'USAGE'
Usage:
  ncbi_best_genome_downloader.sh --taxid TAXID [options]

Required:
  -t, --taxid TAXID_OR_NAME
      NCBI TaxID or taxon name.
      Examples: 33208 (Metazoa), 7711 (Chordata), 7898 (Actinopterygii)

Main options:
  -i, --include TYPES
      Comma-separated file types to download.
      Allowed:
        genome,protein,cds,rna,gff3,gtf,gbff,seq-report,all,none
      Default:
        protein,rna,cds,gff3

  -o, --outdir DIR
      Output directory.
      Default: NCBI_best_genomes

  -j, --threads N
      Concurrent workers for datasets rehydrate.
      NCBI currently allows 1-30.
      Default: 20

  -s, --assembly-source SOURCE
      all | RefSeq | GenBank
      Default: all

  -l, --assembly-level LEVELS
      Optional comma-separated assembly levels:
        complete,chromosome,scaffold,contig
      Default: all levels

Selection options:
  --best-per-species yes|no
      yes: keep one best assembly per species (default)
      no : keep all matching annotated assemblies

  --include-atypical
      Include assemblies marked atypical by NCBI.
      Default: exclude atypical assemblies

  --include-multi-isolate
      Include assemblies from multi-isolate projects.
      Default: exclude multi-isolate assemblies

Execution options:
  --metadata-only
      Stop after generating selected_genomes.tsv and selected_accessions.txt.
      No sequence/annotation files are downloaded.

  --force-metadata
      Re-query NCBI metadata and redo species selection even if cached files exist.

  --force-package
      Recreate the dehydrated download package even if it already exists.

  --api-key KEY
      NCBI API key. Alternatively set environment variable NCBI_API_KEY.

  --no-progress
      Hide NCBI datasets progress bars.

  -h, --help
      Show help.

  --version
      Show script version.

Examples:

  # All annotated animals, one best genome/species; protein + RNA + CDS + GFF3
  ./ncbi_best_genome_downloader.sh \
      --taxid 33208 \
      --include protein,rna,cds,gff3 \
      --outdir NCBI_Metazoa \
      --threads 20

  # Ray-finned fishes, protein only
  ./ncbi_best_genome_downloader.sh \
      --taxid 7898 \
      --include protein \
      --outdir NCBI_Actinopterygii_protein

  # Chordates, genome + protein + CDS + RNA + GFF3 + GBFF
  ./ncbi_best_genome_downloader.sh \
      --taxid 7711 \
      --include genome,protein,cds,rna,gff3,gbff \
      --outdir NCBI_Chordata

  # RefSeq only
  ./ncbi_best_genome_downloader.sh \
      --taxid 33208 \
      --assembly-source RefSeq \
      --include protein,cds \
      --outdir NCBI_Metazoa_RefSeq

  # Only inspect the selected assemblies, do not download sequence files
  ./ncbi_best_genome_downloader.sh \
      --taxid 33208 \
      --metadata-only \
      --outdir NCBI_Metazoa_check

Notes:
  * "rna" means NCBI transcript/RNA FASTA, not strictly protein-coding mRNA only.
  * With --best-per-species yes, assemblies are grouped by NCBI species TaxID.
  * Existing metadata and dehydrated packages are reused by default.
USAGE
}

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 1
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

normalize_yes_no() {
    local x
    x="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
    case "$x" in
        yes|y|true|1) printf 'yes\n' ;;
        no|n|false|0) printf 'no\n' ;;
        *) die "Expected yes/no value, got: $1" ;;
    esac
}

validate_include() {
    local raw token
    raw="$(printf '%s' "$1" | tr -d '[:space:]')"
    [[ -n "$raw" ]] || die "--include cannot be empty"

    IFS=',' read -r -a tokens <<< "$raw"
    for token in "${tokens[@]}"; do
        case "$token" in
            genome|protein|cds|rna|gff3|gtf|gbff|seq-report|all|none) ;;
            *) die "Unsupported --include type: $token" ;;
        esac
    done

    if [[ "$raw" == *"all"* && "$raw" != "all" ]]; then
        die "Use --include all by itself, not together with other types"
    fi
    if [[ "$raw" == *"none"* && "$raw" != "none" ]]; then
        die "Use --include none by itself, not together with other types"
    fi

    INCLUDE="$raw"
}

validate_assembly_level() {
    local raw token
    raw="$(printf '%s' "$1" | tr -d '[:space:]' | tr '[:upper:]' '[:lower:]')"
    [[ -n "$raw" ]] || { ASSEMBLY_LEVEL=""; return; }

    IFS=',' read -r -a levels <<< "$raw"
    for token in "${levels[@]}"; do
        case "$token" in
            complete|chromosome|scaffold|contig) ;;
            *) die "Unsupported assembly level: $token" ;;
        esac
    done
    ASSEMBLY_LEVEL="$raw"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--taxid|--taxon)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            TAXON="$2"
            shift 2
            ;;
        -i|--include)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            validate_include "$2"
            shift 2
            ;;
        -o|--outdir)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            OUTDIR="$2"
            shift 2
            ;;
        -j|--threads)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            THREADS="$2"
            shift 2
            ;;
        -s|--assembly-source)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            ASSEMBLY_SOURCE="$2"
            shift 2
            ;;
        -l|--assembly-level)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            validate_assembly_level "$2"
            shift 2
            ;;
        --best-per-species)
            [[ $# -ge 2 ]] || die "$1 requires yes/no"
            BEST_PER_SPECIES="$(normalize_yes_no "$2")"
            shift 2
            ;;
        --include-atypical)
            EXCLUDE_ATYPICAL="no"
            shift
            ;;
        --include-multi-isolate)
            EXCLUDE_MULTI_ISOLATE="no"
            shift
            ;;
        --metadata-only)
            METADATA_ONLY="yes"
            shift
            ;;
        --force-metadata)
            FORCE_METADATA="yes"
            shift
            ;;
        --force-package)
            FORCE_PACKAGE="yes"
            shift
            ;;
        --api-key)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            API_KEY="$2"
            shift 2
            ;;
        --no-progress)
            NO_PROGRESS="yes"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --version)
            echo "$VERSION"
            exit 0
            ;;
        *)
            die "Unknown option: $1. Use --help."
            ;;
    esac
done

[[ -n "$TAXON" ]] || { usage >&2; die "--taxid/--taxon is required"; }

[[ "$THREADS" =~ ^[0-9]+$ ]] || die "--threads must be an integer"
(( THREADS >= 1 && THREADS <= 30 )) || die "--threads must be between 1 and 30"

case "$ASSEMBLY_SOURCE" in
    all|RefSeq|GenBank) ;;
    refseq) ASSEMBLY_SOURCE="RefSeq" ;;
    genbank) ASSEMBLY_SOURCE="GenBank" ;;
    *) die "--assembly-source must be all, RefSeq, or GenBank" ;;
esac

for cmd in datasets dataformat python3 unzip awk sort find wc; do
    command_exists "$cmd" || die "Required command not found: $cmd"
done

META_DIR="${OUTDIR}/metadata"
PKG_DIR="${OUTDIR}/package"
DATA_DIR="${PKG_DIR}/data"
LOG_DIR="${OUTDIR}/logs"

GENOME_JSONL="${META_DIR}/genomes.jsonl"
GENOME_TSV="${META_DIR}/genomes.raw.tsv"
ORGANISM_TAXIDS="${META_DIR}/organism_taxids.txt"
TAXONOMY_JSONL="${META_DIR}/taxonomy.jsonl"
TAXONOMY_TSV="${META_DIR}/taxonomy.tsv"
SELECTED_TSV="${META_DIR}/selected_genomes.tsv"
SELECTED_ACCESSIONS="${META_DIR}/selected_accessions.txt"
UNRESOLVED_TAXIDS="${META_DIR}/unresolved_taxids.txt"
PACKAGE_ZIP="${PKG_DIR}/selected_genomes.dehydrated.zip"
RUN_CONFIG="${META_DIR}/run_config.txt"

mkdir -p "$META_DIR" "$PKG_DIR" "$LOG_DIR"

DATASETS_API_ARGS=()
DOWNLOAD_COMMON=()
if [[ -n "$API_KEY" ]]; then
    DATASETS_API_ARGS+=(--api-key "$API_KEY")
    DOWNLOAD_COMMON+=(--api-key "$API_KEY")
fi
if [[ "$NO_PROGRESS" == "yes" ]]; then
    DOWNLOAD_COMMON+=(--no-progressbar)
fi

SUMMARY_FILTERS=(
    --annotated
    --assembly-version current
    --assembly-source "$ASSEMBLY_SOURCE"
    --limit all
    --as-json-lines
)

if [[ "$EXCLUDE_ATYPICAL" == "yes" ]]; then
    SUMMARY_FILTERS+=(--exclude-atypical)
fi
if [[ "$EXCLUDE_MULTI_ISOLATE" == "yes" ]]; then
    SUMMARY_FILTERS+=(--exclude-multi-isolate)
fi
if [[ -n "$ASSEMBLY_LEVEL" ]]; then
    SUMMARY_FILTERS+=(--assembly-level "$ASSEMBLY_LEVEL")
fi

SELECTION_CONFIG="${META_DIR}/selection_config.txt"
PACKAGE_CONFIG="${META_DIR}/package_config.txt"

EXPECTED_SELECTION_CONFIG="${META_DIR}/.selection_config.expected.$$"
cat > "$EXPECTED_SELECTION_CONFIG" <<EOF_SELECTION
taxon=${TAXON}
assembly_source=${ASSEMBLY_SOURCE}
assembly_level=${ASSEMBLY_LEVEL:-all}
best_per_species=${BEST_PER_SPECIES}
exclude_atypical=${EXCLUDE_ATYPICAL}
exclude_multi_isolate=${EXCLUDE_MULTI_ISOLATE}
EOF_SELECTION

if [[ -s "$SELECTION_CONFIG" ]] && ! cmp -s "$SELECTION_CONFIG" "$EXPECTED_SELECTION_CONFIG"; then
    log "Selection parameters changed; cached metadata will be rebuilt"
    FORCE_METADATA="yes"
    FORCE_PACKAGE="yes"
fi

if [[ "$FORCE_METADATA" == "yes" ]]; then
    FORCE_PACKAGE="yes"
fi

cp "$EXPECTED_SELECTION_CONFIG" "$SELECTION_CONFIG"
rm -f "$EXPECTED_SELECTION_CONFIG"

cat > "$RUN_CONFIG" <<EOF_CONFIG
script_version=${VERSION}
taxon=${TAXON}
include=${INCLUDE}
outdir=${OUTDIR}
threads=${THREADS}
assembly_source=${ASSEMBLY_SOURCE}
assembly_level=${ASSEMBLY_LEVEL:-all}
best_per_species=${BEST_PER_SPECIES}
exclude_atypical=${EXCLUDE_ATYPICAL}
exclude_multi_isolate=${EXCLUDE_MULTI_ISOLATE}
metadata_only=${METADATA_ONLY}
EOF_CONFIG

log "Taxon: $TAXON"
log "Include: $INCLUDE"
log "Assembly source: $ASSEMBLY_SOURCE"
log "Best per species: $BEST_PER_SPECIES"
log "Output: $OUTDIR"

# ------------------------------------------------------------
# STEP 1: Query annotated assembly metadata
# ------------------------------------------------------------
if [[ "$FORCE_METADATA" == "yes" || ! -s "$GENOME_JSONL" ]]; then
    log "[1/6] Querying annotated genome metadata from NCBI"
    rm -f "$GENOME_JSONL"
    datasets summary genome taxon "$TAXON" \
        "${SUMMARY_FILTERS[@]}" \
        "${DATASETS_API_ARGS[@]}" \
        > "$GENOME_JSONL"
else
    log "[1/6] Reusing genome metadata: $GENOME_JSONL"
fi

[[ -s "$GENOME_JSONL" ]] || die "Genome metadata query returned no records"

# ------------------------------------------------------------
# STEP 2: Convert assembly report to compact TSV
# ------------------------------------------------------------
if [[ "$FORCE_METADATA" == "yes" || ! -s "$GENOME_TSV" ]]; then
    log "[2/6] Converting genome metadata to TSV"
    dataformat tsv genome \
        --inputfile "$GENOME_JSONL" \
        --elide-header \
        --fields accession,organism-name,organism-tax-id,assminfo-level,assminfo-refseq-category,source_database,assminfo-release-date,assmstats-scaffold-n50,assmstats-contig-n50,annotinfo-release-date,annotinfo-name,assminfo-name,assmstats-total-sequence-len \
        > "$GENOME_TSV"
else
    log "[2/6] Reusing genome TSV: $GENOME_TSV"
fi

[[ -s "$GENOME_TSV" ]] || die "Genome TSV is empty"

awk -F '\t' 'NF >= 3 && $3 != "" {print $3}' "$GENOME_TSV" \
    | sort -u \
    > "$ORGANISM_TAXIDS"

[[ -s "$ORGANISM_TAXIDS" ]] || die "No organism TaxIDs were extracted"

# ------------------------------------------------------------
# STEP 3: Resolve organism TaxID -> species TaxID
# ------------------------------------------------------------
if [[ "$FORCE_METADATA" == "yes" || ! -s "$TAXONOMY_JSONL" ]]; then
    log "[3/6] Resolving organism TaxIDs to species TaxIDs"
    rm -f "$TAXONOMY_JSONL"
    datasets summary taxonomy taxon \
        --inputfile "$ORGANISM_TAXIDS" \
        --limit all \
        --as-json-lines \
        "${DATASETS_API_ARGS[@]}" \
        > "$TAXONOMY_JSONL"
else
    log "[3/6] Reusing taxonomy metadata: $TAXONOMY_JSONL"
fi

[[ -s "$TAXONOMY_JSONL" ]] || die "Taxonomy query returned no records"

if [[ "$FORCE_METADATA" == "yes" || ! -s "$TAXONOMY_TSV" ]]; then
    dataformat tsv taxonomy \
        --inputfile "$TAXONOMY_JSONL" \
        --template tax-summary \
        > "$TAXONOMY_TSV"
fi

[[ -s "$TAXONOMY_TSV" ]] || die "Taxonomy TSV is empty"

# ------------------------------------------------------------
# STEP 4: Select best assembly per species
# ------------------------------------------------------------
if [[ "$FORCE_METADATA" == "yes" || ! -s "$SELECTED_TSV" || ! -s "$SELECTED_ACCESSIONS" ]]; then
    log "[4/6] Selecting assemblies"

    python3 - "$GENOME_TSV" "$TAXONOMY_TSV" "$SELECTED_TSV" "$SELECTED_ACCESSIONS" "$UNRESOLVED_TAXIDS" "$BEST_PER_SPECIES" <<'PY'
import csv
import sys
from collections import defaultdict
from datetime import datetime

(
    genome_tsv,
    taxonomy_tsv,
    selected_tsv,
    selected_accessions,
    unresolved_taxids,
    best_per_species,
) = sys.argv[1:]


def clean(x):
    return "" if x is None else str(x).strip()


def to_int(x):
    x = clean(x).replace(",", "")
    if not x:
        return 0
    try:
        return int(float(x))
    except ValueError:
        return 0


def date_score(x):
    x = clean(x)
    if not x:
        return 0
    x = x[:10]
    try:
        d = datetime.strptime(x, "%Y-%m-%d")
        return d.year * 10000 + d.month * 100 + d.day
    except ValueError:
        return 0


def refseq_category_score(x):
    x = clean(x).lower()
    if "reference" in x:
        return 2
    if "representative" in x:
        return 1
    return 0


def assembly_level_score(x):
    x = clean(x).lower()
    return {
        "complete genome": 4,
        "complete": 4,
        "chromosome": 3,
        "scaffold": 2,
        "contig": 1,
    }.get(x, 0)


def source_score(accession, source):
    accession = clean(accession).upper()
    source = clean(source).upper()
    if accession.startswith("GCF_"):
        return 1
    if "REFSEQ" in source:
        return 1
    return 0


# Taxonomy table generated with:
# dataformat tsv taxonomy --template tax-summary
with open(taxonomy_tsv, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    required = {"Taxid", "Tax name", "Rank", "Species name", "Species taxid"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(
            "taxonomy.tsv is missing expected columns: " + ", ".join(sorted(missing))
        )

    taxmap = {}
    for row in reader:
        taxid = clean(row.get("Taxid"))
        rank = clean(row.get("Rank")).upper()
        species_taxid = clean(row.get("Species taxid"))
        species_name = clean(row.get("Species name"))

        # Defensive fallback if a future report leaves species fields empty
        # for records already at species rank.
        if not species_taxid and rank == "SPECIES":
            species_taxid = taxid
            species_name = clean(row.get("Tax name"))

        if taxid:
            taxmap[taxid] = {
                "species_taxid": species_taxid,
                "species_name": species_name,
                "rank": rank,
                "tax_name": clean(row.get("Tax name")),
            }


# genome.raw.tsv column order, fixed by the --fields list in the shell script.
# 0 accession
# 1 organism-name
# 2 organism-tax-id
# 3 assminfo-level
# 4 assminfo-refseq-category
# 5 source_database
# 6 assminfo-release-date
# 7 assmstats-scaffold-n50
# 8 assmstats-contig-n50
# 9 annotinfo-release-date
# 10 annotinfo-name
# 11 assminfo-name
# 12 assmstats-total-sequence-len
species_groups = defaultdict(list)
unresolved = set()
all_records = []

with open(genome_tsv, newline="", encoding="utf-8") as fh:
    reader = csv.reader(fh, delimiter="\t")
    for row in reader:
        if len(row) < 13:
            continue

        (
            accession,
            organism_name,
            organism_taxid,
            assembly_level,
            refseq_category,
            source_database,
            assembly_release_date,
            scaffold_n50,
            contig_n50,
            annotation_release_date,
            annotation_name,
            assembly_name,
            total_sequence_length,
        ) = row[:13]

        accession = clean(accession)
        organism_taxid = clean(organism_taxid)
        if not accession or not organism_taxid:
            continue

        tx = taxmap.get(organism_taxid)
        if tx is None or not tx["species_taxid"]:
            unresolved.add(organism_taxid)
            continue

        species_taxid = tx["species_taxid"]
        species_name = tx["species_name"] or clean(organism_name)

        rec = {
            "species_taxid": species_taxid,
            "species_name": species_name,
            "organism_taxid": organism_taxid,
            "organism_name": clean(organism_name),
            "accession": accession,
            "source_database": clean(source_database),
            "refseq_category": clean(refseq_category),
            "assembly_level": clean(assembly_level),
            "scaffold_n50": to_int(scaffold_n50),
            "contig_n50": to_int(contig_n50),
            "total_sequence_length": to_int(total_sequence_length),
            "annotation_release_date": clean(annotation_release_date),
            "assembly_release_date": clean(assembly_release_date),
            "annotation_name": clean(annotation_name),
            "assembly_name": clean(assembly_name),
        }

        rec["refseq_category_score"] = refseq_category_score(refseq_category)
        rec["assembly_level_score"] = assembly_level_score(assembly_level)
        rec["source_score"] = source_score(accession, source_database)

        rec["_score"] = (
            rec["refseq_category_score"],
            rec["assembly_level_score"],
            rec["source_score"],
            rec["scaffold_n50"],
            rec["contig_n50"],
            date_score(rec["annotation_release_date"]),
            date_score(rec["assembly_release_date"]),
            accession,
        )

        species_groups[species_taxid].append(rec)
        all_records.append(rec)

selected = []

if best_per_species == "yes":
    for species_taxid, candidates in species_groups.items():
        candidates.sort(key=lambda r: r["_score"], reverse=True)
        best = candidates[0]
        best["candidate_count"] = len(candidates)
        selected.append(best)
else:
    for species_taxid, candidates in species_groups.items():
        candidates.sort(key=lambda r: r["_score"], reverse=True)
        n = len(candidates)
        for rec in candidates:
            rec["candidate_count"] = n
            selected.append(rec)

selected.sort(
    key=lambda r: (
        r["species_name"].lower(),
        r["species_taxid"],
        r["accession"],
    )
)

fields = [
    "species_taxid",
    "species_name",
    "organism_taxid",
    "organism_name",
    "accession",
    "source_database",
    "refseq_category",
    "assembly_level",
    "scaffold_n50",
    "contig_n50",
    "total_sequence_length",
    "annotation_release_date",
    "assembly_release_date",
    "annotation_name",
    "assembly_name",
    "candidate_count",
    "refseq_category_score",
    "assembly_level_score",
    "source_score",
]

with open(selected_tsv, "w", newline="", encoding="utf-8") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for rec in selected:
        writer.writerow(rec)

with open(selected_accessions, "w", encoding="utf-8") as out:
    for rec in selected:
        out.write(rec["accession"] + "\n")

with open(unresolved_taxids, "w", encoding="utf-8") as out:
    for taxid in sorted(unresolved):
        if taxid:
            out.write(taxid + "\n")

print(f"Genome records parsed: {len(all_records):,}")
print(f"Species represented:    {len(species_groups):,}")
print(f"Assemblies selected:    {len(selected):,}")
print(f"Unresolved TaxIDs:      {len(unresolved):,}")
PY
else
    log "[4/6] Reusing selected assembly list: $SELECTED_ACCESSIONS"
fi

[[ -s "$SELECTED_ACCESSIONS" ]] || die "No assemblies were selected"

SELECTED_COUNT="$(wc -l < "$SELECTED_ACCESSIONS" | tr -d '[:space:]')"
log "Selected assemblies: $SELECTED_COUNT"
log "Selection table: $SELECTED_TSV"

if [[ "$METADATA_ONLY" == "yes" ]]; then
    log "Metadata-only mode requested; stopping before download"
    exit 0
fi

# ------------------------------------------------------------
# STEP 5: Create dehydrated package
# ------------------------------------------------------------
ACCESSION_CKSUM="$(cksum "$SELECTED_ACCESSIONS" | awk '{print $1 ":" $2}')"
EXPECTED_PACKAGE_CONFIG="${META_DIR}/.package_config.expected.$$"
cat > "$EXPECTED_PACKAGE_CONFIG" <<EOF_PACKAGE
include=${INCLUDE}
accession_cksum=${ACCESSION_CKSUM}
EOF_PACKAGE

if [[ -s "$PACKAGE_CONFIG" ]] && ! cmp -s "$PACKAGE_CONFIG" "$EXPECTED_PACKAGE_CONFIG"; then
    log "Download file types or selected accessions changed; package will be rebuilt"
    FORCE_PACKAGE="yes"
fi

if [[ "$FORCE_PACKAGE" == "yes" || ! -s "$PACKAGE_ZIP" ]]; then
    log "[5/6] Creating dehydrated package"
    rm -f "$PACKAGE_ZIP"

    datasets download genome accession \
        --inputfile "$SELECTED_ACCESSIONS" \
        --include "$INCLUDE" \
        --dehydrated \
        --filename "$PACKAGE_ZIP" \
        "${DOWNLOAD_COMMON[@]}"

    cp "$EXPECTED_PACKAGE_CONFIG" "$PACKAGE_CONFIG"
else
    log "[5/6] Reusing dehydrated package: $PACKAGE_ZIP"
fi
rm -f "$EXPECTED_PACKAGE_CONFIG"

[[ -s "$PACKAGE_ZIP" ]] || die "Failed to create dehydrated package"

# ------------------------------------------------------------
# STEP 6: Extract and rehydrate
# ------------------------------------------------------------
log "[6/6] Extracting and rehydrating selected files"
mkdir -p "$DATA_DIR"
unzip -q -o "$PACKAGE_ZIP" -d "$DATA_DIR"

# Re-running rehydrate is intentional: it provides resume behavior for
# incomplete downloads in the existing extracted package.
if [[ "$INCLUDE" != "none" ]]; then
    datasets rehydrate \
        --directory "$DATA_DIR" \
        --gzip \
        --max-workers "$THREADS" \
        "${DOWNLOAD_COMMON[@]}"
else
    log "--include none: no sequence/annotation files to rehydrate"
fi

# ------------------------------------------------------------
# Build audit file lists.
# Only inspect directories in the CURRENT selected_accessions.txt so that
# files left by an older run in the same output directory cannot contaminate
# downstream file lists.
# ------------------------------------------------------------
ROOT="${DATA_DIR}/ncbi_dataset/data"
if [[ -d "$ROOT" ]]; then
    for list_name in \
        all_downloaded_files.list \
        protein_files.list \
        rna_files.list \
        cds_files.list \
        gff3_files.list \
        gtf_files.list \
        gbff_files.list \
        genome_files.list \
        seq_report_files.list \
        missing_assembly_dirs.list
    do
        : > "${META_DIR}/${list_name}"
    done

    while IFS= read -r accession; do
        [[ -n "$accession" ]] || continue
        assembly_dir="${ROOT}/${accession}"

        if [[ ! -d "$assembly_dir" ]]; then
            printf '%s\n' "$accession" >> "${META_DIR}/missing_assembly_dirs.list"
            continue
        fi

        find "$assembly_dir" -type f >> "${META_DIR}/all_downloaded_files.list"
        find "$assembly_dir" -type f -name 'protein.faa.gz' >> "${META_DIR}/protein_files.list"
        find "$assembly_dir" -type f -name 'rna.fna.gz' >> "${META_DIR}/rna_files.list"
        find "$assembly_dir" -type f \( -name 'cds.fna.gz' -o -name '*cds_from_genomic.fna.gz' -o -name '*cds*.fna.gz' \) >> "${META_DIR}/cds_files.list"
        find "$assembly_dir" -type f -name 'genomic.gff.gz' >> "${META_DIR}/gff3_files.list"
        find "$assembly_dir" -type f -name 'genomic.gtf.gz' >> "${META_DIR}/gtf_files.list"
        find "$assembly_dir" -type f -name 'genomic.gbff.gz' >> "${META_DIR}/gbff_files.list"
        find "$assembly_dir" -type f -name '*_genomic.fna.gz' >> "${META_DIR}/genome_files.list"
        find "$assembly_dir" -type f -name 'sequence_report.jsonl.gz' >> "${META_DIR}/seq_report_files.list"
    done < "$SELECTED_ACCESSIONS"

    for list_file in "${META_DIR}"/*_files.list "${META_DIR}/missing_assembly_dirs.list"; do
        [[ -f "$list_file" ]] || continue
        sort -u "$list_file" -o "$list_file"
    done
fi

log "Finished"
printf '\n%-28s %s\n' "Selected assemblies:" "$SELECTED_COUNT"
printf '%-28s %s\n' "Selected metadata:" "$SELECTED_TSV"
printf '%-28s %s\n' "Selected accessions:" "$SELECTED_ACCESSIONS"
printf '%-28s %s\n' "Downloaded data root:" "$ROOT"
printf '%-28s %s\n' "Run configuration:" "$RUN_CONFIG"

if [[ -d "$ROOT" ]]; then
    printf '\nFile counts (0 is valid if that type was not requested or unavailable):\n'
    for f in \
        protein_files.list \
        rna_files.list \
        cds_files.list \
        gff3_files.list \
        gtf_files.list \
        gbff_files.list \
        genome_files.list \
        seq_report_files.list
    do
        path="${META_DIR}/${f}"
        count=0
        [[ -f "$path" ]] && count="$(wc -l < "$path" | tr -d '[:space:]')"
        printf '  %-24s %s\n' "$f" "$count"
    done
fi
