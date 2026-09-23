#!/usr/bin/env bash

set -euo pipefail

############################################
# Configuration
############################################

OUTDIR="${1:-NCBI_Metazoa_best_proteomes}"
THREADS="${2:-20}"

METAZOA_TAXID=33208

META="${OUTDIR}/metadata"
PACKAGE="${OUTDIR}/package"
PROTEOMES="${OUTDIR}/proteomes"

mkdir -p "${META}"
mkdir -p "${PACKAGE}"
mkdir -p "${PROTEOMES}"


############################################
# Check programs
############################################

for cmd in datasets dataformat python3 unzip
do
    if ! command -v "${cmd}" >/dev/null 2>&1
    then
        echo "[ERROR] Cannot find: ${cmd}" >&2
        exit 1
    fi
done


############################################
# STEP 1
# Get all annotated animal genome metadata
############################################

echo
echo "========================================"
echo "[1/7] Downloading Metazoa genome metadata"
echo "========================================"

datasets summary genome taxon "${METAZOA_TAXID}" \
    --annotated \
    --assembly-version current \
    --exclude-atypical \
    --exclude-multi-isolate \
    --as-json-lines \
    > "${META}/metazoa.genomes.jsonl"


############################################
# STEP 2
# Convert genome metadata to TSV
############################################

echo
echo "========================================"
echo "[2/7] Converting metadata to TSV"
echo "========================================"

dataformat tsv genome \
    --inputfile "${META}/metazoa.genomes.jsonl" \
    --elide-header \
    --fields \
accession,current-accession,organism-name,organism-tax-id,assminfo-level,assminfo-refseq-category,source_database,assminfo-release-date,assmstats-scaffold-n50,assmstats-contig-n50,annotinfo-release-date,annotinfo-name \
    > "${META}/metazoa.genomes.tsv"

echo -n "Genome records: "
wc -l < "${META}/metazoa.genomes.tsv"


############################################
# STEP 3
# Get all organism TaxIDs
############################################

echo
echo "========================================"
echo "[3/7] Collecting taxonomy IDs"
echo "========================================"

awk -F '\t' '
    NF >= 4 && $4 != "" {
        print $4
    }
' "${META}/metazoa.genomes.tsv" \
    | sort -u \
    > "${META}/organism_taxids.txt"

echo -n "Unique organism TaxIDs: "
wc -l < "${META}/organism_taxids.txt"


############################################
# STEP 4
# Convert organism TaxIDs to species TaxIDs
############################################

echo
echo "========================================"
echo "[4/7] Resolving species-level taxonomy"
echo "========================================"

datasets summary taxonomy taxon \
    --inputfile "${META}/organism_taxids.txt" \
    --as-json-lines \
    > "${META}/taxonomy.jsonl"

dataformat tsv taxonomy \
    --inputfile "${META}/taxonomy.jsonl" \
    --template tax-summary \
    > "${META}/taxonomy.tsv"


############################################
# STEP 5
# Select the best assembly for each species
############################################

echo
echo "========================================"
echo "[5/7] Selecting best genome per species"
echo "========================================"

python3 - \
    "${META}/metazoa.genomes.tsv" \
    "${META}/taxonomy.tsv" \
    "${META}/selected_genomes.tsv" \
    "${META}/selected_accessions.txt" \
    "${META}/unresolved_taxids.txt" <<'PY'

import csv
import sys
from collections import defaultdict
from datetime import datetime


genome_file = sys.argv[1]
taxonomy_file = sys.argv[2]
output_file = sys.argv[3]
accession_file = sys.argv[4]
unresolved_file = sys.argv[5]


# ============================================================
# Helpers
# ============================================================

def clean(value):
    if value is None:
        return ""
    return str(value).strip()


def number(value):
    value = clean(value)

    if not value:
        return 0

    value = value.replace(",", "")

    try:
        return int(float(value))
    except ValueError:
        return 0


def date_score(value):
    """
    Convert YYYY-MM-DD or ISO date into YYYYMMDD integer.
    Invalid/missing dates -> 0.
    """

    value = clean(value)

    if not value:
        return 0

    value = value[:10]

    try:
        dt = datetime.strptime(value, "%Y-%m-%d")
        return dt.year * 10000 + dt.month * 100 + dt.day
    except ValueError:
        return 0


def refseq_category_score(value):
    """
    NCBI RefSeq category:

        reference genome      -> 2
        representative genome -> 1
        none                  -> 0
    """

    value = clean(value).lower()

    if "reference" in value:
        return 2

    if "representative" in value:
        return 1

    return 0


def assembly_level_score(value):

    value = clean(value).lower()

    rank = {
        "complete genome": 4,
        "complete": 4,
        "chromosome": 3,
        "scaffold": 2,
        "contig": 1,
    }

    return rank.get(value, 0)


def source_score(accession, source):

    accession = clean(accession).upper()
    source = clean(source).upper()

    if accession.startswith("GCF_"):
        return 1

    if source == "REFSEQ":
        return 1

    return 0


# ============================================================
# Read taxonomy table
# ============================================================

taxonomy = {}

with open(taxonomy_file, newline="", encoding="utf-8") as fh:

    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:

        taxid = clean(row.get("Taxid"))
        rank = clean(row.get("Rank")).upper()

        species_taxid = clean(row.get("Species taxid"))
        species_name = clean(row.get("Species name"))

        # Normally Species taxid is supplied by NCBI.
        # Fallback for records already at species rank.
        if not species_taxid and rank == "SPECIES":
            species_taxid = taxid
            species_name = clean(row.get("Tax name"))

        taxonomy[taxid] = {
            "species_taxid": species_taxid,
            "species_name": species_name,
            "rank": rank,
        }


# ============================================================
# Read genome table
#
# Column order was fixed by dataformat above:
#
# 0 accession
# 1 current-accession
# 2 organism-name
# 3 organism-tax-id
# 4 assembly-level
# 5 refseq-category
# 6 source_database
# 7 assembly-release-date
# 8 scaffold-N50
# 9 contig-N50
# 10 annotation-release-date
# 11 annotation-name
# ============================================================

species_candidates = defaultdict(list)
unresolved_taxids = set()

with open(genome_file, newline="", encoding="utf-8") as fh:

    reader = csv.reader(fh, delimiter="\t")

    for row in reader:

        if len(row) < 12:
            continue

        (
            accession,
            current_accession,
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
        ) = row[:12]

        accession = clean(accession)
        current_accession = clean(current_accession)
        organism_taxid = clean(organism_taxid)

        if not accession:
            continue

        # ----------------------------------------------------
        # NCBI's "current" query can still expose an older
        # counterpart in a current GCA/GCF pair.
        #
        # Keep only the latest accession in each revision chain.
        # ----------------------------------------------------

        if current_accession and accession != current_accession:
            continue

        tax = taxonomy.get(organism_taxid)

        if not tax:
            unresolved_taxids.add(organism_taxid)
            continue

        species_taxid = tax["species_taxid"]
        species_name = tax["species_name"]

        # We want one genome per actual species.
        # Taxa that cannot be resolved to species are excluded.
        if not species_taxid:
            unresolved_taxids.add(organism_taxid)
            continue

        record = {
            "species_taxid": species_taxid,
            "species_name": species_name,
            "organism_taxid": organism_taxid,
            "organism_name": clean(organism_name),
            "accession": accession,
            "source_database": clean(source_database),
            "refseq_category": clean(refseq_category),
            "assembly_level": clean(assembly_level),
            "scaffold_n50": number(scaffold_n50),
            "contig_n50": number(contig_n50),
            "annotation_release_date": clean(annotation_release_date),
            "assembly_release_date": clean(assembly_release_date),
            "annotation_name": clean(annotation_name),
        }

        # ----------------------------------------------------
        # Ranking
        #
        # Higher tuple = better.
        #
        # Priority:
        #
        # 1 reference / representative genome
        # 2 assembly level
        # 3 RefSeq over GenBank
        # 4 scaffold N50
        # 5 contig N50
        # 6 annotation date
        # 7 assembly date
        #
        # ----------------------------------------------------

        record["score"] = (
            refseq_category_score(refseq_category),
            assembly_level_score(assembly_level),
            source_score(accession, source_database),
            number(scaffold_n50),
            number(contig_n50),
            date_score(annotation_release_date),
            date_score(assembly_release_date),
            accession,
        )

        species_candidates[species_taxid].append(record)


# ============================================================
# Select best genome
# ============================================================

selected = []

for species_taxid, candidates in species_candidates.items():

    candidates.sort(
        key=lambda x: x["score"],
        reverse=True
    )

    best = candidates[0]

    best["candidate_count"] = len(candidates)

    selected.append(best)


# Sort by species name
selected.sort(
    key=lambda x: (
        x["species_name"].lower(),
        x["species_taxid"]
    )
)


# ============================================================
# Write metadata table
# ============================================================

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
    "annotation_release_date",
    "assembly_release_date",
    "annotation_name",
    "candidate_count",
]

with open(output_file, "w", newline="", encoding="utf-8") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t",
        extrasaction="ignore",
    )

    writer.writeheader()

    for record in selected:
        writer.writerow(record)


# ============================================================
# Write accession list
# ============================================================

with open(accession_file, "w", encoding="utf-8") as out:

    for record in selected:
        out.write(record["accession"] + "\n")


# ============================================================
# Unresolved taxonomy
# ============================================================

with open(unresolved_file, "w", encoding="utf-8") as out:

    for taxid in sorted(unresolved_taxids):
        if taxid:
            out.write(taxid + "\n")


print(f"Species selected: {len(selected):,}")
print(f"Unresolved TaxIDs: {len(unresolved_taxids):,}")

PY


echo
echo -n "Selected species/genomes: "
wc -l < "${META}/selected_accessions.txt"


############################################
# STEP 6
# Create dehydrated package
############################################

echo
echo "========================================"
echo "[6/7] Creating protein download package"
echo "========================================"

datasets download genome accession \
    --inputfile "${META}/selected_accessions.txt" \
    --include protein \
    --dehydrated \
    --filename "${PACKAGE}/metazoa_best_proteomes.zip"


############################################
# STEP 7
# Rehydrate protein FASTA files
############################################

echo
echo "========================================"
echo "[7/7] Downloading protein FASTA files"
echo "========================================"

mkdir -p "${PACKAGE}/data"

unzip -q -o \
    "${PACKAGE}/metazoa_best_proteomes.zip" \
    -d "${PACKAGE}/data"

datasets rehydrate \
    --directory "${PACKAGE}/data" \
    --gzip \
    --max-workers "${THREADS}"


############################################
# Build protein file list
############################################

find "${PACKAGE}/data/ncbi_dataset/data" \
    -type f \
    -name "protein.faa.gz" \
    | sort \
    > "${META}/protein_files.list"


echo
echo "========================================"
echo "Finished"
echo "========================================"

echo -n "Selected genomes : "
wc -l < "${META}/selected_accessions.txt"

echo -n "Protein FASTA    : "
wc -l < "${META}/protein_files.list"

echo
echo "Metadata:"
echo "  ${META}/selected_genomes.tsv"

echo
echo "Accessions:"
echo "  ${META}/selected_accessions.txt"

echo
echo "Protein files:"
echo "  ${META}/protein_files.list"
