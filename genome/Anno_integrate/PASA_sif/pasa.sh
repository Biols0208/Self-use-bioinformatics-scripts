#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob
export TERM=${TERM:-dumb}

SINGULARITY_BIN=${SINGULARITY_BIN:-/usr/bin/singularity}
SIF=${SIF:-/soft/pasa.sif}
WORK_DIR=${WORK_DIR:-$(pwd -P)}

container_pasahome=/opt/PASApipeline-turbo
container_perl5lib=/opt/PASApipeline-turbo/PerlLib
container_blastmat=/opt/ncbi-blast-2.2.26/data
container_path=/opt/PASApipeline-turbo:/opt/PASApipeline-turbo/scripts:/opt/PASApipeline-turbo/misc_utilities:/opt/PASApipeline-turbo/bin:/opt/PASApipeline-turbo/pasa-plugins/transdecoder:/usr/local/bin:/usr/local/sbin:/usr/bin:/usr/sbin:/bin:/sbin

genome=purged.FINAL.from_pure.fasta
transcripts=all.isoforms.99.fa
transcripts_clean=step1_clean/${transcripts}.clean
annotation_source=Bathypathes_sp._nov._v1.HiC.gff
annotation_gff3=${annotation_source}.gff3

# Optional replacement for the old change_gff_format.pl. Leave empty when annotation_gff3 already exists or annotation_source is already valid GFF3.
# When used, the converter must be accessible under WORK_DIR and accept: perl converter.pl input.gff output.gff3
gff_converter=/soft/change_gff_format.pl

config_align=alignAssembly.config
config_annot=annotCompare.config

cpu=20
max_intron_length=2000000
aligner=minimap2                 # supported: gmap, blat, minimap2

# These two clustering modes are mutually exclusive in PASA.
cluster_mode=stringent           # stringent or gene
stringent_alignment_overlap=30
gene_overlap=50

seqclean_cpus=15
seqclean_min_length=100
univec=/opt/pasa-data/UniVec
run_annotation_round2=true
# --------------------------------------------------------------

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

on_error() {
    local exit_code=$?
    printf 'ERROR: command failed at line %s (exit=%s)\n' \
        "${BASH_LINENO[0]}" "$exit_code" >&2
    exit "$exit_code"
}
trap on_error ERR

[[ -x "$SINGULARITY_BIN" ]] || die "Singularity executable not found: $SINGULARITY_BIN"
[[ -f "$SIF" ]] || die "SIF image not found: $SIF"
[[ -d "$WORK_DIR" ]] || die "WORK_DIR not found: $WORK_DIR"

SIF=$(realpath -e "$SIF")
WORK_DIR=$(realpath -e "$WORK_DIR")
cd "$WORK_DIR"

for required_file in "$genome" "$transcripts" "$config_align" "$config_annot"; do
    [[ -s "$required_file" ]] || die "Required input is missing or empty: $WORK_DIR/$required_file"
done

[[ "$cpu" =~ ^[1-9][0-9]*$ ]] || die "cpu must be a positive integer"
[[ "$max_intron_length" =~ ^[1-9][0-9]*$ ]] || die "max_intron_length must be a positive integer"
[[ "$seqclean_cpus" =~ ^[1-9][0-9]*$ ]] || die "seqclean_cpus must be a positive integer"
((seqclean_cpus <= 50)) || die "seqclean_cpus cannot exceed the seqclean limit of 50"
[[ "$seqclean_min_length" =~ ^[1-9][0-9]*$ ]] || die "seqclean_min_length must be a positive integer"

case "$aligner" in
    gmap|blat|minimap2) ;;
    pblat)
        die "pblat is installed in the image but is not accepted by PASA --ALIGNERS; use blat, gmap or minimap2"
        ;;
    *)
        die "Unsupported aligner: $aligner (choose gmap, blat or minimap2)"
        ;;
esac

case "$cluster_mode" in
    stringent|gene) ;;
    *) die "cluster_mode must be stringent or gene" ;;
esac

run_sif() {
    "$SINGULARITY_BIN" exec \
        --cleanenv \
        --bind "/ldfsqd1/:/ldfsqd1/" \
        --pwd "$PWD" \
        "$SIF" \
        /usr/bin/env \
        "PASAHOME=$container_pasahome" \
        "PASAPIPELINE=$container_pasahome" \
        "PERL5LIB=$container_perl5lib" \
        "BLASTMAT=$container_blastmat" \
        "PATH=$container_path" \
        "LANG=C.UTF-8" \
        "LC_ALL=C.UTF-8" \
        "TERM=$TERM" \
        "USER=${USER:-pasa}" \
        "$@"
}

read_database_value() {
    local config_file=$1
    awk -F= '
        /^[[:space:]]*#/ { next }
        $1 ~ /^[[:space:]]*DATABASE[[:space:]]*$/ {
            value=$0
            sub(/^[^=]*=/, "", value)
            sub(/^[[:space:]]+/, "", value)
            sub(/[[:space:]]+$/, "", value)
            gsub(/^\047|\047$/, "", value)
            gsub(/^\042|\042$/, "", value)
            print value
            exit
        }
    ' "$config_file"
}

looks_like_gff3() {
    local annotation_file=$1
    awk -F '\t' '
        /^##gff-version[[:space:]]+3/ { valid=1; exit }
        /^#/ { next }
        NF >= 9 {
            if ($9 ~ /(^|;)ID=/ || $9 ~ /(^|;)Parent=/) valid=1
            exit
        }
        END { exit(valid ? 0 : 1) }
    ' "$annotation_file"
}

newest_update_file() {
    local files=( *gene_structures_post_PASA_updates.*.gff3 )
    local latest
    local candidate

    ((${#files[@]} > 0)) || return 1
    latest=${files[0]}
    for candidate in "${files[@]:1}"; do
        [[ "$candidate" -nt "$latest" ]] && latest=$candidate
    done
    printf '%s\n' "$latest"
}

printf '[check] validating image and embedded dependencies\n'
run_sif bash -c '
    set -eu
    test "${PASAHOME:-}" = /opt/PASApipeline-turbo
    test -s /opt/pasa-data/UniVec
    test -s /opt/pasa-data/UniVec.nhr
    test -s /opt/pasa-data/UniVec.nin
    test -s /opt/pasa-data/UniVec.nsq
    for program in Launch_PASA_pipeline.pl seqclean psx seqclean.psx \
                   trimpoly mdust cdbfasta cdbyank \
                   blastall megablast formatdb \
                   samtools TransDecoder.LongOrfs
    do
        command -v "$program" >/dev/null || {
            printf "ERROR: required image command is missing: %s\n" "$program" >&2
            exit 1
        }
    done
'

# ------------------------------------------------------------
# Prepare a PASA-compatible GFF3
# ------------------------------------------------------------
if [[ -s "$annotation_gff3" ]]; then
    looks_like_gff3 "$annotation_gff3" || \
        die "Existing annotation_gff3 does not look like GFF3: $annotation_gff3"
elif looks_like_gff3 "$annotation_source"; then
    ln -sfn "$(realpath -e "$annotation_source")" "$annotation_gff3"
elif [[ -n "$gff_converter" ]]; then
    [[ -s "$gff_converter" ]] || die "GFF converter not found: $gff_converter"
    run_sif perl "$gff_converter" "$annotation_source" "$annotation_gff3"
    [[ -s "$annotation_gff3" ]] || die "GFF converter did not create: $annotation_gff3"
    looks_like_gff3 "$annotation_gff3" || \
        die "Converted annotation does not look like GFF3: $annotation_gff3"
else
    die "Annotation is not recognizable as GFF3. Provide annotation_gff3 or set gff_converter to an accessible replacement for change_gff_format.pl"
fi

# ------------------------------------------------------------
# Validate that both PASA configs use one database
# ------------------------------------------------------------
database_align=$(read_database_value "$config_align")
database_annot=$(read_database_value "$config_annot")

[[ -n "$database_align" ]] || die "DATABASE is missing from $config_align"
[[ -n "$database_annot" ]] || die "DATABASE is missing from $config_annot"
[[ "$database_align" == "$database_annot" ]] || \
    die "DATABASE differs between configs: [$database_align] versus [$database_annot]"

if [[ "$database_align" == /* && "$database_align" != "$WORK_DIR"/* ]]; then
    die "SQLite DATABASE must be under WORK_DIR so it is visible in the container: $database_align"
fi

database_name=$(basename -- "$database_align")

# ------------------------------------------------------------
# Step 1: clean transcript sequences
# ------------------------------------------------------------
printf '[step 1] cleaning transcript sequences\n'
step1_dir=$WORK_DIR/step1_clean
transcript_name=$(basename -- "$transcripts")
mkdir -p "$step1_dir"
ln -sfn "$(realpath -e "$transcripts")" "$step1_dir/$transcript_name"

if [[ ! -s "$step1_dir/${transcript_name}.clean" ]]; then
    cd "$step1_dir"
    if ! run_sif seqclean "$transcript_name" \
            -c "$seqclean_cpus" \
            -l "$seqclean_min_length" \
            -v "$univec"; then
        printf 'ERROR: seqclean failed. Diagnostic log tails follow.\n' >&2
        for diagnostic_log in \
            "err_seqcl_${transcript_name}.log" \
            cleaning_*/err_log
        do
            if [[ -s "$diagnostic_log" ]]; then
                printf '\n===== %s =====\n' "$step1_dir/$diagnostic_log" >&2
                tail -n 100 "$diagnostic_log" >&2
            fi
        done
        die "seqclean failed in $step1_dir"
    fi
    cd "$WORK_DIR"
else
    printf '[step 1] existing clean FASTA retained: %s\n' \
        "$step1_dir/${transcript_name}.clean"
fi

[[ -s "$step1_dir/${transcript_name}.clean" ]] || \
    die "seqclean output is missing: $step1_dir/${transcript_name}.clean"
# Turbo can use the genome FASTA index directly. Indexing only the original
# transcript FASTA, as in the old script, did not provide this optimization.
if [[ ! -s "${genome}.fai" ]]; then
    run_sif samtools faidx "$genome"
fi
ln -s step1_clean/${transcripts}.clean ./
ln -s step1_clean/${transcripts}.cln ./

# ------------------------------------------------------------
# Step 2A: alignment and transcript assembly only
# ------------------------------------------------------------
printf '[step 2A] alignment and transcript assembly\n'
alignment_args=(
    Launch_PASA_pipeline.pl
    -c "$config_align"
    -C -r -R
    -g "$genome"
    -t "$transcripts_clean"
    -T
    -u "$transcripts"
    --ALIGNERS "$aligner"
    --CPU "$cpu"
    --MAX_INTRON_LENGTH "$max_intron_length"
    --TRANSDECODER
)

if [[ "$cluster_mode" == stringent ]]; then
    alignment_args+=(--stringent_alignment_overlap "$stringent_alignment_overlap")
else
    alignment_args+=(--gene_overlap "$gene_overlap" --annots "$annotation_gff3")
fi

run_sif "${alignment_args[@]}"

# ------------------------------------------------------------
# Step 2B: comprehensive transcriptome
# ------------------------------------------------------------
printf '[step 2B] building comprehensive transcriptome\n'
run_sif build_comprehensive_transcriptome.dbi \
    -c "$config_align" \
    -t "$transcripts_clean" \
    --min_per_ID 95 \
    --min_per_aligned 30

# ------------------------------------------------------------
# Step 2C: annotation comparison, round 1
# ------------------------------------------------------------
printf '[step 2C] annotation comparison round 1\n'
run_sif Launch_PASA_pipeline.pl \
    -c "$config_annot" \
    -g "$genome" \
    -t "$transcripts_clean" \
    -A -L \
    --annots "$annotation_gff3" \
    --CPU "$cpu"

round1_update=$(newest_update_file) || \
    die "Round 1 did not produce a gene_structures_post_PASA_updates GFF3"
printf '[step 2C] round 1 update: %s\n' "$round1_update"

# ------------------------------------------------------------
# Step 2D: optional annotation comparison, round 2
# ------------------------------------------------------------
if [[ "$run_annotation_round2" == true ]]; then
    printf '[step 2D] annotation comparison round 2\n'
    run_sif Launch_PASA_pipeline.pl \
        -c "$config_annot" \
        -g "$genome" \
        -t "$transcripts_clean" \
        -A -L \
        --annots "$round1_update" \
        --CPU "$cpu"
fi

final_update=$(newest_update_file) || \
    die "No final PASA annotation update GFF3 was found"

# ------------------------------------------------------------
# Step 3: alternative-splicing analysis (must be separate from -A)
# ------------------------------------------------------------
printf '[step 3] alternative-splicing analysis\n'
run_sif Launch_PASA_pipeline.pl \
    -c "$config_annot" \
    -g "$genome" \
    -t "$transcripts_clean" \
    --CPU "$cpu" \
    --ALT_SPLICE

# ------------------------------------------------------------
# Step 4: derive a training set from PASA assemblies
# ------------------------------------------------------------
printf '[step 4] generating PASA training set\n'
assemblies_fasta=${database_name}.assemblies.fasta
pasa_assemblies_gff3=${database_name}.pasa_assemblies.gff3

[[ -s "$assemblies_fasta" ]] || die "Missing PASA assembly FASTA: $assemblies_fasta"
[[ -s "$pasa_assemblies_gff3" ]] || die "Missing PASA assembly GFF3: $pasa_assemblies_gff3"

run_sif pasa_asmbls_to_training_set.dbi \
    --pasa_transcripts_fasta "$assemblies_fasta" \
    --pasa_transcripts_gff3 "$pasa_assemblies_gff3"

printf '\nPASA workflow completed.\n'
printf 'Final updated annotation: %s/%s\n' "$WORK_DIR" "$final_update"
printf 'PASA database: %s\n' "$database_align"
