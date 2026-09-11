#!/bin/bash
set +o posix
set -eo pipefail

ref_2bit=$PWD/make_lastz_chains_out/ref.2bit
query_2bit=$PWD/make_lastz_chains_out/qur.2bit
chain_file=$PWD/make_lastz_chains_out/ref.qur.allfilled.chain
ref_annotation=$PWD/prepare_data/filter.rename.transcripts.bed
isoform_file=$PWD/prepare_data/filter.toga.isoforms.tsv

SIF=/01_soft/singularity/toga2.sif
OUT=$PWD/final_result_toga
TMP=$PWD/tmp_toga

mkdir -p "${TMP}"
rm -rf "${OUT}"

SINGULARITYENV_TMPDIR="${TMP}" \
singularity exec \
    ${SIF} \
    toga2.py run \
    --ref_2bit ${ref_2bit} \
    --query_2bit ${query_2bit} \
    --chain_file ${chain_file} \
    --ref_annotation ${ref_annotation} \
    --isoform_file ${isoform_file} \
    --no_spliceai \
    --no_u12_file \
    --parallel_strategy local \
    --project_name toga2_anno \
    --output "${OUT}"
