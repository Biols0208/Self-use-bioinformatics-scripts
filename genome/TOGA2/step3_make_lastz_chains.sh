#!/bin/bash
set +o posix
set -eo pipefail

export JAVA_HOME=/ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/jdk-24.0.1
export PATH="/ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/jdk-24.0.1/bin/:$PATH"
export NXF_SINGULARITY_CACHEDIR=/ldfsqd1/ST_OCEAN/USER/c-lishuo/09_test/GL_RNAseq/SINGULARITY_CACHEDIR

ref_file=$PWD/prepare_data/ref.fa
ref_name=ref
qur_file=$PWD/prepare_data/qur.fa
qur_name=qur
outdir=make_lastz_chains_out

/ldfsqd1/ST_OCEAN/USER/c-lishuo/09_test/GL_RNAseq/nextflow run /ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/make_lastz_chains/main.nf -params-file /ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/make_lastz_chains/params.json -profile singularity \
    --reference_genome ${ref_file} \
    --reference_name ${ref_name} \
    --query_genome ${qur_file} \
    --query_name ${qur_name} \
    --outdir ${outdir} \
    -resume

#rm -rf .nextflow .nextflow.log* work

echo "All done."
