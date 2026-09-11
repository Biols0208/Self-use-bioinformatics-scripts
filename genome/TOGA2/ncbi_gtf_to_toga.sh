#!/usr/bin/env bash
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: bash $0 genomic.gtf prefix"
    exit 1
fi

gtf=$1
prefix=$2

raw_bed="${prefix}.raw.bed"
tx_meta="${prefix}.transcript_meta.tsv"
gene_locus="${prefix}.gene_locus.tsv"
gene_locus_sorted="${prefix}.gene_locus.sorted.tsv"
gene_unique="${prefix}.gene_unique_name.tsv"
tx_map="${prefix}.transcript_toga_name.tsv"

out_bed="${prefix}.toga.transcripts.bed"
out_iso="${prefix}.toga.isoforms.tsv"

gtfToGenePred="/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/env.sh /ldfsqd1/ /ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/singularity/ubuntu.24.04.sif /ldfsqd1/ST_OCEAN/USER/c-lishuo/09_test/AH_toga/gtfToGenePred"
genePredToBed="/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/env.sh /ldfsqd1/ /ldfsqd1/ST_OCEAN/USER/c-lishuo/01_soft/singularity/ubuntu.24.04.sif /ldfsqd1/ST_OCEAN/USER/c-lishuo/09_test/AH_toga/genePredToBed"

echo "========================================"
echo "NCBI GTF -> TOGA input"
echo "GTF    : ${gtf}"
echo "Prefix : ${prefix}"
echo "========================================"

# ============================================================
# 1. GTF -> genePred -> BED12
# ============================================================

echo "[1/7] Convert GTF to BED12..."

$gtfToGenePred \
    -genePredExt \
    -geneNameAsName2 \
    -ignoreGroupsWithoutExons \
    "${gtf}" stdout | \
$genePredToBed stdin "${raw_bed}"


# ============================================================
# 2. Extract transcript_id / gene_id / gene symbol
#
# output:
# transcript_id    gene_id    gene_symbol
# ============================================================

echo "[2/7] Extract transcript metadata..."

awk -F'\t' '
BEGIN{OFS="\t"}

function get_attr(s,key,   n,a,i,x,prefix){
    n=split(s,a,";");
    prefix=key " ";

    for(i=1;i<=n;i++){
        x=a[i];
        gsub(/^[ \t]+|[ \t]+$/, "", x);

        if(index(x,prefix)==1){
            sub("^[^ \t]+[ \t]+", "", x);
            gsub(/^"|"$/, "", x);
            return x;
        }
    }

    return "";
}

$3=="transcript" {
    tid=get_attr($9,"transcript_id");
    gid=get_attr($9,"gene_id");

    gene=get_attr($9,"gene");

    if(gene=="")
        gene=get_attr($9,"gene_name");

    if(gene=="")
        gene=gid;

    if(tid!="" && gid!="")
        print tid,gid,gene;
}
' "${gtf}" > "${tx_meta}"


# ============================================================
# 3. Extract genomic locus for each gene_id
#
# output:
# gene_id    gene_symbol    chromosome    start    end
#
# Use all GTF records so this also works if explicit
# "gene" features are absent.
# ============================================================

echo "[3/7] Determine gene loci..."

awk -F'\t' '
BEGIN{OFS="\t"}

function get_attr(s,key,   n,a,i,x,prefix){
    n=split(s,a,";");
    prefix=key " ";

    for(i=1;i<=n;i++){
        x=a[i];
        gsub(/^[ \t]+|[ \t]+$/, "", x);

        if(index(x,prefix)==1){
            sub("^[^ \t]+[ \t]+", "", x);
            gsub(/^"|"$/, "", x);
            return x;
        }
    }

    return "";
}

{
    gid=get_attr($9,"gene_id");

    if(gid=="")
        next;

    gene=get_attr($9,"gene");

    if(gene=="")
        gene=get_attr($9,"gene_name");

    if(gene=="")
        gene=gid;

    if(!(gid in seen)){
        seen[gid]=1;
        chr[gid]=$1;
        start[gid]=$4;
        end[gid]=$5;
        symbol[gid]=gene;
    }
    else{
        if($4 < start[gid])
            start[gid]=$4;

        if($5 > end[gid])
            end[gid]=$5;

        if(symbol[gid]==gid && gene!=gid)
            symbol[gid]=gene;
    }
}

END{
    for(gid in seen)
        print gid,symbol[gid],chr[gid],start[gid],end[gid];
}
' "${gtf}" > "${gene_locus}"


# ============================================================
# 4. Sort loci:
# gene_symbol -> chromosome -> genomic position
#
# This makes xxx-1 / xxx-2 assignment reproducible.
# ============================================================

echo "[4/7] Resolve duplicated gene symbols..."

LC_ALL=C sort \
    -t $'\t' \
    -k2,2 \
    -k3,3V \
    -k4,4n \
    -k5,5n \
    "${gene_locus}" \
    > "${gene_locus_sorted}"


# ============================================================
# 5. Make unique gene display names
#
# One locus:
#   dmrt1 -> dmrt1
#
# Multiple independent gene_ids with same symbol:
#   foxl2 geneA -> foxl2-1
#   foxl2 geneB -> foxl2-2
#
# output:
# gene_id    unique_gene_name
# ============================================================

awk -F'\t' '
BEGIN{OFS="\t"}

NR==FNR{
    count[$2]++;
    next;
}

{
    gid=$1;
    gene=$2;

    if(count[gene]==1){
        newname=gene;
    }
    else{
        n[gene]++;
        newname=gene "-" n[gene];
    }

    print gid,newname;
}
' "${gene_locus_sorted}" "${gene_locus_sorted}" \
    > "${gene_unique}"


# ============================================================
# 6. Transcript -> TOGA name
#
# output:
# transcript_id    transcript_id#gene_name    gene_id
# ============================================================

awk -F'\t' '
BEGIN{OFS="\t"}

NR==FNR{
    name[$1]=$2;
    next;
}

{
    tid=$1;
    gid=$2;
    gene=$3;

    if(gid in name)
        gene=name[gid];

    print tid,tid "#" gene,gid;
}
' "${gene_unique}" "${tx_meta}" \
    > "${tx_map}"


# ============================================================
# 7a. Replace BED column 4 and retain coding transcripts only
#
# BED:
# $7 = thickStart
# $8 = thickEnd
#
# thickStart < thickEnd means CDS exists.
# ============================================================

echo "[5/7] Generate TOGA BED12..."

awk -F'\t' '
BEGIN{OFS="\t"}

NR==FNR{
    toga[$1]=$2;
    next;
}

$7 < $8 {
    if($4 in toga){
        $4=toga[$4];
        print;
    }
}
' "${tx_map}" "${raw_bed}" \
    > "${out_bed}"


# ============================================================
# 7b. Generate isoforms file
#
# gene_id    transcript_id#gene_name
#
# Order follows BED exactly.
# ============================================================

echo "[6/7] Generate isoforms file..."

awk -F'\t' '
BEGIN{OFS="\t"}

NR==FNR{
    gid[$2]=$3;
    next;
}

{
    if($4 in gid)
        print gid[$4],$4;
}
' "${tx_map}" "${out_bed}" \
    > "${out_iso}"


# ============================================================
# QC
# ============================================================

echo "[7/7] QC..."

n_bed=$(wc -l < "${out_bed}")
n_iso=$(wc -l < "${out_iso}")

echo
echo "========================================"
echo "Finished"
echo "========================================"
echo "BED12    : ${out_bed}"
echo "Isoforms : ${out_iso}"
echo
echo "BED transcripts     : ${n_bed}"
echo "Isoform transcripts : ${n_iso}"
echo

if [ "${n_bed}" -ne "${n_iso}" ]; then
    echo "WARNING: BED and isoforms line counts differ!"
fi

echo
echo "Duplicated gene symbols:"
awk -F'\t' '
{
    n[$2]++;
}
END{
    for(i in n)
        if(n[i]>1)
            print i,n[i];
}
' "${gene_locus}" | sort -k2,2nr | head -20

echo
echo "BED preview:"
head -3 "${out_bed}"

echo
echo "Isoforms preview:"
head -3 "${out_iso}"

echo
echo "Check BED columns:"
awk 'NF!=12{print "ERROR line",NR,"NF="NF}' "${out_bed}" | head

echo
echo "Check duplicated transcript names:"
cut -f4 "${out_bed}" | sort | uniq -d | head

echo
echo "Done."
