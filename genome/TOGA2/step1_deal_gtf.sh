
input_gtf=GCF_013265735.2_USDA_OmykA_1.1_genomic.gtf

## save file
[ -d prepare_data ] || mkdir prepare_data

ncbi_gtf_to_toga.sh ${input_gtf} deal

awk -F'\t' '$3=="CDS"{split($9,a," "); ID=a[4]; L[ID]+=$5-$4+1} END{for(i in L) print i"\t"L[i]/3}' ${input_gtf} | sed "s/\"//g;s/;//g" | sort -k 2nr,2 > prepare_data/transcript.len
awk '$2>10000' prepare_data/transcript.len | awk '{print $1}' > prepare_data/too_long.id

mv deal* prepare_data
