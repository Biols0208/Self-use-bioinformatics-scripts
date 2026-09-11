
export PATH="/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/:$PATH"

ref_genome=Oncorhynchus_mykiss.softmak.fa
qur_genome=Brachymystax_lenok_tsinlingensis.softmak.fa

[ -d prepare_data ] || mkdir prepare_data

/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/iTools Fatools stat -InPut ${ref_genome} -OutPut prepare_data/ref.stat
/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/iTools Fatools stat -InPut ${qur_genome} -OutPut prepare_data/qur.stat
rm *.chrlist

sed '1d' prepare_data/ref.stat | awk '{print $1"\tchr"NR}' > prepare_data/ref.rename.match.tsv
sed '1d' prepare_data/qur.stat | awk '{print $1"\tCHR"NR}' > prepare_data/qur.rename.match.tsv

## rename ref genome
/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/rename.fa.py --ids prepare_data/ref.rename.match.tsv ${ref_genome} | awk '{print $1}' | seqkit seq -w 100 > prepare_data/ref.fa
/ldfsqd1/ST_OCEAN/USER/c-lishuo/00_tools/rename.fa.py --ids prepare_data/qur.rename.match.tsv ${qur_genome} | awk '{print $1}' | seqkit seq -w 100 > prepare_data/qur.fa

## rename ref bed
awk 'NR==FNR{a[$1]=$2}NR!=FNR{print a[$1]"\t"$0}' prepare_data/ref.rename.match.tsv prepare_data/deal.toga.transcripts.bed | cut -f 1,3- > prepare_data/final.rename.transcripts.bed

## filter too long pep
grep -w -v -f prepare_data/too_long.id prepare_data/final.rename.transcripts.bed > prepare_data/filter.rename.transcripts.bed
grep -w -v -f prepare_data/too_long.id prepare_data/deal.toga.isoforms.tsv > prepare_data/filter.toga.isoforms.tsv
