VISOR HACk -g GRCh38_chromosomes.fa -b snps/HACk.snps.h1.bed snps/HACk.snps.h2.bed -o GRCh38withSNP
VISOR HACk -g GRCh38withSNP/h1.fa -b H1.bed.indel -o H1Indel
VISOR HACk -g GRCh38withSNP/h2.fa -b H2.bed.indel -o H2Indel
VISOR HACk -g H1Indel/h1.fa -b H1.bed -o H1
if [ $? -ne 0 ];then
exit 1
fi
VISOR HACk -g H2Indel/h1.fa -b H2.bed -o H2
if [ $? -ne 0 ];then
exit 1
fi
mkdir $1
cp H1/h1.fa $1/h1.fa
cp H1/h1.fa.fai $1/h1.fa.fai
cp H2/h1.fa $1/h2.fa
cp H2/h1.fa.fai $1/h2.fa.fai
