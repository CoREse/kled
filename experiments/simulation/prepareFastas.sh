wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > GRCh38_chromosomes.fa
samtools faidx GRCh38_chromosomes.fa
bash faitobed.sh GRCh38_chromosomes.fa.fai > GRCh38_chromosomes.fa.bed