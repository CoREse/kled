Ref=$1
Sample=$2
Seed=$3
gunzip -c HG002_ONT_UL.ldist.gz > HG002_ONT_UL.ldist
for h in "h1" "h2"; do
    mkdir $Sample/$h.fa.chrs
    Chrs=`cut -f1 $Sample/$h.fa.fai`
    IFS=$'\n'
    for Chr in $Chrs; do
        samtools faidx $Sample/$h.fa $Chr > $Sample/$h.fa.chrs/$Chr.fa &&\
        samtools faidx $Sample/$h.fa.chrs/$Chr.fa &&\
        lrsim $Sample/$h.fa.chrs/$Chr.fa -e 0.2 -d $Depth -s $Seed -f HG002_ONT_UL.ldist > ${Sample}_bam_30x/$h.$Chr.fastq &
        ((Count++))
        if [ $((Count % 34)) -eq 0 ]; then
            wait
        fi
    done
done
wait
minimap2 -a -L -z 600,200 -x map-ont -t 35 -R '@RG\tID:nanopore\tSM:bulk' $Ref ${Sample}_bam_30x/*.fastq | samtools sort -@ 35 --write-index -o ${Sample}_bam_30x/sim.srt.bam