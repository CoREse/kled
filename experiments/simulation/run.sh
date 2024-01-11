#Seed=`date +"%s"`
Seed=1704683122
SimulatedName="SG"
bash prepareFastas.sh
Ref=GRCh38_chromosomes.fa
cd snps
bash run.sh 2>run.log 1>&2
cd ..
bash genSVBeds.sh $Seed
bash genvcfs.sh
bash genfa.sh $SimulatedName
bash genreads.sh $Ref $SimulatedName $Seed
bash downsamplebam.sh $SimulatedName