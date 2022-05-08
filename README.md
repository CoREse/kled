# Kled: a SV caller
## Compiling
Dependencies: openmp and dependencies of htslib (-lz -lm -lbz2 -llzma -lcurl -lpthread -lcrypto -ldeflate)

You also need a compiler that supports C++ 17 standard.

To get the kled binary:
```
#to get htslib
git submodule update --init --recursive

#prepare the htslib
make htsconf

#make kled
make
```

If you want to install kled:
```
make install
#or if you want to install to a place other than /usr/local:
PREFIX=PATH_YOU_SELECT make install
```
## Usage
Kled need a reference file (fasta) and at least one bam (sam/bam/cram) file that stores the mapped reads to call SVs, and output a VCF file to the standard output.
```
kled -R Refernce.fa Sample.bam > SVs.vcf
```
The default parameters are tuned for ONT data, if your inputs are CLR or CCS data, consider add --CLR or --CCS option to get a better result:
```
kled -R Reference.fa --CCS CCS.bam > SVs.vcf
kled -R Reference.fa --CLR CLR.bam > SVs.vcf
```