# Kled: An ultra-fast and sensitive structural variant detection tool for long-read sequencing data
## Introduction
Kled is designed to call SVs nicely and quickly using long-read sequencing data. It takes mapped reads file (bam) as input and reports SVs to the stdout in the VCF file format. Kled can yield precise and comprehensive SV detection results within minutes and can run on any modern computer without needing of any field knowledge of the user to perform the SV detection.

## Project structure
This project contains the source of the kled, which resides in the root path of the project, and the scripts used in the experiments, which reside in the experiments folder.

## Compiling
Dependencies: openmp and dependencies of htslib (-lz -lm -lbz2 -llzma -lcurl -lpthread -lcrypto -ldeflate)

Here are instructions to get some dependencies from source if systemwide installation is not available:
- zlib: download from https://github.com/madler/zlib/releases, unzip and cd and ./configure && make && make install prefix=/path/to/your/local
- libbz2: download from https://sourceforge.net/projects/bzip2/, unzip and cd and make && make install PREFIX=/path/to/your/local
And export CPATH to make compiler know where the lib is: export CPATH=/path/to/your/local/include:$CPATH

You also need a compiler that supports C++ 17 standard.

To get the kled binary:
```
#to get htslib
git submodule update --init --recursive

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

For the description of all parameters:
```
kled --help
```

## Citation
This work is published in *Briefings in Bioinformatics* (https://academic.oup.com/bib/article/25/2/bbae049/7611936), please visit the site for citations.