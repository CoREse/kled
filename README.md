# Kled: An ultra-fast and sensitive structural variant detection tool for long-read sequencing data
## Introduction
Kled is designed to call SVs nicely and quickly using long-read sequencing data. It takes mapped reads file (bam) as input and reports SVs to the stdout in the VCF file format. Kled can yield precise and comprehensive SV detection results within minutes and can run on any modern computer without needing of any field knowledge of the user to perform the SV detection.

## Compiling
Kled uses cmake build tools to build the project.

Make sure you have the following dependencies and cmake tools (>=3.15), g++ (gxx) (>=15.1):
- zlib >=1.3
- bzip2 >=1.0
- liblzma-devel >=5.8
- libcurl >=8.0
- openssl >=3.5
- libdeflate >=1.24
- libboost-devel >=1.84
- gmp >=6.3

To build the project, run:
```
git clone https://github.com/CoREse/kled
cd kled
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/your/path ..
cmake --build . -j 16
cmake --install .
```
And you will have the kled built and kled and HapKled installed,
### Conda
Kled is now on bioconda! To get kled, simply:
```
conda install kled -c bioconda
```
And that's all! But before that, be sure you will install kled in the right environment, to create a dedicated environment for kled:
```
conda create -n kled kled -c bioconda
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

## HapKled usage
HapKled is a script that helps you handling kled with haplotype-tagged input.

Before running HapKled you should first compile (see compiling the haplotype-aware kled) the haplotype-aware kled and put the path of it to the environment variable HapAwareKled.
```
export HapAwareKled=/path/to/hap-aware-kled
```
And you also need an installed Clair3 and Whatshap, and export the path of the Clair3 models to environment variable Clair3ModelPath.
```
export Clair3ModelPath=/path/to/bin/models
```
HapKled need a reference file (fasta) and at least one bam (sam/bam/cram) file that stores the mapped reads to call SVs, and output a VCF file to the standard output.
```
HapKled -R Refernce.fa Sample.bam > SVs.vcf
```

For the description of all parameters:
```
HapKled --help
```

## Citation
This work is published on [*Briefings in Bioinformatics*](https://academic.oup.com/bib/article/25/2/bbae049/7611936), doi:10.1093/bib/bbae049, please visit the site for citations.

The HapKled is published on [*Frontiers in Genetics*](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2024.1435087/full), doi:10.3389/fgene.2024.1435087, please visit the site for citations.