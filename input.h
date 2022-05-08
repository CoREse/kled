#ifndef KLED_INPUT_H
#define KLED_INPUT_H

#include "contig.h"
#include "signature.h"
#include <vector>
#include "kled.h"
#include "htslib/htslib/sam.h"

Contig * getContigs(const char * ReferenceFN, int &NSeq, int WindowSize=100);

struct Stats
{
    float MedianIS,UpIS,BelowIS,Mean,SD;
};

struct Sam
{
	htsFile* SamFile;//the file of BAM/CRAM
	bam_hdr_t* Header;//header for BAM/CRAM file
	hts_idx_t* BamIndex;
	Sam();
	void close();
};

void collectSignatures(Contig &TheContig, std::vector<Signature> *ContigTypeSignatures, SegmentSet & AllPrimarySegments, Arguments & Args, std::vector<Sam>& SamFiles, std::vector<Stats> AllStats, std::vector<int> AllTechs, double* CoverageWindows, const char * DataSource=0);

std::vector<Stats> getAllStats(const char * ReferenceFileName, const std::vector<const char *> & BamFileNames, std::vector<int> AllTechs);

std::vector<int> getAllTechs(Arguments & Args);

std::vector<Sam> initSam(Arguments & Args);
void closeSam(std::vector<Sam> &SamFiles);

#define read_is_unmapped(b) (((b)->core.flag&BAM_FUNMAP) != 0)
#define read_mate_is_unmapped(b) (((b)->core.flag&BAM_FMUNMAP) != 0)
#define read_is_paired(b) (((b)->core.flag&BAM_FPAIRED) != 0)
#define read_is_read1(b) (((b)->core.flag&BAM_FREAD1) != 0)
#define read_is_read2(b) (((b)->core.flag&BAM_FREAD2) != 0)
#define align_is_primary(b) ((((b)->core.flag&BAM_FSECONDARY) == 0) && (((b)->core.flag&BAM_FSUPPLEMENTARY) == 0))
#define align_is_secondary(b) (((b)->core.flag&BAM_FSECONDARY) != 0)
#define align_is_supplementary(b) (((b)->core.flag&BAM_FSUPPLEMENTARY) != 0)
#define read_is_forward(b) (((b)->core.flag&BAM_FREVERSE) == 0)
#define read_mate_is_forward(b) (((b)->core.flag&BAM_FMREVERSE) == 0)
#endif