#ifndef KLED_REPORT_H
#define KLED_REPORT_H

#include <string>
#include <vector>
#include "signature.h"
#include "htslib/htslib/faidx.h"
#include "contig.h"
#include <map>
#include "kled.h"
#include "input.h"
#include "clustering.h"

struct HeaderEntry
{
    std::string Class;//INFO,ALT,FORMAT
    std::string ID;
    std::string Type;
    std::string Description;
    std::string Number;
    std::string Source;
    std::string Version;
    HeaderEntry(std::string Class, std::string ID, std::string Description, std::string Number="", std::string Type="", std::string Source="", std::string Version="");
    operator std::string() const;
};

class VCFHeader
{
    std::string FileFormat;
    std::string FileDate;
    std::string Reference;
    std::vector<HeaderEntry> HeaderEntries;
    std::vector<std::string> SampleNames;
    std::vector<Contig> Contigs;
    public:
    VCFHeader(const char * Reference);
    void addHeaderEntry(const HeaderEntry & Entry);
    void addSample(const char * SampleName);
    void addContig(const Contig & TheContig);
    std::string genHeader(const Arguments & Args);
};

class VCFRecord
{
    int SVLen;
    std::string SVType;
    int CN;
    int SS;
    int ST;
    int SS2;
    int ST2;
    int LS;
    double CV;
    double CR;
    int MinLength;
    int MaxLength;
    int MediumLength;
    double PSD;

    bool Precise;
    std::string InsConsensus;
    int SVTypeI;
    DEBUG_CODE(std::string MergeStrings;)
    public:
    std::string CHROM;
    int Pos;//0-based reference Pos of the variant, for insertion is the pos after the insertion, otherwise is the 1st base of the variant
    std::string ID;
    std::string REF;
    std::string ALT;
    std::string QUAL;
    std::string FILTER;
    std::string INFO;
    std::map<std::string,std::string> Sample;

    bool Keep;//keep this record
    int HPCounts[3];
    int ConcurrentGT;//Is HapGT concurrent with GT, 0: not detected, 1: concurrent, 2: not concurrent

    int getSVTypeI() const;
    int getPSD() const;
    int getSVLen() const;
    int getST() const;

    VCFRecord();
    VCFRecord(const Contig & TheContig, std::vector<Signature> & SignatureCluster, ClusterCore &Core, SegmentSet & AllPrimarySegments, float* CoverageWindows, double WholeCoverage, Arguments &Args, float* CoverageWindowsSums=NULL, float* CheckPoints=NULL, int CheckPointInterval=0);
    void resolveHPRecord(int* HPCounts, const Contig & TheContig, std::vector<Signature> & SignatureCluster, ClusterCore &Core, SegmentSet & AllPrimarySegments, float* CoverageWindows, double WholeCoverage, Arguments &Args, float* CoverageWindowsSums=NULL, float* CheckPoints=NULL, int CheckPointInterval=0);
    void resolveRef(const Contig & TheContig, faidx_t * Ref, unsigned TypeCount, double CC, Arguments & Args);
    std::string genotype(const Contig & TheContig, SegmentSet & AllPrimarySegments, float * CoverageWindows, float *CoverageWindowsSums, float* Checkpoints, int CheckPointInterval, Arguments & Args);
    void hapGT(unsigned SmallHapCount, unsigned BigHapCount);
    operator std::string() const;
    bool operator<(const VCFRecord& Other) const;
    void calcM3L(std::vector<Signature> & SignatureCluster, bool ExcludeHP0=false);
};

void addKledEntries(VCFHeader & Header);

double getAverageCoverage(int Begin, int End, float * CoverageWindows, Arguments & Args, float* CoverageWindowsSums=NULL, float* CheckPoints=NULL, int CheckPointInterval=0);

void HPClustersDistinction(std::vector<Signature> &Cluster, std::vector<std::vector<Signature>> &HPClusters, Arguments& Args);

#endif