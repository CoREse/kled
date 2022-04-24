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
    std::string genHeader();
};

class VCFRecord
{
    int SVLen;
    std::string SVType;
    int SS;
    int ST;
    int LS;
    double CV;
    bool Precise;
    std::string InsConsensus;
    int SVTypeI;
    //temp
    // double CS;
    // std::vector<Signature> Cluster;
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

    int getSVTypeI() const;

    VCFRecord(const Contig & TheContig, faidx_t * Ref, std::vector<Signature> & SignatureCluster, SegmentSet & AllPrimarySegments, double* CoverageWindows, double WholeCoverage, Arguments &Args, double* CoverageWindowsSums=NULL, double* CheckPoints=NULL, int CheckPointInterval=0);
    void resolveRef(const Contig & TheContig, faidx_t * Ref, unsigned TypeCount, Arguments & Args);
    std::string genotype(const Contig & TheContig, SegmentSet & AllPrimarySegments, double * CoverageWindows, double *CoverageWindowsSums, double* Checkpoints, int CheckPointInterval, Arguments & Args);
    operator std::string() const;
    bool operator<(const VCFRecord& Other) const;
};

void addKledEntries(VCFHeader & Header);

double getAverageCoverage(int Begin, int End, double * CoverageWindows, Arguments & Args, double* CoverageWindowsSums=NULL, double* CheckPoints=NULL, int CheckPointInterval=0);

#endif