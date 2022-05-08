#ifndef KLED_SIGNATURE_H
#define KLED_SIGNATURE_H

#include <string>
#include <vector>
#include <tuple>

struct Segment
{
    int Begin, End;
    Segment(int Begin,int End);
    bool operator==(const Segment &) const;
    bool operator<(const Segment &) const;
};

class Signature
{
    public:
    int Type;//0: base class, for cigar based signatures, 1: DRPSignatures, 2: ClipSignatures, -1: mark of deleted signature object.
    int Tech;//0: SMRT, 1:NGS
    int SupportedSV;//0: DEL, 1: INS, 2: DUP
    static const char* SVTypeNames[];
    std::string TemplateName;
    std::string InsBases;
    //unsigned long long ReadNum;//indicate the signature is from ReadNum-th read this software reads in this contig, count by primary alignment.
    int Begin, End;//0-based, for insertions, end is just for SVLen calculation and clustering convinience.
    int CN;//default -1, 0 for del, >1 for dup if known
    std::vector<Segment> Segments;
    int Length;
    bool InvLeft, InvRight;
    Signature();
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, const char * InsBases="");
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, Segment Read1, Segment Read2, int Length);//for drp Signatures
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, std::vector<Segment> Segments);//for Split reads signatures
    void setCN(int CN);
    void setInvLeft(bool);
    void setInvRight(bool);
    bool operator<(const Signature &Other) const;
    bool operator==(const Signature & Other) const;
};

int precisionLevel(const Signature &A);
int bestPrecision(const Signature &A,const Signature &B);
int worstPrecision(const Signature &A,const Signature &B);

struct SegmentSet
{
    std::vector<Segment> Segments;
    std::vector<unsigned> MaxEnds;
    std::tuple<int,int> getInvolved(int Start, int End);//Start and End is inclusive, return is [)
    void add(int Begin, int End);
    void sortNStat();
    Segment& operator[](unsigned i);
    SegmentSet();
};

#endif