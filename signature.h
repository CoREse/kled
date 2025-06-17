#ifndef KLED_SIGNATURE_H
#define KLED_SIGNATURE_H

#include <string>
#include <vector>
#include <tuple>
#include "kled.h"

#ifdef DEBUG
#include <fstream>

// include headers that implement a archive in simple text format
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#endif
struct Segment
{
    int Begin, End;
    std::string InsBases;
    Segment(int Begin=0,int End=0, const char * CInsBases="");
    bool operator==(const Segment &) const;
    bool operator<(const Segment &) const;
    operator std::string();
    #ifdef DEBUG
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & Begin;
        ar & End;
    }
    #endif
};

class Signature : public Segment
{
    public:
    int Type;//0: base class, for cigar based signatures, 1: DRPSignatures, 2: ClipSignatures, -1: mark of deleted signature object.
    int Tech;//0: SMRT, 1:NGS
    int SupportedSV;//0: DEL, 1: INS, 2: DUP
    static const char* SVTypeNames[];
    std::string TemplateName;
    // std::string InsBases;
    //unsigned long long ReadNum;//indicate the signature is from ReadNum-th read this software reads in this contig, count by primary alignment.
    // int Begin, End;//0-based, for insertions, end is just for SVLen calculation and clustering convinience.
    int CN;//default -1, 0 for del, >1 for dup if known
    bool Covered;//Dup if covered.
    std::vector<Segment> Segments;
    int Length;
    bool InvLeft, InvRight;
    double Quality;//Signature quality, the higher the better
    unsigned ID;
    bool Artificial;//Merged merging template
    int HP;//HP tag
    DEBUG_CODE(std::string MergeString;)
    Signature();
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, double Quality, const char * InsBases="");
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, double Quality, Segment Read1, Segment Read2, int Length);//for drp Signatures
    Signature(int Type, int Tech, int SupportedSV, int Begin, int End, std::string TemplateName, double Quality, std::vector<Segment> Segments);//for Split reads signatures
    void setCN(int CN);
    void setInvLeft(bool);
    void setInvRight(bool);
    bool operator<(const Signature &Other) const;
    bool operator==(const Signature & Other) const;
    void setID(unsigned ID);
    operator std::string();
    #ifdef DEBUG
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & Type;
        ar & Tech;
        ar & SupportedSV;
        ar & TemplateName;
        ar & InsBases;
        ar & Begin;
        ar & End;
        ar & CN;
        ar & Segments;
        ar & Length;
        ar & InvLeft;
        ar & InvRight;
        ar & Quality;
        ar & ID;
        ar & Artificial;
        ar & HP;
        ar & MergeString;
    }
    void setMergeString(std::string MS)
    {
        MergeString=MS;
    }
    #endif
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
    void merge(SegmentSet& Other);//will purge Other.
    void sortNStat();
    Segment& operator[](unsigned i);
    SegmentSet();
    #ifdef DEBUG
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & Segments;
        ar & MaxEnds;
    }
    #endif
};

#endif