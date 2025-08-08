#include "signature.h"
#include <string.h>
#include "crelib/crelib.h"
#include <algorithm>
using namespace std;

Signature::Signature()
:Type(0), Tech(0), SupportedSV(0), Segment(), CN(-1), Covered(false), TemplateName(""), Quality(0), Length(0), InvLeft(false), InvRight(false), Artificial(false), HP(0)DEBUG_CODE(COMMA MergeString(""))
{}
Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, double Quality, const char * InsBases)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Segment({Begin,End,InsBases}), CN(-1), Covered(false), TemplateName(TemplateName), Quality(Quality), Length(End-Begin), InvLeft(false), InvRight(false), Artificial(false), HP(0)DEBUG_CODE(COMMA MergeString("")){}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, double Quality, Segment Read1, Segment Read2, int ALength)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Segment({Begin,End}), CN(-1), Covered(false), TemplateName(TemplateName), Quality(Quality), Segments(), Length(ALength), InvLeft(false), InvRight(false), Artificial(false), HP(0)DEBUG_CODE(COMMA MergeString("")) {Segments.push_back(Read1), Segments.push_back(Read2);}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, double Quality, vector<Segment> Segments)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Segment({Begin,End}), CN(-1), Covered(false), TemplateName(TemplateName), Quality(Quality), Segments(Segments), Length(End-Begin), InvLeft(false), InvRight(false), Artificial(false), HP(0)DEBUG_CODE(COMMA MergeString("")) {}

const char* Signature::SVTypeNames[]={"DEL","INS","DUP","INV"};

void Signature::setCN(int cn)
{
    CN=cn;
}

void Signature::setInvLeft(bool b)
{
    InvLeft=b;
}

void Signature::setInvRight(bool b)
{
    InvRight=b;
}

void Signature::setID(unsigned d)
{
    ID=d;
}

bool Signature::operator<(const Signature &Other) const
{
    if (this->Begin<Other.Begin) return true;
    else if (this->Begin>Other.Begin) return false;
    else if (this->Length<Other.Length) return true;
    else if (this->Length>Other.Length) return false;
    else if (this->Quality>Other.Quality) return true;
    else if (this->Quality<Other.Quality) return false;
    else if (this->InsBases<Other.InsBases) return true;
    else if (this->InsBases>Other.InsBases) return false;
    else return this->TemplateName<Other.TemplateName;
}

Signature::operator std::string()
{
    return TemplateName+", "+std::string(*((Segment*)this))+", "+std::to_string(Type)+" "+std::to_string(Tech)+" "+std::to_string(SupportedSV)+" "+std::to_string(CN)+" "+std::to_string(Covered)+" "+std::to_string(Tech)+" "+std::to_string(Length)+" "+std::to_string(Quality)+" "+std::to_string(Artificial);
}

bool Signature::operator==(const Signature &Other) const
{
    return Tech==Other.Tech && SupportedSV==Other.SupportedSV && TemplateName==Other.TemplateName && Begin==Other.Begin && End==Other.End && Segments==Other.Segments;
}

Segment::Segment(int Begin,int End, const char * CInsBases)
:Begin(Begin),End(End),InsBases(CInsBases){}

bool Segment::operator==(const Segment & Other) const
{
    return Begin==Other.Begin && End==Other.End && InsBases==Other.InsBases;
}

bool Segment::operator<(const Segment & Other) const
{
    return Begin<Other.Begin;
}

Segment::operator std::string()
{
    return "("+std::to_string(Begin)+", "+std::to_string(End)+"), ("+InsBases+")";
}

int precisionLevel(const Signature &A)
{
    if (A.Tech==1 && A.Type==1) return 0;//drp sig, imprecise pricision
    if (A.Tech==0) return 1;//SMRT CIGAR and CLIP, vague precision
    return 2;//NGS CIGAR and clip, precise precision
}

int bestPrecision(const Signature &A,const Signature &B)
{
    return MAX(precisionLevel(A),precisionLevel(B));
}

int worstPrecision(const Signature &A,const Signature &B)
{
    return MIN(precisionLevel(A),precisionLevel(B));
}

std::tuple<int,int> SegmentSet::getInvolved(int Start, int End)
{
    tuple<int,int> Result={0,0};
    int SS=0, SE=Segments.size(), Mid;
    while (SS<SE-1)
    {
        Mid=(SS+SE)/2;
        if (MaxEnds[Mid]<Start) SS=Mid;
        else SE=Mid;
    }
    get<0>(Result)=SS;
    SS=0, SE=Segments.size();
    while (SS<SE-1)
    {
        Mid=(SS+SE)/2;
        if (Segments[Mid].Begin>End) SE=Mid;
        else SS=Mid;
    }
    get<1>(Result)=SE;
    return Result;
}
void SegmentSet::add(int Begin, int End)
{
    Segments.push_back(Segment(Begin,End));
}
void SegmentSet::merge(SegmentSet& Other)
{
    Segments.insert(Segments.end(),make_move_iterator(Other.Segments.begin()),make_move_iterator(Other.Segments.end()));
}
void SegmentSet::sortNStat()
{
	sort(Segments.begin(),Segments.begin()+Segments.size());
    if (Segments.size()>0)
    {
        MaxEnds.resize(Segments.size());
        int MaxEnd=Segments[0].End;
        for (int i=0;i<Segments.size();++i)
        {
            MaxEnd=MAX(MaxEnd,Segments[i].End);
            MaxEnds[i]=MaxEnd;
        }
    }
}
Segment& SegmentSet::operator[](unsigned i)
{
    return Segments[i];
}

SegmentSet::SegmentSet():Segments(),MaxEnds(){}