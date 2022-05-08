#include "signature.h"
#include <string.h>
#include "crelib/crelib.h"
#include <algorithm>
using namespace std;

Signature::Signature()
{}
Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, const char * InsBases)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Length(End-Begin), InsBases(InsBases), InvLeft(false), InvRight(false){}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, Segment Read1, Segment Read2, int ALength)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Segments(), Length(ALength) {Segments.push_back(Read1), Segments.push_back(Read2);}

Signature::Signature(int Type, int Tech, int SupportedSV, int Begin, int End, string TemplateName, vector<Segment> Segments)
:Type(Type), Tech(Tech), SupportedSV(SupportedSV), Begin(Begin), End(End), CN(-1), TemplateName(TemplateName), Segments(Segments), Length(End-Begin) {}

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

bool Signature::operator<(const Signature &Other) const
{
    char AL[100], BL[100];
    sprintf(AL,"%d",this->Length);
    sprintf(BL,"%d",Other.Length);
    int cmp=strcmp(AL,BL);
    if (this->Begin<Other.Begin) return true;
    else if (this->Begin>Other.Begin) return false;
    // else if (this->Length<Other.Length) return true;
    // else if (this->Length>Other.Length) return false;
    else if (cmp<0) return true;
    else if (cmp>0) return false;
    else return this->TemplateName<Other.TemplateName;
}

bool Signature::operator==(const Signature &Other) const
{
    return Tech==Other.Tech && SupportedSV==Other.SupportedSV && TemplateName==Other.TemplateName && Begin==Other.Begin && End==Other.End && Segments==Other.Segments;
}

Segment::Segment(int Begin,int End)
:Begin(Begin),End(End){}

bool Segment::operator==(const Segment & Other) const
{
    return Begin==Other.Begin && End==Other.End;
}

bool Segment::operator<(const Segment & Other) const
{
    return Begin<Other.Begin;
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