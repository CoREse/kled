#include "report.h"
#include <assert.h>
#include "clustering.h"
#include <tuple>
#include <algorithm>
#include <math.h>
#include "crelib/crelib.h"
#include "time.h"
#include <set>
#include <numeric>
#include <stdio.h>
#include <iterator>
#include <gmp.h>

using namespace std;

float stdStat(Signature* Begin,Signature* End)
{
    float Mean=0;
    int k=0;
    for (Signature * i=Begin;i!=End;++i)
    {
        Mean+=i->Begin;
        ++k;
    }
    Mean/=float(k);
    float std=0;
    for (Signature * i=Begin;i!=End;++i)
    {
        std+=pow(i->Begin-Mean,2);
    }
    std/=float(k);
    return pow(std,0.5);
}

tuple<int,int> consensusFixed(vector<Signature> &Signatures, int K, float stat(Signature*,Signature*))//计算标准差最小的连续K个的pos和length的均值
{
    if (Signatures.size()==1) return tuple<int,int>(Signatures[0].Begin,(Signatures[0].End-Signatures[0].Begin));
    float MinStd=0;
    int MinI=-1;
    for (int i=0;i<Signatures.size()-K;++i)
    {
        float Std=stat(Signatures.data()+i,Signatures.data()+i+K);
        if (MinI==-1 || MinStd<Std)
        {
            MinI=i;
            MinStd=Std;
        }
    }
    float Pos=0,Length=0;
    for (int i=MinI;i<MinI+K;++i)
    {
        Pos+=Signatures[i].Begin;
        Length+=Signatures[i].Length;
    }
    return tuple<int,int>(int(Pos/float(K)),int(Length/float(K)));
}

const char * getSVType(const vector<Signature> &Signatures)
{
    return Signature::SVTypeNames[Signatures[0].SupportedSV];
}

int VCFRecord::getSVTypeI() const
{
    return SVTypeI;
}

int VCFRecord::getPSD() const
{
    return PSD;
}

int VCFRecord::getSVLen() const
{
    return SVLen;
}

int VCFRecord::getST() const
{
    return ST;
}

tuple<int,int> getConsensusSite(vector<Signature> &Signatures)
{
    if (Signatures.size()==0) return tuple<int,int>(-1,-1);
    sort(Signatures.data(),Signatures.data()+Signatures.size());
    return consensusFixed(Signatures,Signatures.size()/2,stdStat);
}

tuple<int,int,int> analyzeSignatureCluster(vector<Signature> &SignatureCluster, string SVType, Arguments &Args)
{
    
    unsigned long long llPos=0;
    unsigned long long llSVLen=0;
    int Pos=0, SVLen=0, MeanSVLen=0;
    int CBegin=0, CEnd=SignatureCluster.size();
    // if (SVType=="INS")
    // {
    //     int MinBinSize=10;
    //     int BinNumber=100;
    //     int MaxSVLen=0, MinSVLen=0xffffff;
    //     for (int i=CBegin;i<CEnd;++i)
    //     {
    //         llPos+=SignatureCluster[i].Begin;
    //         MaxSVLen=max(MaxSVLen,SignatureCluster[i].Length);
    //         MinSVLen=min(MinSVLen,SignatureCluster[i].Length);
    //         llSVLen+=SignatureCluster[i].Length;
    //     }
    //     Pos=llPos/(CEnd-CBegin);
    //     MeanSVLen=llSVLen/(CEnd-CBegin);
    //     int BinSize=(MaxSVLen-MinSVLen)/BinNumber+1;
    //     if (BinSize<MinBinSize) SVLen=MeanSVLen;
    //     else
    //     {
    //         vector<vector<int>> SigBins;
    //         for (int i=0;i<BinNumber;++i)
    //         {
    //             SigBins.push_back(vector<int>());
    //         }
    //         for (int i=CBegin;i<CEnd;++i)
    //         {
    //             SigBins[(SignatureCluster[i].Length-MinSVLen)/BinSize].push_back(i);
    //         }
    //         int MaxQualifiledIndex=-1;
    //         for (int i=0;i<BinNumber;++i)
    //         {
    //             if (SigBins[i].size()<3) continue;
    //             MaxQualifiledIndex=i;
    //         }
    //         if (MaxQualifiledIndex!=-1)
    //         {
    //             llSVLen=0;
    //             for (int i=0;i<SigBins[MaxQualifiledIndex].size();++i)
    //             {
    //                 llSVLen+=SignatureCluster[SigBins[MaxQualifiledIndex][i]].Length;
    //             }
    //             SVLen=llSVLen/SigBins[MaxQualifiledIndex].size();
    //         }
    //         else SVLen=MeanSVLen;
    //     }
    //     // for (int i=CBegin;i<CEnd;++i)
    //     // {
    //     //     llPos+=SignatureCluster[i].Begin;
    //     //     llSVLen+=SignatureCluster[i].Length;
    //     // }
    //     // Pos=llPos/(CEnd-CBegin);
    //     // MeanSVLen=llSVLen/(CEnd-CBegin);
    //     // for (int i=CBegin;i<CEnd;++i)
    //     // {
    //     //     if (SignatureCluster[i].Length<MeanSVLen*2)
    //     //     {
    //     //         SVLen=max(SVLen,SignatureCluster[i].Length);
    //     //     }
    //     // }
    // }
    // else
    if (Args.WeightPosLength)
    {
        double WeightSum=0.0;
        for (int i=CBegin;i<CEnd;++i)
        {
            WeightSum+=SignatureCluster[i].Quality;
        }
        double dPos=0, dMeanSVLen=0;
        for (int i=CBegin;i<CEnd;++i)
        {
            dPos+=((double)SignatureCluster[i].Begin)*(SignatureCluster[i].Quality/WeightSum);
            dMeanSVLen+=((double)SignatureCluster[i].Length)*(SignatureCluster[i].Quality/WeightSum);
        }
        Pos=dPos;
        MeanSVLen=dMeanSVLen;
        SVLen=MeanSVLen;
    }
    else
    {
        for (int i=CBegin;i<CEnd;++i)
        {
            llPos+=SignatureCluster[i].Begin;
            llSVLen+=SignatureCluster[i].Length;
        }
        Pos=llPos/(CEnd-CBegin);
        MeanSVLen=llSVLen/(CEnd-CBegin);
        SVLen=MeanSVLen;
    }
    return tuple<int,int,int>(Pos,SVLen,MeanSVLen);
    float L2Weight=0.75;
    int Length;
    vector<Signature> PSigs[3];
    for (int i=0;i<SignatureCluster.size();++i)
    {
        PSigs[precisionLevel(SignatureCluster[i])].push_back(SignatureCluster[i]);
    }
    tuple<int,int> Sites[3];
    for (int i=0;i<3;++i) Sites[i]=getConsensusSite(PSigs[i]);
    if (get<0>(Sites[2])!=-1)//has pricise level signatures
    {
        Pos=get<0>(Sites[2]);
        Length=get<1>(Sites[2]);
    }
    else
    {
        if (get<0>(Sites[1])==-1)//doesn't have vague level signatures
        {
            Pos=get<0>(Sites[0]);
            Length=get<1>(Sites[0]);
        }
        else if (get<0>(Sites[0])==-1)//doesn't have imprecise level signatures
        {
            Pos=get<0>(Sites[1]);
            Length=get<1>(Sites[1]);
        }
        else
        {
            Pos=get<0>(Sites[0])*(1.0-L2Weight)+get<0>(Sites[1])*(L2Weight);
            Length=get<1>(Sites[0])*(1.0-L2Weight)+get<1>(Sites[1])*(L2Weight);
        }
    }
    return tuple<int,int,int>(Pos,Length,MeanSVLen);
}

template <typename IterType> double calcSD(IterType Begin, IterType End)//Not estimate, just calc
{
    int Size=0;
    double Mean=0;
    for (IterType it=Begin;it!=End;++it)
    {
        Mean+=*it;
        ++Size;
    }
    Mean/=(double)Size;
    double SD=0;
    for (IterType it=Begin;it!=End;++it)
    {
        SD+=pow((*it)-Mean,2.0);
    }
    SD/=(double)Size;
    return pow(SD,0.5);
}

template <typename ValueType,typename IterType> class MemberIter
{
    protected:
    IterType BaseIter;
    public:
    MemberIter(IterType A):BaseIter(A){}
    MemberIter& operator++()
    {
        ++BaseIter;
        return *this;
    }
    bool operator!=(const MemberIter & Other)
    {
        return BaseIter!=Other.BaseIter;
    }
    bool operator!=(const IterType & Other)
    {
        return BaseIter!=Other;
    }
    // virtual ValueType operator* ();//Can't use virtual in template!
};

template <typename ValueType,typename IterType> class BeginIter : public MemberIter<ValueType,IterType>
{
    public:
    // using MemberIter<ValueType,IterType>::MemberIter<ValueType,IterType>;
    BeginIter(IterType A):MemberIter<ValueType,IterType>(A){}
    ValueType operator* ()
    {
        return this->BaseIter->Begin;
    }
};

template <typename ValueType,typename IterType> class EndIter : public MemberIter<ValueType,IterType>
{
    public:
    EndIter(IterType A):MemberIter<ValueType,IterType>(A){}
    ValueType operator* ()
    {
        return this->BaseIter->End;
    }
};
template <typename ValueType,typename IterType> class LengthIter : public MemberIter<ValueType,IterType>
{
    public:
    LengthIter(IterType A):MemberIter<ValueType,IterType>(A){}
    ValueType operator* ()
    {
        return this->BaseIter->Length;
    }
};

vector<double> scoring(vector<Signature> &SignatureCluster,int SVLen)
{
    SVLen=abs(SVLen);
    vector<double> Scores;//supported reads, supported signatures, (begin sd+end sd)/2, length sd, number of signature source
    int SS=0;
    set<string> SupportTemps;;
    for (int i =0;i<SignatureCluster.size();++i)
    {
        ++SS;
        SupportTemps.insert(SignatureCluster[i].TemplateName);
    }
    int ST=SupportTemps.size();
    double STS=ST>50?100.0:double(ST)/50.0*100.0;
    double SSS=SS>50?100.0:double(SS)/50.0*100.0;
    double ASS=(STS+SSS)/2.0;
    Scores.push_back(ASS);
    double BESD=0.0;//calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    double LengthSD=calcSD(LengthIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),LengthIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    double BESDS=0.0;//BESD>50?100.0:double(BESD)/50.0*100.0;BESDS=100.0-BESDS;
    double LengthSDS=0.0;//LengthSD>50?100.0:double(LengthSD)/50.0*100.0;LengthSDS=100.0-LengthSDS;
    Scores.push_back(BESDS);
    double BESDRatioS=0.0;//BESDRatio>0.5?100:BESDRatio/0.5*100.0;BESDRatioS=100.0-BESDRatioS;
    Scores.push_back(BESDRatioS);
    Scores.push_back(LengthSDS);
    double LengthSDRatio=LengthSD/(double)SVLen;
    double LengthSDRatioS=LengthSDRatio>0.5?100:LengthSDRatio/0.5*100.0;LengthSDRatioS=100.0-LengthSDRatioS;
    Scores.push_back(LengthSDRatioS);
    Scores.push_back(0.0);//double(Sources.size())/3.0*100.0);
    return Scores;
}

int getLengthSDRatioScore(vector<Signature> &SignatureCluster, int SVLen, double * SD=NULL)
{
    double LengthSD=calcSD(LengthIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),LengthIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    if (SD!=NULL) *SD=LengthSD;
    return MAX(0.0,100.0-(LengthSD/(double)SVLen*100.0));
    // return MAX(0.0,100.0-(LengthSD/(double)(max(SVLen,500))*100.0));
    int SizePenalty=0;
    // if (SignatureCluster.size()<4) SizePenalty=4-SignatureCluster.size();
    // return MAX(0.0,100.0-(LengthSD/(double)SVLen*500.0-SizePenalty));
}

double getAverageCoverage(int Begin, int End, float * CoverageWindows, Arguments & Args, float* CoverageWindowsSums, float* CheckPoints, int CheckPointInterval)
{
    if (End<=Begin) return -1;
    double Cov=0;
    int WBegin=Begin/Args.CoverageWindowSize;
    int WEnd=End/Args.CoverageWindowSize+1;
    if (false and CoverageWindowsSums!=NULL)//buggy, guess related to checkpoint memory access cross boundary
    {
        if ((WBegin/CheckPointInterval+1)*CheckPointInterval>WEnd)
        {
            Cov=(CoverageWindowsSums[WEnd]-CoverageWindowsSums[WBegin])/(double)(WEnd-WBegin);
        }
        else
        {
            int Range=0;
            for (int i=WBegin/CheckPointInterval+2;i<=WEnd/CheckPointInterval;++i)
            {
                Cov*=(double)(Range)/(double)(Range+CheckPointInterval);
                Cov+=CheckPoints[i]/(double)(Range+CheckPointInterval);
                Range+=CheckPointInterval;
            }
            Cov*=(double)(Range)/(Range+((double)(CheckPointInterval-(WBegin%CheckPointInterval))));
            Cov+=(CheckPoints[(int)WBegin/CheckPointInterval+1]-CoverageWindowsSums[WBegin])/(double)(Range+(CheckPointInterval-(WBegin%CheckPointInterval)));
            Range+=(CheckPointInterval-(WBegin%CheckPointInterval));
            if (WEnd%CheckPointInterval!=0)
            {
                Cov*=(double)(Range)/(double)(Range+(WEnd%CheckPointInterval));
                Cov+=CoverageWindowsSums[WEnd]/(double)(Range+(WEnd%CheckPointInterval));
                Range+=(WEnd%CheckPointInterval);
            }
        }
    }
    else
    {
        for (int i=WBegin;i<WEnd;++i) 
        {
            Cov*=((float)(i-WBegin))/(i-WBegin+1);//in case too big
            Cov+=CoverageWindows[i]*1/(i-WBegin+1);
        }
    }
    return Cov;
}

bool statCluster(vector<Signature> &SignatureCluster, int & SS, int &ST, int& SS2, int& ST2, bool ExcludeHP0=false)
{
    SS=0;
    SS2=0;
    set<string> SupportTemps,SupportTemps2;//,SupportTempsCigar;
    if (ExcludeHP0)
    {
        for (int i =0;i<SignatureCluster.size();++i)
        {
            if (SignatureCluster[i].HP==0) continue;
            ++SS;
            SupportTemps.insert(SignatureCluster[i].TemplateName);
            if (SignatureCluster[i].Tech==1)
            {
                ++SS2;
                SupportTemps2.insert(SignatureCluster[i].TemplateName);
            }
        }
    }
    else
    {
        for (int i =0;i<SignatureCluster.size();++i)
        {
            ++SS;
            SupportTemps.insert(SignatureCluster[i].TemplateName);
            if (SignatureCluster[i].Tech==1)
            {
                ++SS2;
                SupportTemps2.insert(SignatureCluster[i].TemplateName);
            }
        }
    }
    ST=SupportTemps.size();
    ST2=SupportTemps2.size();
    return true;
}

bool isClipOnly(vector<Signature> & Cluster)
{
    bool HasClip=false, OnlyClip=true;
    for (int i=0;i<Cluster.size();++i)
    {
        if (Cluster[i].Type==2) HasClip=true;
        else OnlyClip=false;
    }
    return HasClip && OnlyClip;
}

double err=0.1;
double prior=1.0/3.0;
string calc_GL_cute(int c0, int c1)
{
    int MaxAllowed=100;
    if (c0+c1>MaxAllowed)
    {
        c0=((double)MaxAllowed)*((double)c0)/((double)(c0+c1));
        c1=MaxAllowed-c0;
    }
    double ori_GL00=pow((1.0-err),c0)*pow(err,c1)*(1.0-prior)/2.0;
    double ori_GL11=pow(err,c0)*pow(1.0-err,c1)*(1.0-prior)/2.0;
    double ori_GL01=pow(0.5, c0+c1)*prior;
    // if (ori_GL11>ori_GL01 && ori_GL11>ori_GL00) return "1/1";
    // return "0/1";
    double log10Ps[3]={log10(ori_GL00),log10(ori_GL01),log10(ori_GL11)};
    double MaxLog10=-log10Ps[0];
    for (int i=1;i<3;++i)
    {
        if (MaxLog10<log10Ps[i]) MaxLog10=log10Ps[i];
    }
    double sum=0;
    for (int i=0;i<3;++i)
    {
        sum+=pow(10,log10Ps[i]-MaxLog10);
    }
    double lse=MaxLog10+log10(sum);
    double np[3];
    for (int i=0;i<3;++i) np[i]=MIN(log10Ps[i]-lse,0.0);
    // fprintf(stderr, "%lf,%lf,%lf\n,",np[0],np[1],np[2]);
    if (np[2]>np[1] && np[2]>np[0]) return "1/1";
    return "0/1";
}

unsigned long long combination(unsigned long long Total, unsigned long long Selected)
{
    unsigned long long Result=1;
    for (unsigned long long i=Total; i>Selected;--i)
    {
        Result*=i;
    }
    return Result;
}

const double pi=3.14159265358979323846;
double normal(double mu, double sigma, double x)
{
    return exp(-pow(x-mu,2)/(2*sigma*sigma))/(2*pi*sigma);
}

double poisson(double lambda, int x)
{
    double p=1;
    for (int i=1;i<=x;++i) p*=lambda/i;
    p*=exp(-lambda);
    return p;
}

int getNAllTemplates(vector<Signature> & Cluster,const Contig & TheContig, vector<Sam>& SamFiles, int Start, int End)
{
    // for (auto s : Cluster) {Start=MIN(s.Begin,Start);End=MAX(s.Begin,End);}
    set<string> Temps;
    for (int k=0;k<SamFiles.size();++k)
    {
        bam1_t *br=bam_init1();
        hts_itr_t* RegionIter=sam_itr_querys(SamFiles[k].BamIndex,SamFiles[k].Header,(TheContig.Name+":"+to_string(Start)+"-"+to_string(End)).c_str());
        while(sam_itr_next(SamFiles[k].SamFile, RegionIter, br) >=0)//read record
        {
	        if (!align_is_primary(br)) continue;
            // for (auto s: Cluster)
            if (br->core.pos<Start && br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br))>=End) Temps.insert(bam_get_qname(br));
        }
        bam_destroy1(br);
    }
    return Temps.size();
}

int getNAllTemplates(const Contig & TheContig, SegmentSet & AllPrimarySegments, int Start, int End)
{
    // for (auto s : Cluster) {Start=MIN(s.Begin,Start);End=MAX(s.Begin,End);}
    // set<string> Temps;
    tuple<int, int> SearchRange=AllPrimarySegments.getInvolved(Start,End);
    int Result=0;
    for (int i=get<0>(SearchRange);i<get<1>(SearchRange);++i)
    {
        if (AllPrimarySegments[i].Begin<=Start && AllPrimarySegments[i].End>=End) ++Result;
    }
    return Result;
}

string VCFRecord::genotype(const Contig & TheContig, SegmentSet & AllPrimarySegments, float * CoverageWindows, float *CoverageWindowsSums, float* Checkpoints, int CheckPointInterval, Arguments & Args)
{
    double HomoThreshold=0.8;
    int Scope=0;
    int CountStart=Pos;
    int CountEnd=Pos+SVLen;
    if (SVType=="INS")
    {
        Scope=400;
        CountStart=MAX(Pos-Scope,0);
        CountEnd=MIN(Pos+Scope,TheContig.Size-1);
        HomoThreshold=0.75;
    }
    int All;
    if (SVType=="INS") All=getNAllTemplates(TheContig, AllPrimarySegments,CountStart,CountEnd);
    else All=getAverageCoverage(CountStart,CountEnd,CoverageWindows,Args, CoverageWindowsSums, Checkpoints, CheckPointInterval);
    CV=All;
    if (All==0) Sample["GT"]="1/1";
    else if (double(ST)/All<HomoThreshold) Sample["GT"]="0/1";
    else Sample["GT"]="1/1";
    return Sample["GT"];
}

void resizeCluster(vector<Signature> &Cluster, int MaxSize)
{
    if (Cluster.size()<=MaxSize) return;
    int NewBegin=Cluster.size()/2-(MaxSize/2);
    int NewEnd=NewBegin+MaxSize;
    Cluster.insert(Cluster.begin(),make_move_iterator(vector<Signature>::iterator(&Cluster[NewBegin])),make_move_iterator(vector<Signature>::iterator(&Cluster[NewEnd])));
    Cluster.resize(MaxSize);
    // if (Cluster.size()<=2*MaxSize) return;
    // int Step=Cluster.size()/MaxSize;
    // auto InitialFirst=Cluster.begin();
    // int NewSize=0;
    // for (int i=0;i<Cluster.size();i+=Step)
    // {
    //     Cluster.insert(vector<Signature>::iterator(&Cluster[NewSize]),make_move_iterator(vector<Signature>::iterator(&Cluster[i])),make_move_iterator(i<Cluster.size()?vector<Signature>::iterator(&Cluster[i+1]):Cluster.end()));
    //     ++NewSize;
    // }
    // Cluster.resize(NewSize);
}

const short Alphabet[128]={//NATGC=01234
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,4,0,0,0,3,0,0,0,0,0,0,0,0,
    0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,4,0,0,0,3,0,0,0,0,0,0,0,0,
    0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0
};
const char * Betalpha="NATGC";
string getInsConsensus(int SVLen, vector<Signature> & SignatureCluster, double Endurance=0.01)
{
    char Consensus[SVLen+1];
    Consensus[0]='\0';
    int Counts[5];
    int ConsensusCount=0;
    if (ConsensusCount==0)
    {
        for (int i=0;i<SignatureCluster.size();++i)
        {
            if (SignatureCluster[i].InsBases.length()>=SVLen)
            {
                strcpy(Consensus,SignatureCluster[i].InsBases.substr(0,SVLen).c_str());
                break;
            }
        }
    }
    else
    {
        for (int i=0;i<SVLen;++i)
        {
            for (int j=0;j<5;++j) Counts[j]=0;
            for (int j=0;j<SignatureCluster.size();++j)
            {
                if (abs(int(SignatureCluster[j].InsBases.length())-SVLen)>SVLen*Endurance) continue;
                // if (SignatureCluster[j].Length!=SignatureCluster[j].InsBases.length()) fprintf(stderr,"warning! not the same!%d,%d\n",SignatureCluster[j].Length,SignatureCluster[j].InsBases.length());
                ++Counts[
                    Alphabet[
                        SignatureCluster[j].InsBases[(int)(((double)i)/((double)SVLen)*SignatureCluster[j].InsBases.length())]
                        ]
                    ];
            }
            int MaxI=0;
            for (int j=1;j<5;++j) if (Counts[MaxI]<Counts[j]) MaxI=j;
            Consensus[i]=Betalpha[MaxI];
        }
    }
    Consensus[SVLen]='\0';
    return Consensus;
}

void VCFRecord::calcM3L(vector<Signature> & SignatureCluster, bool ExcludeHP0)
{
    vector<int> Lengths;
    if (ExcludeHP0)
    {
        for (auto s : SignatureCluster)
        {
            if (s.HP==0) continue;
            Lengths.push_back(s.Length);
        }
    }
    else
    {
        for (auto s : SignatureCluster)
        {
            Lengths.push_back(s.Length);
        }
    }
    sort(Lengths.data(),Lengths.data()+Lengths.size());
    MinLength=Lengths[0];
    MediumLength=Lengths[(int)(Lengths.size()/2)];
    MaxLength=Lengths[Lengths.size()-1];
}

double getWeightSum(vector<Signature> & SignatureCluster)
{
    map<string,double> HighstTempQuality;
    for (int i=0;i<SignatureCluster.size();++i)
    {
        if (HighstTempQuality.count(SignatureCluster[i].TemplateName)==0)
        {
            HighstTempQuality[SignatureCluster[i].TemplateName]=SignatureCluster[i].Quality;
        }
        else
        {
            HighstTempQuality[SignatureCluster[i].TemplateName]=max(HighstTempQuality[SignatureCluster[i].TemplateName],SignatureCluster[i].Quality);
        }
    }
    double WeightSum=0;
    for (auto iter=HighstTempQuality.begin();iter!=HighstTempQuality.end();++iter)
    {
        WeightSum+=iter->second;
    }
    return WeightSum;
}

inline void biCdf(mpq_t Result, unsigned n, const mpq_t p, unsigned Start, unsigned End)
{
    mpq_t One,Pq,Subp;
    mpz_t P;
    mpq_set_ui(Result,0,1);
    mpz_init(P);
    mpq_init(One);
    mpq_init(Pq);
    mpq_init(Subp);
    mpq_set_ui(One,1,1);
    mpq_sub(Subp,One,p);
    for (int i=Start;i<End;++i)
    {
        mpz_set_ui(P,1);
        for (int j=n-i+1;j<=n;++j) mpz_mul_ui(P,P,j);
        for (int j=1;j<=i;++j) mpz_tdiv_q_ui(P,P,j);
        mpq_set_z(Pq,P);
        for (int j=0;j<i;++j) mpq_mul(Pq,Pq,p);
        for (int j=i;j<n;++j) mpq_mul(Pq,Pq,Subp);
        mpq_add(Result,Result,Pq);
    }
}

int VN=0;
VCFRecord::VCFRecord()
: SVLen(0),SVType(),SS(0),ST(0),SS2(0),ST2(0),LS(0),CV(0),CR(0),MinLength(0),MaxLength(0),MediumLength(0),Precise(0),InsConsensus(""),SVTypeI(0),CHROM(),Pos(0),ID(),REF(),ALT(),QUAL(),FILTER(),INFO(),Sample(),Keep(false)
{}
VCFRecord::VCFRecord(const Contig & TheContig,vector<Signature> & SignatureCluster, ClusterCore &Core, SegmentSet & AllPrimarySegments, float* CoverageWindows, double WholeCoverage, Arguments& Args, float * CoverageWindowsSums, float * CheckPoints, int CheckPointInterval)
: SVLen(0),SVType(),CN(0),SS(0),ST(0),SS2(0),ST2(0),LS(0),CV(0),CR(0),MinLength(0),MaxLength(0),MediumLength(0),Precise(0),InsConsensus(""),SVTypeI(0),CHROM(),Pos(0),ID(),REF(),ALT(),QUAL(),FILTER(),INFO(),Sample(),PSD(-1),Keep(false)
{
    SVType=getSVType(SignatureCluster);
    SVTypeI=SignatureCluster[0].SupportedSV;
    double HPRatio=Args.HPRatio[SVTypeI], HomoRatio=Args.HomoRatio[SVTypeI], HomoMinusRatio=Args.HomoMinusRatio[SVTypeI];
    // if (SVTypeI==1) {
        // HPRatio=0.6, HomoRatio=0.8, HomoCutoffRatio=0.9;
    // }
    double CutOffRatio=1.0;
    assert(SignatureCluster.size()>=0);
    resizeCluster(SignatureCluster,Args.MaxClusterSize);
    Keep=true;
    // int HPCounts[3]={0,0,0};
    HPCounts[0]=0;HPCounts[1]=0;HPCounts[2]=0;
    for (int i=0;i<SignatureCluster.size();++i)
    {
        ++HPCounts[SignatureCluster[i].HP];
    }
    //HP GT hypothesis testing
    bool HPGTKeep=false;
    // int HPCount=HPCounts[1]+HPCounts[2];
    // double HPPortion=((double)HPCount)/((double)(HPCounts[0]+HPCount));
    // double Alpha=0.05;
    // double PRandomOrHomo=0;
    // int MinHP=min(HPCounts[1],HPCounts[2]);
    // double MultipleP=pow(0.5,HPCount);
    // for (int i=0;i<=MinHP;++i)
    // {
    //     // tgamma(i+1)==i!
    //     //binominal distribution
    //     unsigned long long P=1;
    //     for (int j=HPCount-i+1;j<=HPCount;++j) P*=j;
    //     for (int j=1;j<=i;++j) P/=j;
    //     double dP=((double)P)*MultipleP;
    //     PRandomOrHomo+=dP;
    // }
    // if (PRandomOrHomo<Alpha)
    // {
    //     // HPGTKeep=true;
    //     Sample["GT"]="0/1";
    //     // if (HPCounts[1]>HPCounts[2]) Sample["GT"]="1|0";
    //     // else Sample["GT"]="0|1";
    // }

    //HP Ratio
    double SHPRatio=float(HPCounts[1]+HPCounts[2])/float(HPCounts[0]+HPCounts[1]+HPCounts[2]);
    double GTHPRatio=0.8, GTHomoRatio=0.6;
    float HPRatios[3]={0,0,0};
    if (SHPRatio>HPRatio || SHPRatio>GTHPRatio)
    {
        for (int i=0;i<3;++i) HPRatios[i]=float(HPCounts[i])/float(SignatureCluster.size());
    }
    // if (SHPRatio>GTHPRatio)
    // {
    //     // double Alpha=0.0005;
    //     mpq_t Alpha;
    //     mpq_init(Alpha);
    //     mpq_set_str(Alpha,"5/100000",10);
    //     // fprintf(stderr, mpq_get_str(NULL,10,Alpha));
    //     // exit(0);
    //     // double PRandomOrHomo=0;
    //     int HPCount=HPCounts[1]+HPCounts[2];
    //     int MinHP=min(HPCounts[1],HPCounts[2]);
    //     mpq_t MinorP,MajorP,One;
    //     mpq_init(MinorP);
    //     mpq_init(MajorP);
    //     mpq_init(One);
    //     mpq_set_ui(One,1,1);
    //     // mpq_t SubAlpha;
    //     // mpq_init(SubAlpha);
    //     // mpq_sub(SubAlpha,One,Alpha);
    //     mpq_set_str(MinorP,"15/1000",10);
    //     mpq_sub(MajorP,One,MinorP);
    //     // if (HPCount>=30 && MinHP>0.3*HPCount) Sample["GT"]="1/1";
    //     // else if (HPCount>=30)
    //     {
    //         // mpq_t MultipleP;
    //         // mpz_t MPBase;
    //         // mpz_init(MPBase);
    //         // mpq_init(MultipleP);
    //         // mpz_set_ui(MPBase,2);
    //         // mpz_pow_ui(MPBase,MPBase,HPCount);
    //         // mpq_set_z(MultipleP,MPBase);
    //         // mpq_inv(MultipleP,MultipleP);
    //         // mpz_t P;
    //         // mpz_init(P);
    //         // mpq_t Pq;
    //         // mpq_init(Pq);
    //         mpq_t PRandom;
    //         mpq_init(PRandom);
    //         mpq_set_ui(PRandom,0,1);
    //         // double MultipleP=pow(0.5,HPCount);
    //         // for (int i=0;i<=MinHP;++i)
    //         biCdf(PRandom,HPCount,MinorP,MinHP,HPCount);
    //         // for (int i=MinHP;i<HPCount;++i)
    //         // {
    //         //     // tgamma(i+1)==i!
    //         //     //binominal distribution

    //         //     mpz_set_ui(P,1);
    //         //     for (int j=HPCount-i+1;j<=HPCount;++j) mpz_mul_ui(P,P,j);
    //         //     for (int j=1;j<=i;++j) mpz_tdiv_q_ui(P,P,j);
    //         //     //Now P=C(HPCount,i)
    //         //     mpq_set_z(Pq,P);
    //         //     // mpq_mul(Pq,Pq,MultipleP);
    //         //     for (int j=0;j<i;++j) mpq_mul(Pq,Pq,MinorP);
    //         //     for (int j=i;j<HPCount;++j) mpq_mul(Pq,Pq,MajorP);
    //         //     mpq_add(PRandom,PRandom,Pq);
    //         //     // fprintf(stderr,"PRandom:%s, Alpha:%s, MajorP:%s, MinorP:%s\n",mpq_get_str(NULL,10,PRandom),mpq_get_str(NULL,10,Alpha),mpq_get_str(NULL,10,MajorP),mpq_get_str(NULL,10,MinorP));
    //         //     // unsigned long long P=1;
    //         //     // for (int j=HPCount-i+1;j<=HPCount;++j) P*=j;
    //         //     // for (int j=1;j<=i;++j) P/=j;
    //         //     // double dP=((double)P)*MultipleP;
    //         //     // PRandomOrHomo+=dP;
    //         // }
    //         // if (MinHP!=0) fprintf(stderr,"PRandom:%s, Alpha:%s, MajorP:%s, MinorP:%s, MinHP:%d, HP:%d\n",mpq_get_str(NULL,10,PRandom),mpq_get_str(NULL,10,Alpha),mpq_get_str(NULL,10,MajorP),mpq_get_str(NULL,10,MinorP),MinHP,HPCount);
    //         // mpq_sub(PRandom,One,PRandom);
    //         // if (MinHP!=0) fprintf(stderr,"PRandom:%s, Alpha:%s, MajorP:%s, MinorP:%s, MinHP:%d, HP:%d, Success:%d\n",mpq_get_str(NULL,10,PRandom),mpq_get_str(NULL,10,Alpha),mpq_get_str(NULL,10,MajorP),mpq_get_str(NULL,10,MinorP),MinHP,HPCount,mpq_cmp(PRandom,SubAlpha)>0);
    //         // if (PRandomOrHomo<Alpha)
    //         // if (mpq_cmp(PRandom,Alpha)<0)
    //         if (mpq_cmp(PRandom,Alpha)<0)
    //         {
    //             Sample["GT"]="1/1";
    //         }
    //     }
    //     // if (HPRatios[1]/(HPRatios[1]+HPRatios[2])> GTHomoRatio)
    //     // {
    //     //     // HPGTKeep=true;
    //     //     Sample["GT"]="0/1";
    //     // }
    //     // else if (HPRatios[2]/(HPRatios[1]+HPRatios[2])>GTHomoRatio)
    //     // {
    //     //     // HPGTKeep=true;
    //     //     Sample["GT"]="0/1";
    //     // }
    //     // else
    //     // {
    //     //     // HPGTKeep=true;
    //     //     Sample["GT"]="1/1";
    //     // }
    // }
    if (SHPRatio>HPRatio)
    {
        if (HPRatios[1]/(HPRatios[1]+HPRatios[2])>HomoRatio || HPRatios[2]/(HPRatios[1]+HPRatios[2])>HomoRatio)
        {
            sort(SignatureCluster.begin(),SignatureCluster.end(),[](Signature &a, Signature&b){
                return a.HP==b.HP? a<b : a.HP<b.HP;
            });
            int End0=0;
            for (;End0<SignatureCluster.size();++End0) if (SignatureCluster[End0].HP!=0) break;
            SignatureCluster.erase(SignatureCluster.begin(),SignatureCluster.begin()+End0);
            // CutOffRatio=HomoCutoffRatio;
            CutOffRatio-=Args.HomoMinus[SVTypeI];
            CutOffRatio-=HomoMinusRatio*(SHPRatio-HPRatio);
            // CutOffRatio=0.9;
            // if (CutOffRatio<HomoCutoffRatio) fprintf(stderr,"CutOffRatio:%lf, HomoCutoffRatio:%lf, SHPRatio:%lf, HPRatio:%lf",CutOffRatio, HomoCutoffRatio, SHPRatio, HPRatio);
        }
        // else CutOffRatio-=0.01*(SHPRatio-HPRatio);
        else
        {
            CutOffRatio-=Args.NonHomoMinus[SVTypeI];
            CutOffRatio-=Args.NonHomoMinusRatio[SVTypeI]*(SHPRatio-HPRatio);
        }
        // else CutOffRatio=10000;
        // resolveHPRecord(HPCounts, TheContig, SignatureCluster, Core, AllPrimarySegments, CoverageWindows, WholeCoverage, Args, CoverageWindowsSums, CheckPoints, CheckPointInterval);
        // return;
    }
    else
    {
        // CutOffRatio+=float(HPCounts[1]+HPCounts[2])/float(HPCounts[0]+HPCounts[1]+HPCounts[2]);
        CutOffRatio+=Args.LowHPPlus[SVTypeI]*(HPRatio-SHPRatio);
        CutOffRatio+=Args.LowHPPlusRatio[SVTypeI]*(HPRatio-SHPRatio);
    }
    statCluster(SignatureCluster,SS,ST,SS2,ST2);
    // if (SVTypeI==0 or SVTypeI==1) if (!cflag) {Keep=false;return;}
    calcM3L(SignatureCluster);
    tuple<int,int,int> Site=analyzeSignatureCluster(SignatureCluster, SVType, Args);
    int MeanSVLen=get<2>(Site);
    SVLen=get<1>(Site);
    Pos=get<0>(Site);//0-bsed now, after ref and alt then transform to 1-based, but should be the base before variantion. End should be the last base, but also should be transform to 1-based. So they don't change.
    //VCF version 4.2 says if alt is <ID>, the pos is the base preceding the polymorphism. No mention of the "base 1" rule.
    
    if (SVLen<=0) {Keep=false;return;}

    int CBegin=0, CEnd=SignatureCluster.size();
    CR=0;
    if (Core.End!=0)
    {
        CBegin=Core.Begin;
        CEnd=Core.End;
        CR=(CEnd-CBegin)/(SignatureCluster.size());
    }

    double LengthSD;
    LS=getLengthSDRatioScore(SignatureCluster,MeanSVLen,&LengthSD);

    if (ST>=Args.MinimumPreciseTemplates && LengthSD<Args.PreciseStandard && calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()))<Args.PreciseStandard && calcSD(EndIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),EndIter<int,vector<Signature>::iterator>(SignatureCluster.end()))<Args.PreciseStandard) Precise=true;
    else Precise=false;

    bool HasLeft=false, HasRight=false;
    if (SVType=="INV")
    {
        for (int i=0;i<SignatureCluster.size();++i)
        {
            HasLeft|=SignatureCluster[i].InvLeft;
            HasRight|=SignatureCluster[i].InvRight;
        }
        INFO="LR=";
        if (HasLeft) INFO+="L";
        if (HasRight) INFO+="R";
    }
    // bool Covered=false;
    // if (SVType=="DUP")
    // {
    //     for (int i=0;i<SignatureCluster.size();++i) if (SignatureCluster[i].Covered) {Covered=true;break;}
    // }
    // if (SVType=="DUP" && !Covered) {Keep=false;return;}
    double Score=ST;
    if (Args.WeightFilter)
    {
        Score=getWeightSum(SignatureCluster);
    }
    if (Args.NoFilter) Keep=true;
    else if (ST<2) {Keep=false;return;}
    else if (HPGTKeep) Keep=true;
    else if (Args.Filter2ST) Keep=true;
    else if (ST2>30) Keep=true;
    else
    {
        if (Args.CalcPosSTD || Args.MinPosSTD[SVTypeI]!=-1 || SVTypeI==1)
        {
            PSD=calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
            if (Args.MinPosSTD[SVTypeI]!=-1 && PSD>Args.MinPosSTD[SVTypeI]) {Keep=false;return;}
        }
        double (*ASSBases)[2]=Args.ASSBases;
        double (*ASSCoverageMulti)[2]=Args.ASSCoverageMulti;
        double (*LSDRSs)[2]=Args.LSDRSs;
        if (SVType=="DUP" or SVType=="INV") if (SS-ST>ST) {Keep=false;return;}
        if (SVType=="INV")
        {
            Keep=true;
            if (!(HasLeft&&HasRight)) {Keep=false;return;}
        }
        // if (Score>=ASSBases[SVTypeI][0]+WholeCoverage*ASSCoverageMulti[SVTypeI][0] && LS>=LSDRSs[SVTypeI][0]) {Keep=true;}
        // else
        // {
        //     if (Score<ASSBases[SVTypeI][1]+WholeCoverage*ASSCoverageMulti[SVTypeI][1]) {Keep=false;return;}
        //     if (LS<LSDRSs[SVTypeI][1]) {Keep=false;return;}
        // }
        if (Score>=(ASSBases[SVTypeI][0]+WholeCoverage*ASSCoverageMulti[SVTypeI][0])*CutOffRatio && LS>=LSDRSs[SVTypeI][0]*CutOffRatio) {Keep=true;}
        else
        {
            if (Score<(ASSBases[SVTypeI][1]+WholeCoverage*ASSCoverageMulti[SVTypeI][1])*CutOffRatio) {Keep=false;return;}
            if (LS<LSDRSs[SVTypeI][1]*CutOffRatio) {Keep=false;return;}
        }
    }
    if (SVType=="INS") InsConsensus=getInsConsensus(SVLen,SignatureCluster);

    if (!HPGTKeep && Sample.count("GT")==0) genotype(TheContig,AllPrimarySegments,CoverageWindows,CoverageWindowsSums,CheckPoints,CheckPointInterval,Args);
    ID=".";
    QUAL=".";
    FILTER="PASS";
    CHROM=TheContig.Name;
    if (SVTypeI==2)
    {
        // double CCN=0;
        // for (int i=0;i<SignatureCluster.size();++i)
        // {
        //     CCN+=SignatureCluster[i].CN;
        // }
        // CCN/=(double)(SignatureCluster.size());
        // CN=int(CCN+0.5)
        for (int i=0;i<SignatureCluster.size();++i)
        {
            CN=max(CN, SignatureCluster[i].CN);
        }
        // CN=int((double)(SS)/(double)(ST)+0.5);
    }
    #ifdef DEBUG
    MergeStrings="";
    if (SVType=="DEL")
    {
        int MSMaxOut=3;
        int MSCount=0;
        for (int i=0;i<SignatureCluster.size();++i)
        {
            if (SignatureCluster[i].MergeString!="")
            {
                if (MSCount>=MSMaxOut)
                {
                    MergeStrings+="...";
                    break;
                }
                string MSub=SignatureCluster[i].MergeString.substr(SignatureCluster[i].MergeString.find('['),SignatureCluster[i].MergeString.find(']'));
                int MCount=count(MSub.begin(),MSub.end(),'D');
                MergeStrings+=(MCount==1?"":"*")+SignatureCluster[i].MergeString+":";
                ++MSCount;
            }
        }
    }
    #endif
}
void VCFRecord::resolveHPRecord(int * HPCounts, const Contig & TheContig,vector<Signature> & SignatureCluster, ClusterCore &Core, SegmentSet & AllPrimarySegments, float* CoverageWindows, double WholeCoverage, Arguments& Args, float * CoverageWindowsSums, float * CheckPoints, int CheckPointInterval)
{
    float HomoRatio=0.8f, HomoCutoffRatio=0.5f;
    float HPRatios[3]={0,0,0};
    for (int i=0;i<3;++i) HPRatios[i]=float(HPCounts[i])/float(SignatureCluster.size());
    sort(SignatureCluster.begin(),SignatureCluster.end(),[](Signature &a, Signature&b){
        return a.HP==b.HP? a<b : a.HP<b.HP;
    });
    int End0=0;
    for (;End0<SignatureCluster.size();++End0) if (SignatureCluster[End0].HP!=0) break;
    SignatureCluster.erase(SignatureCluster.begin(),SignatureCluster.begin()+End0);
    if (HPRatios[1]/(HPRatios[1]+HPRatios[2])>HomoRatio || HPRatios[2]/(HPRatios[1]+HPRatios[2])>HomoRatio)
    {
        statCluster(SignatureCluster,SS,ST,SS2,ST2);
        SVType=getSVType(SignatureCluster);
        SVTypeI=SignatureCluster[0].SupportedSV;
        // if (SVTypeI==0 or SVTypeI==1) if (!cflag) {Keep=false;return;}
        calcM3L(SignatureCluster,true);
        tuple<int,int,int> Site=analyzeSignatureCluster(SignatureCluster, SVType, Args);
        int MeanSVLen=get<2>(Site);
        SVLen=get<1>(Site);
        Pos=get<0>(Site);//0-bsed now, after ref and alt then transform to 1-based, but should be the base before variantion. End should be the last base, but also should be transform to 1-based. So they don't change.
        //VCF version 4.2 says if alt is <ID>, the pos is the base preceding the polymorphism. No mention of the "base 1" rule.
        
        if (SVLen<=0) {Keep=false;return;}

        int CBegin=0, CEnd=SignatureCluster.size();
        CR=0;
        if (Core.End!=0)
        {
            CBegin=Core.Begin;
            CEnd=Core.End;
            CR=(CEnd-CBegin)/(SignatureCluster.size());
        }

        double LengthSD;
        LS=getLengthSDRatioScore(SignatureCluster,MeanSVLen,&LengthSD);

        if (ST>=Args.MinimumPreciseTemplates && LengthSD<Args.PreciseStandard && calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()))<Args.PreciseStandard && calcSD(EndIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),EndIter<int,vector<Signature>::iterator>(SignatureCluster.end()))<Args.PreciseStandard) Precise=true;
        else Precise=false;

        bool HasLeft=false, HasRight=false;
        if (SVType=="INV")
        {
            for (int i=0;i<SignatureCluster.size();++i)
            {
                HasLeft|=SignatureCluster[i].InvLeft;
                HasRight|=SignatureCluster[i].InvRight;
            }
            INFO="LR=";
            if (HasLeft) INFO+="L";
            if (HasRight) INFO+="R";
        }
        double Score=ST;
        if (Args.WeightFilter)
        {
            Score=getWeightSum(SignatureCluster);
        }
        if (Args.NoFilter) Keep=true;
        else if (ST<2) {Keep=false;return;}
        else if (Args.Filter2ST) Keep=true;
        else if (ST2>30) Keep=true;
        else
        {
            if (Args.CalcPosSTD || Args.MinPosSTD[SVTypeI]!=-1 || SVTypeI==1)
            {
                PSD=calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
                if (Args.MinPosSTD[SVTypeI]!=-1 && PSD>Args.MinPosSTD[SVTypeI]) {Keep=false;return;}
            }
            double (*ASSBases)[2]=Args.ASSBases;
            double (*ASSCoverageMulti)[2]=Args.ASSCoverageMulti;
            double (*LSDRSs)[2]=Args.LSDRSs;
            if (SVType=="DUP" or SVType=="INV") if (SS-ST>ST) {Keep=false;return;}
            if (SVType=="INV")
            {
                Keep=true;
                if (!(HasLeft&&HasRight)) {Keep=false;return;}
            }
            if (Score>=(ASSBases[SVTypeI][0]+WholeCoverage*ASSCoverageMulti[SVTypeI][0])*HomoCutoffRatio && LS>=LSDRSs[SVTypeI][0]*HomoCutoffRatio) {Keep=true;}
            else
            {
                if (Score<(ASSBases[SVTypeI][1]+WholeCoverage*ASSCoverageMulti[SVTypeI][1])*HomoCutoffRatio) {Keep=false;return;}
                if (LS<LSDRSs[SVTypeI][1]*HomoCutoffRatio) {Keep=false;return;}
            }
        }
        if (SVType=="INS") InsConsensus=getInsConsensus(SVLen,SignatureCluster);

        // Sample["GT"]="0/1";
        genotype(TheContig,AllPrimarySegments,CoverageWindows,CoverageWindowsSums,CheckPoints,CheckPointInterval,Args);
        ID=".";
        QUAL=".";
        FILTER="PASS";
        CHROM=TheContig.Name;
    }
}

void VCFRecord::resolveRef(const Contig & TheContig, faidx_t * Ref, unsigned TypeCount, double CC, Arguments & Args)
{
    ID="kled."+Args.RunHash.substr(0,8)+"."+SVType+"."+to_string(TypeCount);
    int TLen;
    int End;
    bool OutTag=false;
    if (SVLen>1000000) OutTag=true;
    //According to VCF v4.2, Pos should be 1 base before the event, unless event happens at the 1st pos, then ref should contain the base after the event
    //End should contain the last base, and either pos is 1 base before the event or ref contains 1 base after the event, End is always Pos+SVLen (excpet for SVs like INS)
    if (Pos!=0) --Pos;
    if (SVType=="INS") End=Pos;
    else End=Pos+SVLen;//End should contain the last base, but 
    char * TSeq=faidx_fetch_seq(Ref,TheContig.Name.c_str(),Pos,End,&TLen);
    if (OutTag)
    {
        REF=TSeq[0];
        // if (Pos==0) ALT="<"+SVType+">"+REF;
        ALT="<"+SVType+">";//as described in VCF v4.2 example, no need to add ref
    }
    else
    {
        if (SVType=="DUP" || SVType=="INV") REF=TSeq[0];
        else REF=TSeq;
        if (SVType=="DEL")
        {
            // if (Pos==0) ALT="<DEL>"+REF[REF.length()-1];
            ALT=REF[0];
        }
        else if (SVType=="INS")
        {
            if (InsConsensus=="") ALT="<INS>";
            else ALT=REF[0]+InsConsensus;
        }
        else
        {
            ALT="<"+SVType+">";
        }
    }
    free(TSeq);
    ++Pos;++End;//trans to 1-based
    if (INFO!="") INFO+=";";
    INFO+=(Precise?"PRECISE;":"IMPRECISE;")+string("SVTYPE=")+SVType+";END="+to_string(End)+";SVLEN="+to_string(SVType=="DEL"?-SVLen:SVLen)+(SVType=="DUP"?";CN="+to_string(CN):"")+";SS="+to_string(SS)+";ST="+to_string(ST)+";LS="+to_string(LS)+";CV="+to_string(CV)+";SS2="+to_string(SS2)+";ST2="+to_string(ST2)+";CC="+to_string(CC)+";CR="+to_string(CR)+";M3L="+to_string(MinLength)+","+to_string(MediumLength)+","+to_string(MaxLength)+";PSTD="+to_string(PSD)+";HPC="+to_string(HPCounts[0])+","+to_string(HPCounts[1])+","+to_string(HPCounts[2])DEBUG_CODE(+(MergeStrings==""?"":(";MSs="+MergeStrings)));
}

VCFRecord::operator std::string() const
{
    string Output=CHROM+"\t"+to_string(Pos)+"\t"+ID+"\t"+REF+"\t"+ALT+"\t"+QUAL+"\t"+FILTER+"\t"+INFO+"\t";
    bool FirstFormat=true;
    string SampleString=("");
    for (auto i=Sample.begin();i!=Sample.end();++i)
    {
        if (!FirstFormat) Output+=":";
        Output+=i->first;
        SampleString+=i->second;
        FirstFormat=false;
    }
    Output+="\t"+SampleString;
    return Output;
}

bool VCFRecord::operator<(const VCFRecord& Other) const
{
    return Pos<Other.Pos;
}

HeaderEntry::HeaderEntry(string AClass, string AID, string ADescription, string ANumber, string AType, string ASource, string AVersion)
:Class(AClass),ID(AID),Description(ADescription),Type(AType),Number(ANumber),Source(ASource),Version(AVersion)
{
}
HeaderEntry::operator string() const
{
    if (Class=="INFO" || Class=="FORMAT")
        return "##"+Class+"=<ID="+ID+",Number="+Number+",Type="+Type+",Description=\""+Description+"\""+(Source!=""?(",Source=\""+Source+"\""):"")+(Version!=""?",Version=\""+Version+"\"":"")+">";
    else
        return "##"+Class+"=<ID="+ID+",Description=\""+Description+"\">";
}

inline string getDate()
{
    time_t t = time(0);   // get time now
    tm* now = localtime(&t);
    return to_string(now->tm_year + 1900)+(now->tm_mon<9?"0":"")+to_string((now->tm_mon + 1))+to_string(now->tm_mday);
}

VCFHeader::VCFHeader(const char * AReference)
{
    FileFormat="##fileformat=VCFv4.2";
    FileDate="##fileDate="+getDate();
    Reference=string("##reference=")+string(AReference);
}
void VCFHeader::addHeaderEntry(const HeaderEntry & Entry)
{
    HeaderEntries.push_back(Entry);
}
void VCFHeader::addSample(const char * SampleName)
{
    SampleNames.push_back(SampleName);
}

void VCFHeader::addContig(const Contig & TheContig)
{
    Contigs.push_back(TheContig);
}

std::string VCFHeader::genHeader(const Arguments & Args)
{
    string Output;
    Output+=FileFormat+"\n"+FileDate+"\n"+Reference;
    for (int i=0;i<Contigs.size();++i)
    {
        Output+=string("\n")+"##contig=<ID="+Contigs[i].Name+",length="+to_string(Contigs[i].Size)+">";
    }
    for (int i=0;i<HeaderEntries.size();++i)
    {
        Output+=string("\n")+string(HeaderEntries[i]);
    }
    Output+="\n##SVCaller=\"Kled version "+string(Args.Version)+"\"";
    Output+="\n##CommandLine=\""+Args.CommandLine+"\"";
    Output+="\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (int i=0;i<SampleNames.size();++i) Output+="\t"+SampleNames[i];
    return Output;
}

void addKledEntries(VCFHeader & Header)
{
    Header.addHeaderEntry(HeaderEntry("INFO","SVTYPE","Type of structural variant","1","String"));
    Header.addHeaderEntry(HeaderEntry("INFO","END","End position of the variant described in this record","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","SVLEN","Difference in length between REF and ALT alleles","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","PRECISE","Precise structural variant","0","Flag"));
    Header.addHeaderEntry(HeaderEntry("INFO","IMPRECISE","Imprecise structural variant","0","Flag"));
    Header.addHeaderEntry(HeaderEntry("INFO","SS","Number of supported signatures","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","ST","Number of supported templates","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","SS2","Number of supported NGS signatures","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","ST2","Number of supported NGS templates","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","LS","Length SD ratio score (100-ratio*100) of the cluster.","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","LR","L,R,LR","1","String"));
    Header.addHeaderEntry(HeaderEntry("INFO","CV","Nearby coverage for genotyping.","1","Float"));
    Header.addHeaderEntry(HeaderEntry("INFO","CC","Average coverage of this contig.","1","Float"));
    Header.addHeaderEntry(HeaderEntry("INFO","CN","Copy number for duplication alleles(normal=1).","1","Float"));//As defined in VCFv4.4, the info copy number shall be the CN for the allele
    Header.addHeaderEntry(HeaderEntry("INFO","CR","Core ratio.","1","Float"));
    Header.addHeaderEntry(HeaderEntry("INFO","M3L","(min length, medium length, max length).","3","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","PSTD","Position STD, -1 if not caclulated.","1","Float"));
    Header.addHeaderEntry(HeaderEntry("INFO","HPC","HP counts:HP0,HP1,HP2.","3","Integer"));
    DEBUG_CODE(Header.addHeaderEntry(HeaderEntry("INFO","MSs","MergeStrings.","1","String"));)
    // Header.addHeaderEntry(HeaderEntry("INFO","CS","Nearby coverage for genotyping.(by windows)","1","Float"));
    Header.addHeaderEntry(HeaderEntry("ALT","DEL","Deletion"));
    Header.addHeaderEntry(HeaderEntry("ALT","DUP","Duplication"));
    Header.addHeaderEntry(HeaderEntry("FORMAT","GT","Genotype","1","String"));
    // Header.addHeaderEntry(HeaderEntry("INFO","SCORES","Scores, supported reads, supported signatures, (begin sd+end sd)/2, length sd, number of signature source","6","Integer"));
}

void HPClustersDistinction(vector<Signature> &Cluster, vector<vector<Signature>> &HPClusters, Arguments& Args)
{
    HPClusters[1].clear();
    HPClusters[2].clear();
    double HP1Mean=0, HP2Mean=0;
    int H1Size=0,H2Size=0;
    double HP1SD=0, HP2SD=0;
    
    int Size=0;
    double Mean=0;
    for (int i=0;i<Cluster.size();++i)
    {
        if (Cluster[i].HP==1)
        {
            HP1Mean+=Cluster[i].Length;
            ++H1Size;
        }
        if (Cluster[i].HP==2)
        {
            HP2Mean+=Cluster[i].Length;
            ++H2Size;
        }
    }
    HP1Mean/=(double)H1Size;
    HP2Mean/=(double)H2Size;
    double SD=0;
    for (int i=0;i<Cluster.size();++i)
    {
        if (Cluster[i].HP==1)
        {
            HP1SD+=pow(((double)Cluster[i].Length-HP1Mean),2.0);
        }
        if (Cluster[i].HP==2)
        {
            HP2SD+=pow(((double)Cluster[i].Length-HP1Mean),2.0);
        }
    }
    HP1SD/=(double)H1Size;
    HP2SD/=(double)H2Size;
    double MaxSD=max(HP1SD,HP2SD);
    if (abs(HP1Mean-HP2Mean)>3.0*MaxSD)
    {
        for (int i=0;i<Cluster.size();++i)
        {
            if (Cluster[i].HP==0)
            {
                HPClusters[1].push_back(Cluster[i]);
                HPClusters[2].push_back(Cluster[i]);
            }
            else
            {
                HPClusters[Cluster[i].HP].push_back(Cluster[i]);
            }
        }
    }
}

void VCFRecord::hapGT(unsigned SmallHapCount, unsigned BigHapCount)
{
    double SHPRatio=float(HPCounts[1]+HPCounts[2])/float(HPCounts[0]+HPCounts[1]+HPCounts[2]);
    double GTHPRatio=0.8;
    if (SHPRatio>GTHPRatio)
    {
        mpq_t Alpha, MaxAlpha, MinAlpha;
        mpq_init(Alpha);
        mpq_init(MaxAlpha);
        mpq_init(MinAlpha);
        mpq_set_str(Alpha,"2/100000",10);
        mpq_set_str(MaxAlpha,"5/100000",10);
        mpq_set_str(MinAlpha,"1/1000000",10);
        int HPCount=HPCounts[1]+HPCounts[2];
        int MinHP=min(HPCounts[1],HPCounts[2]);
        mpq_t MinorP,One,DefaultMinorP,MPRatio;
        mpq_init(MinorP);
        mpq_init(MPRatio);
        mpq_init(One);
        mpq_init(DefaultMinorP);
        mpq_set_ui(One,1,1);
        mpq_set_str(MinorP,"15/1000",10);
        mpq_set_str(DefaultMinorP,"15/1000",10);
        if (BigHapCount>100)
        {
            mpq_set_ui(MinorP,SmallHapCount,BigHapCount+SmallHapCount);
            // mpq_div(MPRatio,DefaultMinorP,MinorP);
            // mpq_div(MPRatio,MinorP,DefaultMinorP);
            // mpq_mul(Alpha,Alpha,MPRatio);
            // if (mpq_cmp(Alpha,MaxAlpha)>0) mpq_set(Alpha,MaxAlpha);
            // else if (mpq_cmp(Alpha,MinAlpha)<0) mpq_set(Alpha,MinAlpha);
        }
        mpq_t PRandom;
        mpq_init(PRandom);
        mpq_set_ui(PRandom,0,1);
        biCdf(PRandom,HPCount,MinorP,MinHP,HPCount);
        if (mpq_cmp(PRandom,Alpha)<0)
        {
            if (Sample["GT"]!="1/1") ConcurrentGT=2;
            else ConcurrentGT=1;
            Sample["GT"]="1/1";
        }
        // else
        // {
        //     if (ConcurrentGT==2)//If changed last round, change back
        //     {
        //         ConcurrentGT=1;
        //         Sample["GT"]="0/1";
        //     }
        // }
    }
}