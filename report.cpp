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

tuple<int,int> getConsensusSite(vector<Signature> &Signatures)
{
    if (Signatures.size()==0) return tuple<int,int>(-1,-1);
    sort(Signatures.data(),Signatures.data()+Signatures.size());
    return consensusFixed(Signatures,Signatures.size()/2,stdStat);
}

tuple<int,int> analyzeSignatureCluster(vector<Signature> &SignatureCluster)
{
    float L2Weight=0.75;
    int Pos,Length;
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
    return tuple<int,int>(Pos,Length);
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
    // Scores.push_back(STS);
    // Scores.push_back(SSS);
    double BESD=0.0;//calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    // BESD+=calcSD(EndIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),EndIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    // BESD/=2.0;
    double LengthSD=calcSD(LengthIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),LengthIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    double BESDS=0.0;//BESD>50?100.0:double(BESD)/50.0*100.0;BESDS=100.0-BESDS;
    double LengthSDS=0.0;//LengthSD>50?100.0:double(LengthSD)/50.0*100.0;LengthSDS=100.0-LengthSDS;
    Scores.push_back(BESDS);
    // double BESDRatio=BESD/(double)SVLen;
    double BESDRatioS=0.0;//BESDRatio>0.5?100:BESDRatio/0.5*100.0;BESDRatioS=100.0-BESDRatioS;
    Scores.push_back(BESDRatioS);
    Scores.push_back(LengthSDS);
    double LengthSDRatio=LengthSD/(double)SVLen;
    double LengthSDRatioS=LengthSDRatio>0.5?100:LengthSDRatio/0.5*100.0;LengthSDRatioS=100.0-LengthSDRatioS;
    Scores.push_back(LengthSDRatioS);
    // set<int> Sources;
    // for (int i =0;i<SignatureCluster.size();++i)
    // {
    //     if (SignatureCluster[i].Type!=-1) Sources.insert(SignatureCluster[i].Type);
    // }
    Scores.push_back(0.0);//double(Sources.size())/3.0*100.0);
    return Scores;
}

int getLengthSDRatioScore(vector<Signature> &SignatureCluster,int SVLen, double * SD=NULL)
{
    double LengthSD=calcSD(LengthIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),LengthIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    if (SD!=NULL) *SD=LengthSD;
    return MAX(0.0,100.0-(LengthSD/(double)SVLen*100.0));
}

double getAverageCoverage(int Begin, int End, double * CoverageWindows, Arguments & Args, double* CoverageWindowsSums, double* CheckPoints, int CheckPointInterval)
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
            Cov*=((double)(i-WBegin))/(i-WBegin+1);//in case too big
            Cov+=CoverageWindows[i]*1/(i-WBegin+1);
        }
    }
    return Cov;
}

bool keepCluster(vector<Signature> &SignatureCluster, int & SS, int &ST)
{
    SS=0;
    set<string> SupportTemps;;
    for (int i =0;i<SignatureCluster.size();++i)
    {
        ++SS;
        SupportTemps.insert(SignatureCluster[i].TemplateName);
    }
    ST=SupportTemps.size();
    return true;
    // if (ST>=1) return true;
    // return false;
    // double LengthSD=calcSD(LengthIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),LengthIter<int,vector<Signature>::iterator>(SignatureCluster.end()));
    // // if (ST>=12 && LengthSD<20) return true;
    // if (ST>=10 && LengthSD<10) return true;
    // if (ST>=5 && LengthSD<5 ) return true;
    // set<int> Sources;
    // for (int i =0;i<SignatureCluster.size();++i)
    // {
    //     if (SignatureCluster[i].Type!=-1) Sources.insert(SignatureCluster[i].Type);
    // }
    // // if (ST>=3 && LengthSD<5 && Sources.size()>1) return true;
    // return false;
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

string VCFRecord::genotype(const Contig & TheContig, SegmentSet & AllPrimarySegments, double * CoverageWindows, double *CoverageWindowsSums, double* Checkpoints, int CheckPointInterval, Arguments & Args)
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
    // else All=getAverageCoverage(Pos,End,CoverageWindows,Args, CoverageWindowsSums, Checkpoints, CheckPointInterval);
    // else
    // {
    //     Pos=MAX(Pos-Scope,0);
    //     End+=Scope;
    // }
    // double ErrorRate=0.5;//Error rate of support/deny wrongly seqed/mapped to deny/support.
    // double EAF=0.5;//Estimated Allele Frequencey.
    CV=All;
    if (double(ST)/All<HomoThreshold) Sample["GT"]="0/1";
    else Sample["GT"]="1/1";
    return Sample["GT"];
    // if (All<=ST) return "1/1";
    // double AF=0.5;

    // double lambda=All/2.0;
    // double Slack=1.6;
    // if (SVType=="INS")
    // {
    //     AF=0.57;
    //     Slack=0.8;
    // }
    // int DT=All-int(ST);
    // if (DT<0) DT=0;

    // double GL[2]={0,0};//0/1,1/1
    // double post=2.0*AF*(1.0-AF)+AF*AF;
    // GL[0]=2.0*AF*(1.0-AF)/post*poisson(lambda,DT)*poisson(lambda,ST);
    // GL[1]=AF*AF/post*pow(poisson(lambda,ST/2.0),2);
    // if (GL[1]>Slack*GL[0])
    // {
    //     Sample["GT"]="1/1";
    // }
    // else Sample["GT"]="0/1";
    // return Sample["GT"];
    // // double GL[3]={0,0,0};//0/0,0/1,1/1
    // // // GL[0]=(1.0-AF)*(1.0-AF)*pow(normal(lambda,pow(lambda,0.5),DT/2.0),2);
    // // // GL[1]=2.0*AF*(1.0-AF)*normal(lambda,pow(lambda,0.5),DT)*normal(lambda,pow(lambda,0.5),ST);
    // // // GL[2]=AF*AF*pow(normal(lambda,pow(lambda,0.5),ST/2.0),2);
    // // GL[0]=pow(poisson(lambda,DT/2.0),2);
    // // GL[1]=poisson(lambda,DT)*poisson(lambda,ST);
    // // GL[2]=pow(poisson(lambda,ST/2.0),2);
    // // if (GL[2]>Slack*GL[1]) return "1/1";
    // // return "0/1";

    // // fprintf(stderr, "%d,%d,",DT, int(ST));
    // return calc_GL_cute(DT,ST);
    // double PrioPs[3];
    // PrioPs[0]=pow((1.0-EAF),2)*pow(ErrorRate,ST);
    // PrioPs[1]=(1.0-EAF)*EAF*2*combination(All, DT)*pow(0.5,All);
    // PrioPs[2]=EAF*EAF*pow(ErrorRate,DT);
    // double PrioP=0;
    // for (int i=0;i<3;++i) PrioP+=PrioPs[i];
    // double PostPs[3];
    // double MaxP=0;
    // for (int i=0;i<3;++i) {PostPs[i]=PrioPs[i]/PrioP;if (MaxP<PostPs[i]) MaxP=PostPs[i];}
    // if (MaxP==PostPs[3]) return "1/1";
    // return "0/1";
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
    // for (int i=0;i<SignatureCluster.size();++i)
    // {
    //     if (abs(int(SignatureCluster[i].InsBases.length())-SVLen)<=SVLen*Endurance) ++ConsensusCount;
    // }
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

int VN=0;

VCFRecord::VCFRecord(const Contig & TheContig, faidx_t * Ref,vector<Signature> & SignatureCluster, SegmentSet & AllPrimarySegments, double* CoverageWindows, double WholeCoverage, Arguments& Args, double * CoverageWindowsSums, double * CheckPoints, int CheckPointInterval)
{
    assert(SignatureCluster.size()>=0);
    resizeCluster(SignatureCluster,Args.MaxClusterSize);
    if (keepCluster(SignatureCluster,SS,ST)) Keep=true;
    else {Keep=false; return;}
    SVType=getSVType(SignatureCluster);
    SVTypeI=SignatureCluster[0].SupportedSV;
    tuple<int,int> Site=analyzeSignatureCluster(SignatureCluster);
    SVLen=get<1>(Site);
    Pos=get<0>(Site);//0-bsed now, after ref and alt then transform to 1-based, but should be the base before variantion. End should be the last base, but also should be transform to 1-based. So they don't change.
    //VCF version 4.2 says if alt is <ID>, the pos is the base preceding the polymorphism. No mention of the "base 1" rule.

    unsigned long long llPos=0;
    unsigned long long llSVLen=0;
    for (int i=0;i<SignatureCluster.size();++i)
    {
        llPos+=SignatureCluster[i].Begin;
        llSVLen+=SignatureCluster[i].Length;
    }
    Pos=llPos/SignatureCluster.size();
    SVLen=llSVLen/SignatureCluster.size();
    
    if (SVLen<=0) {Keep=false;return;}

    double LengthSD;
    LS=getLengthSDRatioScore(SignatureCluster,SVLen,&LengthSD);

    if (ST>=Args.MinimumPreciseTemplates && LengthSD<Args.PreciseStandard && calcSD(BeginIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),BeginIter<int,vector<Signature>::iterator>(SignatureCluster.end()))<Args.PreciseStandard && calcSD(EndIter<int,vector<Signature>::iterator>(SignatureCluster.begin()),EndIter<int,vector<Signature>::iterator>(SignatureCluster.end()))<Args.PreciseStandard) Precise=true;
    else Precise=false;

    if (Args.NoFilter) Keep=true;
    else
    {
        double (*ASSBases)[2]=Args.ASSBases;
        double (*ASSCoverageMulti)[2]=Args.ASSCoverageMulti;
        double (*LSDRSs)[2]=Args.LSDRSs;
        if (Args.AllCLR)
        {
            ASSBases=Args.CLRASSBases;
            ASSCoverageMulti=Args.CLRASSCoverageMulti;
            LSDRSs=Args.CLRLSDRSs;
        }
        else if (Args.AllCCS)
        {
            ASSBases=Args.CCSASSBases;
            ASSCoverageMulti=Args.CCSASSCoverageMulti;
            LSDRSs=Args.CCSLSDRSs;
        }
        if (SVType=="DUP" or SVType=="INV") if (SS-ST>ST) {Keep=false;return;}
        if (SVType=="INV")
        {
            Keep=true;
            bool HasLeft=false, HasRight=false;
            for (int i=0;i<SignatureCluster.size();++i)
            {
                HasLeft|=SignatureCluster[i].InvLeft;
                HasRight|=SignatureCluster[i].InvRight;
            }
            INFO="LR=";
            if (HasLeft) INFO+="L";
            if (HasRight) INFO+="R";
            if (HasLeft&&HasRight)
            {
                if (ST<ASSBases[SVTypeI][0]+WholeCoverage*ASSCoverageMulti[SVTypeI][0]) {Keep=false;return;}
                if (LS<LSDRSs[SVTypeI][0]) {Keep=false;return;}
            }
            else
            {
                if (ST<ASSBases[SVTypeI][1]+WholeCoverage*ASSCoverageMulti[SVTypeI][1]) {Keep=false;return;}
                if (LS<LSDRSs[SVTypeI][1]) {Keep=false;return;}
            }
        }
        else if (SS+ST>=ASSBases[SVTypeI][0]+WholeCoverage*ASSCoverageMulti[SVTypeI][0] && LS>=LSDRSs[SVTypeI][0]) {Keep=true;}
        else
        {
            if (SS+ST<ASSBases[SVTypeI][1]+WholeCoverage*ASSCoverageMulti[SVTypeI][1]) {Keep=false;return;}
            if (LS<LSDRSs[SVTypeI][1]) {Keep=false;return;}
            // if (abs(SVLen)>10000)
            // {
            //     double BeforeCovergae=-1;
            //     if (Pos>10000) BeforeCovergae=getAverageCoverage(Pos-10000,Pos,CoverageWindows,Args);
            //     double AfterCoverage=-1;
            //     if (Pos+abs(SVLen)<TheContig.Size-10000) AfterCoverage=getAverageCoverage(Pos+abs(SVLen),Pos+abs(SVLen)+10000,CoverageWindows,Args);
            //     double Ambient=0;
            //     double Valid=0;
            //     if (BeforeCovergae!=-1) {Valid+=1;Ambient+=BeforeCovergae;}
            //     if (AfterCoverage!=-1) {Valid+=1;Ambient+=AfterCoverage;}
            //     if (Valid!=0) Ambient/=Valid;
            //     if (Ambient!=0)
            //     {
            //         double Center=getAverageCoverage(Pos+abs(SVLen)/4.0,Pos+abs(SVLen)*0.75,CoverageWindows,Args);
            //         if (Center>Ambient*0.6) {Keep=false;return;}
            //     }
            // }
            // if (abs(SVLen)>5000 && AmbientCoverage>WholeCoverage*1.0) {Keep=false;return;}
            // else if (Scores[0]>=10 && (Scores[2]+Scores[3]+Scores[4])>=240) Keep=true;
            // else if (Scores[0]>=6 && ((Scores[2]+Scores[3])>=190 || (Scores[2]+Scores[4])>=190 || (Scores[3]+Scores[4])>=190 )) Keep=true;
            // else {Keep=false;return;}
        }
    }
    if (SVType=="INS") InsConsensus=getInsConsensus(SVLen,SignatureCluster);

    genotype(TheContig,AllPrimarySegments,CoverageWindows,CoverageWindowsSums,CheckPoints,CheckPointInterval,Args);
   
    ID=".";
    QUAL=".";
    FILTER="PASS";
    CHROM=TheContig.Name;
    // Cluster=SignatureCluster;
}

void VCFRecord::resolveRef(const Contig & TheContig, faidx_t * Ref, unsigned TypeCount, Arguments & Args)
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
        // if (Pos==0) REF=TSeq[End-Pos];
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
    INFO+=(Precise?"PRECISE;":"IMPRECISE;")+string("SVTYPE=")+SVType+";END="+to_string(End)+";SVLEN="+to_string(SVType=="DEL"?-SVLen:SVLen)+";SS="+to_string(SS)+";ST="+to_string(ST)+";LS="+to_string(LS)+";CV="+to_string(CV);
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
    Header.addHeaderEntry(HeaderEntry("INFO","LS","Length SD ratio score (100-ratio*100) of the cluster.","1","Integer"));
    Header.addHeaderEntry(HeaderEntry("INFO","LR","L,R,LR","1","String"));
    Header.addHeaderEntry(HeaderEntry("INFO","CV","Nearby coverage for genotyping.","1","Float"));
    // Header.addHeaderEntry(HeaderEntry("INFO","CS","Nearby coverage for genotyping.(by windows)","1","Float"));
    Header.addHeaderEntry(HeaderEntry("ALT","DEL","Deletion"));
    Header.addHeaderEntry(HeaderEntry("ALT","DUP","Duplication"));
    Header.addHeaderEntry(HeaderEntry("FORMAT","GT","Genotype","1","String"));
    // Header.addHeaderEntry(HeaderEntry("INFO","SCORES","Scores, supported reads, supported signatures, (begin sd+end sd)/2, length sd, number of signature source","6","Integer"));
}