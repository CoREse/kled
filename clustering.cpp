#include "clustering.h"
#include <math.h>
#include "crelib/crelib.h"
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <list>
#include <set>
#include <omp.h>
using namespace std;
using namespace cre;

float calcOverlap(float B1, float E1, float B2, float E2)
{
    if (B1<=E2 && B1>=B2)
        return MIN(E2-B1,E1-B1);
    if (E1>=B2 && E1<=E2)
        return E1-B2;
    if (B2>=B1 && B2<=E1)//1 covers 2
        return E2-B2;
    else
        return -1;
}

float distance(const Signature &A, const Signature &B, bool Partial, float * PPD, Stats BamStats)
{
    float SDW, PDW, ODW;
    float PDN=900;//TODO: should be realated by tech and type, should make pd>0.5 means not the same
    int BP=bestPrecision(A,B),WP=worstPrecision(A,B);
    if (WP==2) PDN=30;
    else if (WP==0 && BP!=1) PDN=2*(BamStats.Mean+3*BamStats.SD);//should be related with the drp stats
    else if (BP==1 && WP==0) PDN=MAX(PDN,2*(BamStats.Mean+3*BamStats.SD));
    SDW=1;
    PDW=1;
    ODW=1;
    float SD,PD,OD;
    PD=MIN(abs(A.Begin-B.Begin),abs(A.End-B.End));
    PD=MIN(PD,abs((A.Begin+A.End)/2-(B.Begin+B.End)/2));
    PD/=PDN;
    float Overlap=calcOverlap(A.Begin,A.End,B.Begin,B.End);
    if (!Partial || precisionLevel(A)==precisionLevel(B))
    {
        SD=abs(A.Length-B.Length)/MAX(A.Length,B.Length);
        OD=Overlap/(MAX((A.End-A.Begin),(B.End-B.Begin)));
        OD=1-OD;
    }
    else
    {
        const Signature * HighPre, *LowPre;
        if (precisionLevel(A)>precisionLevel(B)) {HighPre=&A;LowPre=&B;}
        else {HighPre=&B;LowPre=&A;}
        SD=abs(HighPre->Length-LowPre->Length)/(HighPre->Length);
        OD=Overlap/(HighPre->End-HighPre->Begin);
        OD=abs(1-OD);
    }
    if (PPD!=NULL) *PPD=PD;
    if (SD>1) SD=10;
    if (PD>1) PD=10;
    if (OD>1) OD=10;
    return (SD*SDW+PDW*PD+ODW*OD);//+2.0-float(BP);
}

inline int first0(short * A, int S, int B=0)
{
    for (int i=B;i<S;++i) if (A==0) return i;
}

bool SigLengthLess(const Signature & a, const Signature &b)
{
    if (a.Length<b.Length) return true;
    if (a.Length>b.Length) return false;
    if (a.Begin<b.Begin) return true;
    if (a.Begin>b.Begin) return false;
    return a.TemplateName<b.TemplateName;
}

void keepLongestPerRead(vector<Signature> & SignatureCluster)
{
    unordered_map<string, Signature> TempsSig;
    for (int i=0;i<SignatureCluster.size();++i)
    {
        if(TempsSig.count(SignatureCluster[i].TemplateName)==0)
        {
            TempsSig[SignatureCluster[i].TemplateName]=SignatureCluster[i];
        }
        else if (TempsSig[SignatureCluster[i].TemplateName].Length<SignatureCluster[i].Length)
        {
            TempsSig[SignatureCluster[i].TemplateName]=SignatureCluster[i];
        }
    }
    SignatureCluster.clear();
    for (auto it=TempsSig.begin();it!=TempsSig.end();++it)
    {
        SignatureCluster.push_back(it->second);
    }
}

void simpleClustering(vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, Stats BamStats, bool CuteVer=false)//like jcrd and cuteSV, SortedSignatures may have deleted ones marked by Type=-1
{
    if (CuteVer)
    {
        float MaxDis=200;
        int MinSupport=0;
        for (int i=0;i<SortedSignatures.size();++i)
        {
            if (SortedSignatures[i].Type==-1) continue;
            if (Clusters.size()==0)
            {
                Clusters.push_back(vector<Signature>());
            }
            if (Clusters.back().size()==0) Clusters.back().push_back(SortedSignatures[i]);
            else
            {
                if (SortedSignatures[i].Begin-Clusters.back().back().Begin>MaxDis)
                {
                    if (Clusters.back().size()<MinSupport) Clusters.pop_back();
                    Clusters.push_back(vector<Signature>());
                }
                Clusters.back().push_back(SortedSignatures[i]);
            }
        }
        if (Clusters.size()>0 && Clusters.back().size()<MinSupport) Clusters.pop_back();
    }
    else
    {
        float MaxDis=0.7;
        int Size=SortedSignatures.size();
        short *Clustered=(short*) calloc(sizeof(short),Size);
        int F0=0;
        while (F0<Size)
        {
            for (int i=F0;i<SortedSignatures.size();++i)
            {
                if (Clustered[i]) continue;
                if (Clusters.size()==0)
                {
                    Clusters.push_back(vector<Signature>());
                }
                bool AllLink=true;
                for (int j=0;j<Clusters.back().size();++j)
                {
                    float PD;
                    if (distance(SortedSignatures[i],Clusters.back()[j],true,&PD, BamStats)>MaxDis)//+bestPrecision(SortedSignatures[i],Clusters.back()[j])-2>MaxDis)
                    {
                        AllLink=false;
                        if (PD+bestPrecision(SortedSignatures[i],Clusters.back()[j])-2>MaxDis)
                        {
                            Clusters.push_back(vector<Signature>());
                            Clusters.back().push_back(SortedSignatures[i]);
                            break;
                        }
                    }
                }
                if (AllLink)
                {
                    Clusters.back().push_back(SortedSignatures[i]);
                    Clustered[i]=1;
                }
            }
            F0=first0(Clustered,Size,F0);
        }
        free(Clustered);
    }
}

struct ClusterStats
{
    double ABegin,AEnd,ALength;
    set<string> Templates;
    ClusterStats():ABegin(0),AEnd(0),ALength(0),Templates(){}
};

bool isBrother(const Signature &A, const Signature &B, const Arguments *Args, float PosRatio=0.1, int ForceBrother=5, float LengthRatio=0.1, int LengthMinEndurance=5, int ForceBrother2=5, float LengthRatio2=0.1, int NearRange=500, const ClusterStats * AS=nullptr, const ClusterStats* BS=nullptr)
{
    if (A.TemplateName==B.TemplateName) return false;
    double CompareRatio=1.0;
    //DEL: 1.2,0.3
    if (A.HP!=0 && A.HP==B.HP) CompareRatio=Args->HPSameCompareRatio[A.SupportedSV];
    if (A.HP!=0 && B.HP!=0 && A.HP!=B.HP) CompareRatio=Args->HPDiffCompareRatio[A.SupportedSV];
    if (abs(A.Begin-B.Begin)<=ForceBrother*CompareRatio && abs(A.End-B.End)<=ForceBrother*CompareRatio) return true;
    int MinLength=min(A.Length,B.Length);
    int LengthEndurance=max(LengthMinEndurance,int(MinLength*LengthRatio))*CompareRatio;
    if (abs(A.Begin-B.Begin)<=MinLength*PosRatio && abs(A.End-B.End)<=MinLength*PosRatio && abs(A.Length-B.Length)<=LengthEndurance) return true;
    return false;
}

bool isBrother(const Signature* A, const Signature* B, const Arguments *Args, float PosRatio=0.1, int ForceBrother=5, float LengthRatio=0.1, int LengthMinEndurance=5, int ForceBrother2=5, float LengthRatio2=0.1, int NearRange=500, const ClusterStats * AS=nullptr, const ClusterStats* BS=nullptr)
{
    return isBrother(*A, *B, Args, PosRatio, ForceBrother, LengthRatio, LengthMinEndurance, ForceBrother2, LengthRatio2, NearRange, AS, BS);
}

bool isBrother(const vector<Signature> &A, const vector<Signature> &B, const Arguments *Args, float PosRatio=0.1, int ForceBrother=5, float LengthRatio=0.1, int LengthMinEndurance=5, int ForceBrother2=5, float LengthRatio2=0.1, int NearRange=500, const ClusterStats * AS=nullptr, const ClusterStats* BS=nullptr)
{
    if (A[0].SupportedSV!=B[0].SupportedSV) return false;
    NearRange=200;
    if (A.size()==1 && B.size()==1)
    {
        if (abs(A[0].Length-B[0].Length)<ForceBrother2 && abs(A[0].Begin-B[0].Begin)<NearRange && abs(A[0].End-B[0].End)<NearRange) return true;
        return false;
    }
    const set<string> & ATemplates=AS->Templates;
    double MeanLength=(AS->ALength+BS->ALength)/2.0;
    double LengthDiff=ForceBrother2;
    double LengthRatioDiff=LengthRatio2;
    double MinSize=MIN(A.size(),B.size());
    LengthDiff/=MinSize;
    LengthRatioDiff/=MinSize;
    if (abs(AS->ABegin-BS->ABegin)<NearRange && abs(AS->AEnd-BS->AEnd)<NearRange && (((abs(AS->ALength-MeanLength)/MeanLength)<LengthRatioDiff)||(abs(AS->ALength-BS->ALength)<LengthDiff))) 
    {
        for (int i=0;i<B.size();++i)
        {
            if (ATemplates.count(B[i].TemplateName)!=0) return false;//has sigs from same template, no brother
        }
        return true;
    }
    return false;
}

inline int getSupportedSV(const vector<Signature> & Cluster)
{
    return Cluster[0].SupportedSV;
}

inline int getSupportedSV(const vector<Signature*> & Cluster)
{
    return Cluster[0]->SupportedSV;
}

inline int getSupportedSV(const vector<vector<Signature>> & Cluster)
{
    return Cluster[0][0].SupportedSV;
}

void setValues(Signature & S, unsigned& MinBegin, unsigned& MaxBegin, unsigned& MinEnd, unsigned & MaxEnd, unsigned& SignatureCount, vector<ClusterStats>& CStats)
{
    MinBegin=S.Begin;
    MaxBegin=S.Begin;
    MinEnd=S.End;
    MaxEnd=S.End;
    SignatureCount=1;
    CStats.push_back(ClusterStats());
}

void setValues(Signature* & S, unsigned& MinBegin, unsigned& MaxBegin, unsigned& MinEnd, unsigned & MaxEnd, unsigned& SignatureCount, vector<ClusterStats>& CStats)
{
    setValues(*S, MinBegin, MaxBegin, MinEnd, MaxEnd, SignatureCount, CStats);
}

void setValues(vector<Signature> & Ss, unsigned& MinBegin, unsigned& MaxBegin, unsigned& MinEnd, unsigned & MaxEnd, unsigned& SignatureCount, vector<ClusterStats>& CStats)
{
    ClusterStats CS;
    MinBegin=Ss[0].Begin;
    MaxBegin=Ss[0].Begin;
    MinEnd=Ss[0].End;
    MaxEnd=Ss[0].End;
    for (int i=0;i<Ss.size();++i)
    {
        Signature & S=Ss[i];
        MinBegin=MIN(MinBegin,S.Begin);
        MaxBegin=MAX(MaxBegin,S.Begin);
        MinEnd=MIN(MinEnd,S.End);
        MaxEnd=MAX(MaxEnd,S.End);
        CS.ABegin=(CS.ABegin*i)/(i+1)+Ss[i].Begin/(double(i+1));
        CS.AEnd=(CS.AEnd*i)/(i+1)+Ss[i].End/(double(i+1));
        CS.ALength=(CS.ALength*i)/(i+1)+Ss[i].Length/(double(i+1));
        CS.Templates.insert(Ss[i].TemplateName);
    }
    CStats.push_back(CS);
    SignatureCount=Ss.size();
}

template<typename T> class Brotherhood
{
    public:
    vector<T> Cluster;
    bool CCS;
    unsigned int MinBegin, MaxBegin, MinEnd, MaxEnd, SignatureCount;
    vector<ClusterStats> CStats;
    Brotherhood(bool IsCCS=false):Cluster(),MinBegin(0),MaxEnd(0),CCS(IsCCS) {};
    Brotherhood(T &S, bool IsCCS=false):Cluster(),MinBegin(0),MaxEnd(0),CCS(IsCCS) {setCluster(S);};
    void setCluster(T & S)
    {
        Cluster.clear();
        Cluster.push_back(S);
        setValues(S,MinBegin, MaxBegin, MinEnd, MaxEnd, SignatureCount,CStats);
    }
    bool canMerge(const Brotherhood<T> & Other, Arguments &Args) const
    {
        int SupportedSV=getSupportedSV(Cluster);
        int ForceBrother=Args.BrotherhoodTypeForceBrothers[SupportedSV];
        float Ratio=Args.BrotherhoodTypeRatios[SupportedSV];
        int LengthMinEndurance=Args.BrotherhoodTypeLengthMinEndurance[SupportedSV];
        float LengthRatio=Args.BrotherhoodTypeLengthRatios[SupportedSV];
        int NearRange=Args.BrotherhoodNearRanges[SupportedSV];
        int ForceBrother2=Args.BrotherhoodTypeForceBrothers2[SupportedSV];
        float LengthRatio2=Args.BrotherhoodTypeLengthRatios2[SupportedSV];
        if (CCS)
        {
            ForceBrother=Args.BrotherhoodCCSTypeForceBrothers[SupportedSV];
            Ratio=Args.BrotherhoodCCSTypeRatios[SupportedSV];
        }
        if (((Other.MinBegin>MaxEnd+ForceBrother+NearRange) || ((MinBegin>Other.MaxEnd+ForceBrother+NearRange)))) return false;
        for (int i=0;i<Cluster.size();++i)
        {
            const T & A=Cluster[i];
            for (int j=0;j<Other.Cluster.size();++j)
            {
                const T & B=Other.Cluster[j];
                if (isBrother(A,B,&Args,Ratio,ForceBrother,LengthRatio,LengthMinEndurance,ForceBrother2,LengthRatio2,NearRange,&(CStats[i]),&(Other.CStats[j])))
                {
                    return true;
                }
            }
        }
        return false;
    }
    bool merge(Brotherhood<T> & Other, Arguments &Args)
    {
        if (canMerge(Other,Args))
        {
            for (T & B:Other.Cluster)
            {
                Cluster.push_back(B);
            }
            MinBegin=min(MinBegin,Other.MinBegin);
            MaxBegin=MAX(MaxBegin,Other.MaxBegin);
            MinEnd=MIN(MinEnd,Other.MinEnd);
            MaxEnd=max(MaxEnd,Other.MaxEnd);
            SignatureCount+=Other.SignatureCount;
            CStats.insert(CStats.end(),make_move_iterator(Other.CStats.begin()),make_move_iterator(Other.CStats.end()));
            return true;
        }
        return false;
    }
};

template<typename T> void brotherClusteringList(list<Brotherhood<T>> &Brotherhoods, Arguments &Args)
{
    int Round=0;
    while (1)
    {
        ++Round;
        Brotherhoods.sort([](Brotherhood<T> & a, Brotherhood<T> &b)-> bool {return a.MinBegin<b.MinBegin;});
        bool Next=false;
        for (auto Ai=Brotherhoods.begin();Ai!=Brotherhoods.end();++Ai)
        {
            for (auto Bi=next(Ai);Bi!=Brotherhoods.end();++Bi)
            {
                if (Ai->MaxEnd+MIN(Args.ClusteringMaxMergeRange,0.5*(Ai->MaxEnd-Ai->MinBegin))<Bi->MinBegin) break;
                if (Ai->merge(*Bi,Args))
                {
                    Brotherhoods.erase(Bi);
                    Next=true;
                    break;
                }
            }
        }
        if (!Next) break;
    }
    Brotherhoods.sort([](Brotherhood<T> & a, Brotherhood<T> &b)-> bool {return a.MinBegin<b.MinBegin;});
}

template<typename T> void brotherClustering(vector<T> & SortedTs, list<Brotherhood<T>>& Brotherhoods, Stats BamStats, Arguments &Args)
{
    //Need improve 1 thread efficiency for large size
    if (/*Args.ThreadN==1 || */SortedTs.size()<=Args.ClusteringBatchSize*1.1)
    {
        for (T & S:SortedTs) Brotherhoods.push_back(Brotherhood<T>(S,Args.AllCCS));
        brotherClusteringList(Brotherhoods,Args);
    }
    else
    {
        vector<list<Brotherhood<T>>> BatchBrotherhoods;
        for (int i=0;i<SortedTs.size();++i)
        {
            if (i%Args.ClusteringBatchSize==0) BatchBrotherhoods.push_back(list<Brotherhood<T>>());
            BatchBrotherhoods[BatchBrotherhoods.size()-1].push_back(Brotherhood<T>(SortedTs[i],Args.AllCCS));
        }
        for (int i=0;i<BatchBrotherhoods.size();++i)
        {
            brotherClusteringList(BatchBrotherhoods[i],Args);
        }
        vector<list<Brotherhood<T>>> Edges;
        for (int i=0;i<BatchBrotherhoods.size()-1;++i)
        {
            int Left,Right;
            Left=BatchBrotherhoods[i+1].front().MinBegin-Args.ClusteringMaxMergeRange;
            Right=BatchBrotherhoods[i].back().MaxBegin+Args.ClusteringMaxMergeRange;
            auto LeftMost=BatchBrotherhoods[i].begin();
            if (i!=0) advance(LeftMost,unsigned(BatchBrotherhoods[i].size()/2.0));
            auto RightMost=BatchBrotherhoods[i+1].end();
            if (i+1!=BatchBrotherhoods.size()-1)
            {
                RightMost=BatchBrotherhoods[i+1].begin();
                advance(RightMost,unsigned(BatchBrotherhoods[i+1].size()/2.0));
            }
            Edges.push_back(list<Brotherhood<T>>());
            typename list<Brotherhood<T>>::iterator LeftIter=BatchBrotherhoods[i].end();
            --LeftIter;
            typename list<Brotherhood<T>>::iterator RightIter=BatchBrotherhoods[i+1].begin();
            for (;LeftIter!=LeftMost;--LeftIter)
            {
                if (LeftIter->MaxBegin<Left) break;
                if (LeftIter==LeftMost)
                {
                    if (i!=0 && !(LeftIter->MaxBegin<Left))
                    Args.Log.verbose("Risk of not fully clustered detected.");//Change batch size to sig count related to avoid? Or, can be avoided by extend the boundaries. The influence is trivial, so it's not urgent to solve it.
                    break;
                }
            }
            for (;RightIter!=RightMost;++RightIter)
            {
                if (RightIter->MinBegin>Right) break;
            }
            if (RightIter==RightMost && i+1!=BatchBrotherhoods.size()-1)
            Args.Log.verbose("Risk of not fully clustered detected.");
            Edges[Edges.size()-1].splice(Edges[Edges.size()-1].end(),BatchBrotherhoods[i],LeftIter,BatchBrotherhoods[i].end());
            Edges[Edges.size()-1].splice(Edges[Edges.size()-1].end(),BatchBrotherhoods[i+1],BatchBrotherhoods[i+1].begin(),RightIter);
        }
        for (int i=0;i<Edges.size();++i)
        {
            brotherClusteringList(Edges[i],Args);
        }
        for (int i=0;i<BatchBrotherhoods.size();++i)
        {
            Brotherhoods.insert(Brotherhoods.end(),make_move_iterator(BatchBrotherhoods[i].begin()),make_move_iterator(BatchBrotherhoods[i].end()));
        }
        for (int i=0;i<Edges.size();++i)
        {
            Brotherhoods.insert(Brotherhoods.end(),make_move_iterator(Edges[i].begin()),make_move_iterator(Edges[i].end()));
        }
        Brotherhoods.sort([](Brotherhood<T> & a, Brotherhood<T> &b)-> bool {return a.MinBegin<b.MinBegin;});
    }
}

//To reduce memory usage
template<typename T> void brotherClustering(vector<T> & SortedTs, list<Brotherhood<T*>>& Brotherhoods, Stats BamStats, Arguments &Args)
{
    T* PT;
    //Need improve 1 thread efficiency for large size
    if (/*Args.ThreadN==1 || */SortedTs.size()<=Args.ClusteringBatchSize*1.1)
    {
        for (T & S:SortedTs)
        {
            PT=&S;
            Brotherhoods.push_back(Brotherhood<T*>(PT,Args.AllCCS));
        }
        brotherClusteringList(Brotherhoods,Args);
    }
    else
    {
        vector<list<Brotherhood<T*>>> BatchBrotherhoods;
        for (int i=0;i<SortedTs.size();++i)
        {
            if (i%Args.ClusteringBatchSize==0) BatchBrotherhoods.push_back(list<Brotherhood<T*>>());
            PT=&SortedTs[i];
            BatchBrotherhoods[BatchBrotherhoods.size()-1].push_back(Brotherhood<T*>(PT,Args.AllCCS));
        }
        for (int i=0;i<BatchBrotherhoods.size();++i)
        {
            brotherClusteringList(BatchBrotherhoods[i],Args);
        }
        vector<list<Brotherhood<T*>>> Edges;
        for (int i=0;i<BatchBrotherhoods.size()-1;++i)
        {
            int Left,Right;
            Left=BatchBrotherhoods[i+1].front().MinBegin-Args.ClusteringMaxMergeRange;
            Right=BatchBrotherhoods[i].back().MaxBegin+Args.ClusteringMaxMergeRange;
            auto LeftMost=BatchBrotherhoods[i].begin();
            if (i!=0) advance(LeftMost,unsigned(BatchBrotherhoods[i].size()/2.0));
            auto RightMost=BatchBrotherhoods[i+1].end();
            if (i+1!=BatchBrotherhoods.size()-1)
            {
                RightMost=BatchBrotherhoods[i+1].begin();
                advance(RightMost,unsigned(BatchBrotherhoods[i+1].size()/2.0));
            }
            Edges.push_back(list<Brotherhood<T*>>());
            typename list<Brotherhood<T*>>::iterator LeftIter=BatchBrotherhoods[i].end();
            --LeftIter;
            typename list<Brotherhood<T*>>::iterator RightIter=BatchBrotherhoods[i+1].begin();
            for (;LeftIter!=LeftMost;--LeftIter)
            {
                if (LeftIter->MaxBegin<Left) break;
                if (LeftIter==LeftMost)
                {
                    if (i!=0 && !(LeftIter->MaxBegin<Left))
                    Args.Log.verbose("Risk of not fully clustered detected.");//Change batch size to sig count related to avoid? Or, can be avoided by extend the boundaries. The influence is trivial, so it's not urgent to solve it.
                    break;
                }
            }
            for (;RightIter!=RightMost;++RightIter)
            {
                if (RightIter->MinBegin>Right) break;
            }
            if (RightIter==RightMost && i+1!=BatchBrotherhoods.size()-1)
            Args.Log.verbose("Risk of not fully clustered detected.");
            Edges[Edges.size()-1].splice(Edges[Edges.size()-1].end(),BatchBrotherhoods[i],LeftIter,BatchBrotherhoods[i].end());
            Edges[Edges.size()-1].splice(Edges[Edges.size()-1].end(),BatchBrotherhoods[i+1],BatchBrotherhoods[i+1].begin(),RightIter);
        }
        for (int i=0;i<Edges.size();++i)
        {
            brotherClusteringList(Edges[i],Args);
        }
        for (int i=0;i<BatchBrotherhoods.size();++i)
        {
            Brotherhoods.insert(Brotherhoods.end(),make_move_iterator(BatchBrotherhoods[i].begin()),make_move_iterator(BatchBrotherhoods[i].end()));
        }
        for (int i=0;i<Edges.size();++i)
        {
            Brotherhoods.insert(Brotherhoods.end(),make_move_iterator(Edges[i].begin()),make_move_iterator(Edges[i].end()));
        }
        Brotherhoods.sort([](Brotherhood<T*> & a, Brotherhood<T*> &b)-> bool {return a.MinBegin<b.MinBegin;});
    }
}

void brotherClustering(vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, Stats BamStats, Arguments &Args)
{
    list<Brotherhood<Signature*>> Brotherhoods;
    brotherClustering(SortedSignatures,Brotherhoods,BamStats,Args);
    for (Brotherhood<Signature*>& B :Brotherhoods)
    {
        Clusters.push_back(vector<Signature>());
        for (int i=0;i<B.Cluster.size();++i) Clusters.back().push_back(move(*B.Cluster[i]));
    }
}

//May duplicate some signatures.
void fitInClusters(vector<vector<Signature>> &Clusters, vector<Signature> & OtherSignatures, vector<Signature> &Leftovers, Arguments& Args)
{
    vector<bool> Merged;
    for (int i=0;i<OtherSignatures.size();++i) Merged.push_back(false);
    for (vector<Signature>& Cluster : Clusters)
    {
        int ForceBrother=Args.BrotherhoodTypeForceBrothers[Cluster[0].SupportedSV];
        float Ratio=Args.BrotherhoodTypeRatios[Cluster[0].SupportedSV];
        int LengthMinEndurance=Args.BrotherhoodTypeLengthMinEndurance[Cluster[0].SupportedSV];
        float LengthRatio=Args.BrotherhoodTypeLengthRatios[Cluster[0].SupportedSV];
        double ClusterSize=Cluster.size();
        double SupportedBrother=0;
        for (int i=0;i<OtherSignatures.size();++i)
        {
            SupportedBrother=0;
            for (const Signature & s: Cluster)
            {
                if (isBrother(s, OtherSignatures[i],&Args,Ratio,ForceBrother,LengthRatio,LengthMinEndurance))
                {
                    SupportedBrother+=1;
                    if (SupportedBrother/ClusterSize>=0.8)
                    {
                        Cluster.push_back(OtherSignatures[i]);
                        Merged[i]=true;
                        break;
                    }
                }
            }
        }
    }
    for (int i=0;i<OtherSignatures.size();++i)
    {
        if (!Merged[i]) Leftovers.push_back(OtherSignatures[i]);
    }
}

void batchClustering(vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, Stats BamStats, Arguments& Args)
{
    vector<Signature> BaseSignatures,OtherSignatures;
    vector<double> Qualities;
    for (int i=0;i<SortedSignatures.size();++i)
    {
        Qualities.push_back(SortedSignatures[i].Quality);
    }
    sort(Qualities.data(),Qualities.data()+Qualities.size());
    double Quantile=0;
    if (Qualities.size()!=0) Quantile=Qualities[int(Qualities.size()*0.2)];
    for (int i=0;i<SortedSignatures.size();++i)
    {
        if (SortedSignatures[i].Quality>=Quantile) BaseSignatures.push_back(SortedSignatures[i]);
        else OtherSignatures.push_back(SortedSignatures[i]);
    }
    brotherClustering(BaseSignatures,Clusters,BamStats,Args);
    vector<Signature> Leftovers;
    fitInClusters(Clusters, OtherSignatures, Leftovers, Args);
    fprintf(stderr, "Quantile: %lf, Base: %ld, Other: %ld, Leftover: %ld\n",Quantile, BaseSignatures.size(),OtherSignatures.size(),Leftovers.size());
}


void mergeClusters(vector<vector<Signature>> &Clusters, Stats BamStats, Arguments& Args)
{
    for (int i=0;i<Clusters.size();++i)
    {
        if (Clusters[i].size()==0) continue;
        for (int j=0;j<Clusters.size();++j)
        {
            if (j==i) continue;
            if (Clusters[j].size()==0) continue;
            if(isBrother(Clusters[i],Clusters[j],&Args))
            {
                Clusters[i].insert(Clusters[i].end(),make_move_iterator(Clusters[j].begin()),make_move_iterator(Clusters[j].end()));
                if(Clusters[j].size()!=0) Clusters[j].clear();
            }
        }
    }
}

inline void addSigMap(unordered_map<string,list<list<Brotherhood<Signature>>::iterator>> &PotentialMerges, list<Brotherhood<Signature>>::iterator & Iter)
{
    if (PotentialMerges.count(Iter->Cluster[0].TemplateName)==0) PotentialMerges[Iter->Cluster[0].TemplateName]=list<list<Brotherhood<Signature>>::iterator>();
    list<list<Brotherhood<Signature>>::iterator>& TheList=PotentialMerges[Iter->Cluster[0].TemplateName];
    TheList.push_back(Iter);
}

inline void addSigMap(unordered_map<string,list<list<Signature>::iterator>> &PotentialMerges, list<Signature>::iterator & Iter)
{
    if (PotentialMerges.count(Iter->TemplateName)==0) PotentialMerges[Iter->TemplateName]=list<list<Signature>::iterator>();
    list<list<Signature>::iterator>& TheList=PotentialMerges[Iter->TemplateName];
    TheList.push_back(Iter);
}

inline unsigned getBegin(list<list<Brotherhood<Signature>>::iterator>::iterator & Iter)
{
    return (*Iter)->MinBegin;
}

inline unsigned getEnd(list<list<Brotherhood<Signature>>::iterator>::iterator & Iter)
{
    return (*Iter)->MaxEnd;
}

inline unsigned getBegin(list<Brotherhood<Signature>>::iterator & Iter)
{
    return Iter->MinBegin;
}

inline unsigned getEnd(list<Brotherhood<Signature>>::iterator & Iter)
{
    return Iter->MaxEnd;
}

inline unsigned getBegin(list<list<Signature>::iterator>::iterator & Iter)
{
    return (*Iter)->Begin;
}

inline unsigned getEnd(list<list<Signature>::iterator>::iterator & Iter)
{
    return (*Iter)->End;
}

inline unsigned getBegin(list<Signature>::iterator & Iter)
{
    return Iter->Begin;
}

inline unsigned getEnd(list<Signature>::iterator & Iter)
{
    return Iter->End;
}

inline unsigned getSupportedSV(const list<list<Brotherhood<Signature>>::iterator>::iterator & Iter)
{
    return (*Iter)->Cluster[0].SupportedSV;
}

inline unsigned getSupportedSV(const list<list<Signature>::iterator>::iterator & Iter)
{
    return (*Iter)->SupportedSV;
}

inline unsigned getSupportedSV(const list<Signature>::iterator & Iter)
{
    return Iter->SupportedSV;
}

//All Brotherhoods in Sigs have only one Signature.
//Can be DP optimized.
template<typename TypeIterType> bool BestFit(list<TypeIterType> &Sigs, double Begin, double End, typename list<TypeIterType>::iterator &FitBegin, typename list<TypeIterType>::iterator &FitEnd, unsigned &FitLength)
{
    int SupportedSV=getSupportedSV(Sigs.begin());
    int MaxBeginDiff=1000000;
    double LengthEndureMin=10;
    double LengthEndureRatio=0.1;
    if (SupportedSV==1)
    {
        MaxBeginDiff=200;
        LengthEndureRatio=1.0;
    }
    double Length=End-Begin;
    double LengthMin=MIN(Length-LengthEndureMin,Length*(1.0-LengthEndureRatio));
    double LengthMax=MAX(Length+LengthEndureMin,Length*(1.0+LengthEndureRatio));
    Sigs.sort([](TypeIterType & a, TypeIterType &b)-> bool {return getBegin(a)<getBegin(b);});
    unsigned MinDiff=0xffffffff;
    unsigned MinDiffMergeCount=0;
    unsigned MergeCount=0;
    for (auto Iter=Sigs.begin();Iter!=Sigs.end();++Iter)
    {
        double LengthAdd=getEnd(Iter)-getBegin(Iter);
        double ThisBegin=getBegin(Iter);
        double ThisEnd=getEnd(Iter);
        MergeCount=1;
        for (auto Jter=Iter;Jter!=Sigs.end();++Jter)
        {
            if (Jter==Iter) continue;
            if (abs(getBegin(Jter)-ThisBegin)<MaxBeginDiff)
            {
                if (getBegin(Jter)>ThisBegin) break;
                continue;
            }
            LengthAdd+=getEnd(Jter)-getBegin(Jter);
            ThisEnd=ThisBegin+LengthAdd;
            ++MergeCount;
            if (LengthAdd>=LengthMin && LengthAdd<=LengthMax)
            {
                unsigned Diff=ThisBegin-Begin+ThisEnd-End;
                if (Diff<MinDiff)
                {
                    MinDiff=Diff;
                    MinDiffMergeCount=MergeCount;
                    FitBegin=Iter;
                    FitEnd=Jter;
                    FitLength=LengthAdd;
                    ++FitEnd;
                }
            }
        }
    }
    if (MinDiff<200) return true;
    return false;
}

void mergeAdjacent(list<Brotherhood<Signature>> &SignatureBrotherhoods, Arguments & Args)//deprecated
{
    int Range=500;
    int MaxMergeLength=2000000;
    list<Brotherhood<Signature>> MergeCandidates;
    for (auto Iter=SignatureBrotherhoods.begin();Iter!=SignatureBrotherhoods.end();++Iter)
    {
        if ((Iter->MinEnd-Iter->MinBegin)<=MaxMergeLength)
        {
            list<Brotherhood<Signature>>::iterator PreviousIter=Iter;
            --PreviousIter;
            MergeCandidates.splice(MergeCandidates.end(),SignatureBrotherhoods,Iter);
            Iter=PreviousIter;
        }
    }
    fprintf(stderr,"Brotherhoods:%lu , merge candidates:%lu\n",SignatureBrotherhoods.size()+MergeCandidates.size(),MergeCandidates.size());
    list<Brotherhood<Signature>> NewBros;
    unsigned Merged=0;
    for (auto Iter=MergeCandidates.begin();Iter!=MergeCandidates.end();++Iter)
    {
        if ((Iter->MinEnd-Iter->MinBegin)<50) continue;
        unsigned SearchLeft=Iter->MinBegin-Range;
        unsigned SearchRight=Iter->MaxEnd+Range;
        unordered_map<string,list<list<Brotherhood<Signature>>::iterator>> PotentialMerges;
        for (auto Jter=Iter;Jter!=MergeCandidates.end();++Jter)
        {
            if (Jter==Iter) continue;
            if (Jter->MinBegin>SearchRight) break;
            if (Jter->SignatureCount==1)
            {
                // list<Brotherhood<Signature>>::iterator PreviousIter=Jter;
                // --PreviousIter;
                // addSigMap(PotentialMerges,MergeCandidates,Jter,PreviousIter);
                addSigMap(PotentialMerges,Jter);
            }
        }
        for (auto Jter=Iter;true;--Jter)
        {
            if (Jter==Iter)
            {
                if (Jter==MergeCandidates.begin()) break;
                continue;
            }
            if (Jter->MinBegin<SearchLeft) break;
            if (Jter->SignatureCount==1)
            {
                addSigMap(PotentialMerges,Jter);
            }
            if (Jter==MergeCandidates.begin()) break;
        }
        list<list<Brotherhood<Signature>>::iterator>::iterator FitBegin, FitEnd;
        unsigned FitLength;
        for (auto mi=PotentialMerges.begin();mi!=PotentialMerges.end();++mi)
        {
            if (BestFit(mi->second,Iter->MinBegin,Iter->MinEnd,FitBegin,FitEnd,FitLength))
            {
                Signature NewSig((*FitBegin)->Cluster[0]);
                NewSig.End=NewSig.Begin+FitLength;
                NewSig.Length=FitLength;
                if (NewSig.SupportedSV==1)
                {
                    NewSig.InsBases="";
                }
                for (auto I=FitBegin;I!=FitEnd;++I)
                {
                    NewSig.InsBases+=(*I)->Cluster[0].InsBases;
                    MergeCandidates.erase((*I));
                }
                Brotherhood<Signature> NewBro(NewSig,(*FitBegin)->CCS);
                if (!Iter->merge(NewBro,Args))
                {
                    NewBros.push_back(NewBro);
                }
                else
                {
                    ++Merged;
                }
            }
        }
    }
    fprintf(stderr,"NewBros:%lu, Merged:%lu\n",NewBros.size(),Merged);
    MergeCandidates.splice(MergeCandidates.end(),NewBros);
    SignatureBrotherhoods.splice(SignatureBrotherhoods.end(),MergeCandidates);
    SignatureBrotherhoods.sort([](Brotherhood<Signature> & a, Brotherhood<Signature> &b)-> bool {return a.MinBegin<b.MinBegin;});
}

//Merge adjacent sigs that not clustered in the first clustering and add to brotherhoods for the next step
//Will change ClusteredBrotherhoods. Will not change SortedSignatures.
void genMergedSigs(list<Brotherhood<Signature>> &ClusteredBrotherhoods, vector<Signature> & SortedSignatures, int *ClusterCount, Stats BamStats, Arguments & Args)
{
    for (Brotherhood<Signature>& C : ClusteredBrotherhoods)
    {
        for (int i=0;i<C.Cluster.size();++i) ClusterCount[C.Cluster[i].ID]=C.Cluster.size();
    }
    unordered_map<string,vector<vector<Signature>::iterator>> TemplateSigs;
    vector<Signature> Sigs=SortedSignatures;
    int SupportedSV=0;
    if (Sigs.size()>0) SupportedSV=Sigs[0].SupportedSV;
    for (int i=0;i<SortedSignatures.size();++i)
    {
        if (ClusterCount[SortedSignatures[i].ID]>1) continue;
        if (TemplateSigs.count(SortedSignatures[i].TemplateName)==0) TemplateSigs[SortedSignatures[i].TemplateName]=vector<vector<Signature>::iterator>();
        TemplateSigs[SortedSignatures[i].TemplateName].push_back(SortedSignatures.begin()+i);
    }
    if (SupportedSV!=1)
    {
        int Endurance=200;
        for (auto mi=TemplateSigs.begin();mi!=TemplateSigs.end();++mi)
        {
            sort(mi->second.begin(),mi->second.end(),[](vector<Signature>::iterator & a, vector<Signature>::iterator &b)-> bool {return a->Begin<b->Begin;});
            int NewBegin=mi->second[0]->Begin;
            int Merged=0;
            string NewInsBases="";
            for (int i=1;i<mi->second.size();++i)
            {
                Signature &Former=*mi->second[i-1];
                Signature &Latter=*mi->second[i];
                if (Latter.Begin>=Former.End && Latter.Begin<Former.End+Endurance)
                {
                    ++Merged;
                    NewInsBases+=Latter.InsBases;
                }
                else
                {
                    Signature NewSig(Latter.Type,Latter.Tech,Latter.SupportedSV,NewBegin,Latter.End,Latter.TemplateName,Latter.Quality,NewInsBases.c_str());
                    NewSig.Artificial=true;
                    if (Merged>0) Sigs.push_back(NewSig);
                    Merged=0;
                    NewInsBases="";
                    NewBegin=Latter.Begin;
                }
            }
        }
    }
    else
    {
        int Endurance=500;
        for (auto mi=TemplateSigs.begin();mi!=TemplateSigs.end();++mi)
        {
            sort(mi->second.begin(),mi->second.end(),[](vector<Signature>::iterator & a, vector<Signature>::iterator &b)-> bool {return a->Begin<b->Begin;});
            int NewBegin=mi->second[0]->Begin;
            int Merged=0;
            string NewInsBases="";
            for (int i=1;i<mi->second.size();++i)
            {
                Signature &Former=*mi->second[i-1];
                Signature &Latter=*mi->second[i];
                if (Latter.Begin>=NewBegin && Latter.Begin<NewBegin+Endurance)
                {
                    ++Merged;
                    NewInsBases+=Latter.InsBases;
                }
                else
                {
                    Signature NewSig(Latter.Type,Latter.Tech,Latter.SupportedSV,NewBegin,NewBegin+NewInsBases.length(),Latter.TemplateName,Latter.Quality,NewInsBases.c_str());
                    NewSig.Artificial=true;
                    if (Merged>0) Sigs.push_back(NewSig);
                    Merged=0;
                    NewInsBases="";
                    NewBegin=Latter.Begin;
                }
            }
        }

    }
    sort(Sigs.begin(),Sigs.end(),[](Signature& a, Signature&b)-> bool {return a.Begin<b.Begin;});
    ClusteredBrotherhoods.clear();
    brotherClustering(Sigs,ClusteredBrotherhoods,BamStats,Args);
}

bool allArti(Brotherhood<Signature> &B)
{
    for (Signature & s : B.Cluster)
    {
        if (!s.Artificial) return false;
    }
    return true;
}

void mergeAdjacent(list<Brotherhood<Signature>> &ClusteredBrotherhoods, vector<Signature> & SortedSignatures, Stats BamStats, Arguments & Args)
{
    int Range=500;
    int MaxMergeLength=2000;
    int MinKeep=1;
    int CanMergeMaxClusterCount=5;
    
    int ClusterCount[SortedSignatures.size()];
    genMergedSigs(ClusteredBrotherhoods,SortedSignatures,ClusterCount,BamStats,Args);

    list<Brotherhood<Signature>> MergeCandidates;
    for (auto Iter=ClusteredBrotherhoods.begin();Iter!=ClusteredBrotherhoods.end();++Iter)
    {
        if ((Iter->MinEnd-Iter->MinBegin)<=MaxMergeLength)
        {
            list<Brotherhood<Signature>>::iterator PreviousIter=Iter;
            --PreviousIter;
            MergeCandidates.splice(MergeCandidates.end(),ClusteredBrotherhoods,Iter);
            Iter=PreviousIter;
        }
    }
    
    list<Signature> Sigs, NewSigs;
    for (int i=0;i<SortedSignatures.size();++i) Sigs.push_back(SortedSignatures[i]);
    
    list<Signature>::iterator SearchBegin=Sigs.begin();
    for (auto Iter=MergeCandidates.begin();Iter!=MergeCandidates.end();++Iter)
    {
        if ((Iter->MinEnd-Iter->MinBegin)<50) continue;
        unsigned SearchLeft=Iter->MinBegin-Range;
        unsigned SearchRight=Iter->MaxEnd+Range;
        unordered_map<string,list<list<Signature>::iterator>> PotentialMerges;
        for (auto Jter=SearchBegin;Jter!=Sigs.end();++Jter)
        {
            if (Jter->Begin>SearchRight) break;
            if (Jter->Begin<SearchLeft)
            {
                SearchBegin=Jter;
                continue;
            }
            if (ClusterCount[Jter->ID]>CanMergeMaxClusterCount) continue;
            addSigMap(PotentialMerges,Jter);
        }
        list<list<Signature>::iterator>::iterator FitBegin, FitEnd;
        unsigned FitLength;
        vector<tuple<list<list<Signature>::iterator>::iterator,list<list<Signature>::iterator>::iterator,unsigned>> Fits;
        for (auto mi=PotentialMerges.begin();mi!=PotentialMerges.end();++mi)
        {
            if (BestFit(mi->second,Iter->MinBegin,Iter->MinEnd,FitBegin,FitEnd,FitLength))
            {
                Fits.push_back(make_tuple(FitBegin,FitEnd,FitLength));
            }
        }
        if (Fits.size()<MinKeep+allArti(*Iter)?1:0) continue;
        for (auto Fit : Fits)
        {
            FitLength=get<2>(Fit);
            FitBegin=get<0>(Fit);
            FitEnd=get<1>(Fit);
            Signature NewSig((*(*FitBegin)));
            NewSig.End=NewSig.Begin+FitLength;
            NewSig.Length=FitLength;
            if (NewSig.SupportedSV==1)
            {
                NewSig.InsBases="";
            }
            for (auto I=FitBegin;I!=FitEnd;++I)
            {
                NewSig.InsBases+=(*I)->InsBases;
                if (*I==SearchBegin)
                {
                    ++SearchBegin;
                }
                Sigs.erase((*I));
            }
            NewSigs.push_back(NewSig);
        }
    }
    fprintf(stderr,"NewSigs:%lu\n",NewSigs.size());
    Sigs.splice(Sigs.end(),NewSigs);
    Sigs.sort([](Signature & a, Signature &b)-> bool {return a.Begin<b.Begin;});
    SortedSignatures.clear();
    for (Signature & S : Sigs) SortedSignatures.push_back(S);
    
    ClusteredBrotherhoods.splice(ClusteredBrotherhoods.end(),MergeCandidates);
    ClusteredBrotherhoods.sort([](Brotherhood<Signature> & a, Brotherhood<Signature> &b)-> bool {return a.MinBegin<b.MinBegin;});
}

//no merging adjacent sigs in input phase, then clustering by two layers here
//first layer clusters all similar sigs, and merge nearby clusters
//second layer cluster all 1st layer clusters.
void brotherClustering2(vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, vector<ClusterCore> &Cores, Stats BamStats, Arguments &Args)
{
    fprintf(stderr,"First layer clustering... %lu %lu:\n",SortedSignatures.size(), Clusters.size());
    list<Brotherhood<Signature>> SignatureBrotherhoods;
    brotherClustering(SortedSignatures,SignatureBrotherhoods,BamStats,Args);
    mergeAdjacent(SignatureBrotherhoods,SortedSignatures,BamStats,Args);
    SignatureBrotherhoods.clear();
    brotherClustering(SortedSignatures,SignatureBrotherhoods,BamStats,Args);
    for (Brotherhood<Signature>& B :SignatureBrotherhoods)
    {
        Clusters.push_back(B.Cluster);
    }
    list<Brotherhood<vector<Signature>>> ClusterBrotherhoods;
    fprintf(stderr,"Second layer clustering... %lu %lu:\n",SortedSignatures.size(), Clusters.size());
    brotherClustering(Clusters, ClusterBrotherhoods,BamStats,Args);
    Clusters.clear();
    //core generation
    for (Brotherhood<vector<Signature>>& B :ClusterBrotherhoods)
    {
        int BiggestSize=0, BiggestBegin=0, BiggestEnd=0, FullSize=0;
        for (int i=0;i<B.Cluster.size();++i)
        {
            if (BiggestSize<B.Cluster[i].size())
            {
                BiggestSize=B.Cluster[i].size();
                BiggestBegin=FullSize;
                BiggestEnd=BiggestBegin+BiggestSize;
            }
            FullSize+=B.Cluster[i].size();
            if (i==0) Clusters.push_back(B.Cluster[i]);
            else Clusters[Clusters.size()-1].insert(Clusters[Clusters.size()-1].end(),make_move_iterator(B.Cluster[i].begin()),make_move_iterator(B.Cluster[i].end()));
        }
        if (BiggestSize>=FullSize*0.7) Cores.push_back(ClusterCore(BiggestBegin,BiggestEnd));
        else Cores.push_back(ClusterCore());
    }
}

void clustering(int SVTypeI, string & ContigName, vector<Signature> & SortedSignatures, vector<vector<Signature>> &Clusters, vector<ClusterCore> &Cores, Stats BamStats, Arguments& Args)
{
    if (SortedSignatures.size()>0 && Args.BrotherhoodNearRanges[SVTypeI]==-1)
    {
        brotherClustering(SortedSignatures,Clusters,BamStats,Args);
    }
    else if (SortedSignatures.size()>0) brotherClustering2(SortedSignatures,Clusters,Cores,BamStats,Args);
    Args.Log.debug("%s of %s: number of signatures: %lu, number of clusters:%lu:",SVTypeNames[SVTypeI],ContigName.c_str(),SortedSignatures.size(), Clusters.size());
}