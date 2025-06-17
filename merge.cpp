#include "merge.h"
#include <limits.h>
#include "htslib/htslib/thread_pool.h"
#include <algorithm>
#include <set>
#include <pthread.h>

using namespace std;

//==========Omni Merge Here==========

vector<int> getAllProperEdges(vector<Signature> &Signatures, int MaxDis=500)// a 0 b 1 c 2 d, return the number of edge that can connect, sorted by edge length
{
	vector<int> Edges;
	for (int i=0;i<Signatures.size()-1;++i)
	{
		// Edges.push_back(i);
		if (Signatures[i+1].Begin-Signatures[i].End<MaxDis)
		{
			Edges.push_back(i);
		}
	}
    sort(Edges.begin(),Edges.end(),[Signatures](int &e1, int&e2){
        return Signatures[e1+1].Begin-Signatures[e1].End<Signatures[e2+1].Begin-Signatures[e2].End;
    });
	return Edges;
}

inline void generateSegsFromSuccessiveLink(vector<Signature> &Sigs, vector<Segment> &Segs, int LinkStart, int LinkEnd, int &ScoreA, int LimitStartI=0, int LimitEndI=-1, int A=500)
{
    Segs.clear();
	ScoreA=0;
	int a=0,b=0;
	for (int i=LimitStartI;i<LinkStart;++i)
	{
        Segs.push_back(Segment(Sigs[i].Begin, Sigs[i].End, Sigs[i].InsBases.c_str()));
	}
    int Length=Sigs[LinkStart].Length;
	string InsBases=Sigs[LinkStart].InsBases;
	for (int i=LinkStart+1;i<=LinkEnd;++i)
	{
        Length+=Sigs[i].Length;
		InsBases+=Sigs[i].InsBases;
	}
    Segs.push_back(Segment(Sigs[LinkStart].Begin,Sigs[LinkStart].Begin+Length,InsBases.c_str()));
    if (LimitEndI==-1) LimitEndI=Sigs.size();
	for (int i=LinkEnd+1;i<LimitEndI;++i)
	{
        Segs.push_back(Segment(Sigs[i].Begin, Sigs[i].End, Sigs[i].InsBases.c_str()));
	}
	ScoreA=A*a-b;
}

inline vector<Segment> generateSegsFromLink(vector<Signature> &Sigs, vector<int> &Link, int &ScoreA, int LimitStartI=0, int LimitEndI=-1, int A=500)
{
	vector<Segment> Segs;
	ScoreA=0;
	int a=0,b=0;
	sort(Link.begin(),Link.end());
	set<int> Linked;
	for (int i=0;i<Link.size();++i)
	{
		for (int j=i+1;j<=Link.size();++j)
		{
			if (j==Link.size() || Link[j]-Link[i]!=j-i)
			{
				int Length=Sigs[Link[i]].Length;
				string InsBases=Sigs[Link[i]].InsBases;
                Linked.insert(Link[i]);
				for (int k=i;k<j;++k)
				{
					Linked.insert(Link[k]+1);
					Length+=Sigs[Link[k]+1].Length;
					InsBases+=Sigs[Link[k]+1].InsBases;
					b+=Sigs[Link[k]+1].Begin-Sigs[Link[k]].End;
					b=max(0,b);
				}
				a+=j-i;
				Segs.push_back(Segment(Sigs[Link[i]].Begin,Sigs[Link[i]].Begin+Length,InsBases.c_str()));
				i=j-1;
				break;
			}
		}
	}
    if (LimitEndI==-1) LimitEndI=Sigs.size();
	for (int i=LimitStartI;i<LimitEndI;++i)
	{
		if (Linked.count(i)==0)
		{
			Segs.push_back(Segment(Sigs[i].Begin, Sigs[i].End,Sigs[i].InsBases.c_str()));
		}
	}
	sort(Segs.begin(),Segs.end(),[](const Segment& a, const Segment& b) {
        return a.Begin < b.Begin;
    });
	ScoreA=A*a-b;
	return Segs;
}

int scoreB(vector<Segment>& Segs1, vector<Segment>& Segs2, int B)//Segs1 and Segs2 shall be sorted
{
	int ScoreB=0;
    //already sorted.
	// sort(Segs1.begin(), Segs1.end(), [](const Segment& a, const Segment& b) {
    //     return a.Begin < b.Begin;
    // });
    // sort(Segs2.begin(), Segs2.end(), [](const Segment& a, const Segment& b) {
    //     return a.Begin < b.Begin;
    // });
    // for(int j = 0; j < Segs2.size(); ++j) {
    //     int i = 0;
    //     while(i < Segs1.size() && Segs1[i].Begin <= Segs2[j].Begin + 10*B) {
    //         int dist = abs(Segs1[i].Begin - Segs2[j].Begin) + abs(Segs1[i].End - Segs2[j].End);
	// 		int difflen=abs((Segs1[i].End-Segs1[1].Begin)-(Segs2[j].End-Segs2[j].Begin));
    //         if(dist <= 10*B && difflen<=B) {
	// 			ScoreB+= B - dist;
	// 			// ScoreB+=B-difflen;
    //             ++i;
    //         } else {
    //             ++i;
    //         }
    //     }
    // }
	for (int i=0;i<Segs1.size();++i)
	{
		int j=0;
		// int Best=0;
		while (j<Segs2.size() && Segs2[j].Begin <=Segs1[i].Begin+ B)
		{
            int dist = abs(Segs1[i].Begin - Segs2[j].Begin) + abs(Segs1[i].End - Segs2[j].End);
            if (dist >B) {++j;continue;}
			int difflen=abs((Segs1[i].End-Segs1[i].Begin)-(Segs2[j].End-Segs2[j].Begin));
			// difflen=max(difflen,abs(int(Segs1[i].InsBases.length())-int(Segs2[j].InsBases.length())));
            if(difflen<=B) {
				// Best=max(Best, B - dist);
				ScoreB+=B-difflen;
                ++j;
            } else {
                ++j;
            }
		}
		// ScoreB+=Best;
	}
    return ScoreB;
}

//get the best score for the combine of a given alignment
int getBestScore(vector<Segment> &Segs, vector<Signature> &Sigs, int B=100)
{
	if (Sigs.size()==0) return 0;
	// if (Sigs.size()>1) return 0;
    bool UseEdge=false;
    if (UseEdge)
    {
        int BestScores[Sigs.size()+1][Sigs.size()+1];//DP scores, [A][B], for Score for Sigs[A,B)
        // fprintf(stderr,"DP space: %lfMB.\n" ,(double)((Sigs.size()+1)*(Sigs.size()+1))*4.0/1024.0/1024.0);
        set<int> AllEdges;
        vector<int> Edges;
        Edges=getAllProperEdges(Sigs);
        for (int i=0;i<Edges.size();++i) AllEdges.insert(Edges[i]);
        for (int len=1; len<=Sigs.size(); ++len)
        {
            for (int i=0; i<=Sigs.size()-len; ++i)
            {
                int j = i+len;
                if ((i!=0 && AllEdges.count(i)==0) || (j!=Sigs.size() && AllEdges.count(j-1)==0)) continue;
                int best = 0;
                for (int k=i+1; k<j; k++){
                    // break at k
                    if (AllEdges.count(k)==1)
                    best = max(best, BestScores[i][k]+BestScores[k][j]);
                }
                vector<int> Link;
                int ScoreA;
                for (int k=i;k<j-1;++k) if (AllEdges.count(k)==1) Link.push_back(k);
                vector<Segment> Segs2=generateSegsFromLink(Sigs,Link,ScoreA,i,j);
                BestScores[i][j] = max(best, scoreB(Segs,Segs2,100));
            }
        }
        return BestScores[0][Edges.size()];
    }
	int BestScores[Sigs.size()+1][Sigs.size()+1];//DP scores, [A][B], for Score for Sigs[A,B)
	// fprintf(stderr,"DP space: %lfMB.\n" ,(double)((Sigs.size()+1)*(Sigs.size()+1))*4.0/1024.0/1024.0);
    vector<Segment> Segs2;
	for (int len=1; len<=Sigs.size(); ++len)
	{
		for (int i=0; i<=Sigs.size()-len; ++i)
		{
			int j = i+len;
			int best = 0;
			for (int k=i+1; k<j; k++){
				best = max(best, BestScores[i][k]+BestScores[k][j]);
			}
			int ScoreA;
            generateSegsFromSuccessiveLink(Sigs,Segs2,i,j-1,ScoreA,i,j);
			BestScores[i][j] = max(best, scoreB(Segs,Segs2,B));
    	}
	}
	return BestScores[0][Sigs.size()];
}

//for every candidate link and return a best merged segments
void forAllCandidateLinks(int Type, int Index, vector<AlignmentSigs> * pAlignmentsSigs, const int *SigIndexes, int RelevantBeginIndexI, int RelevantEndIndexI, vector<Segment> & BestSegs, int & BestBestScore, vector<int> &Edges, Arguments& Args, int MaxEdges=10, int CurrentIndex=0, vector<int> CurrentLink=vector<int>())//shall be recursive
{
	int End=min(int(Edges.size()),MaxEdges);
	vector<Signature> &Sigs=(*pAlignmentsSigs)[Index].TypeSignatures[Type];
	if (CurrentIndex==End)
	{
		int ScoreA=0;
		int ScoreB=0;
		//generate the corresponding segs according to the CurrentLink
		vector<Segment> Segs=generateSegsFromLink(Sigs,CurrentLink,ScoreA,0,-1,Args.OmniA[Type]);
        int CountLimit=Args.OmniBCountLimit;
        int Step=(RelevantEndIndexI-RelevantBeginIndexI)/CountLimit;
        Step=max(1,Step);
        int Count=0;
		for (int i=RelevantBeginIndexI;i<RelevantEndIndexI;i+=Step)
		{
			if (Index==SigIndexes[i]) continue;
            if (Count>CountLimit) break;
			ScoreB+=getBestScore(Segs, (*pAlignmentsSigs)[SigIndexes[i]].TypeSignatures[Type],Args.OmniB[Type]);
            ++Count;
		}
		int BestScore=ScoreA;
		BestScore+=ScoreB*Args.OmniBScoreBRatio;//(RelevantEndIndexI-RelevantBeginIndexI-1);
		if (BestScore>BestBestScore)
		{
			BestBestScore=BestScore;
			BestSegs=Segs;
		}
		return;
	}
	forAllCandidateLinks(Type, Index, pAlignmentsSigs, SigIndexes, RelevantBeginIndexI, RelevantEndIndexI, BestSegs, BestBestScore, Edges, Args, MaxEdges, CurrentIndex+1, CurrentLink);
	CurrentLink.push_back(Edges[CurrentIndex]);
	forAllCandidateLinks(Type, Index, pAlignmentsSigs, SigIndexes, RelevantBeginIndexI, RelevantEndIndexI, BestSegs, BestBestScore, Edges, Args, MaxEdges, CurrentIndex+1, CurrentLink);
}

void getSigsFromSegs(int Type, int HP, vector<Segment> & BestSegs, vector<Signature> &Signatures, int Tech, const char * qname, double Quality)
{
	int End;
	for (int i=0;i<BestSegs.size();++i)
	{
		End=BestSegs[i].End;
		if (BestSegs[i].InsBases.size()!=0) End=BestSegs[i].Begin+BestSegs[i].InsBases.size();
		Signature Temp(0,Tech,Type,BestSegs[i].Begin,End,qname,Quality, BestSegs[i].InsBases.c_str());
		Temp.HP=HP;
		Signatures.push_back(Temp);
	}
}

int searchBegin(const int Begin, const int *MaxEnds, const int Size, int s=0, int S=-1)
{
	if (s==S) return s;
	if (S==-1) S=Size;
	int m=(s+S)/2;
	if (MaxEnds[m]<=Begin) return searchBegin(Begin,MaxEnds,Size, m==s?s+1:m, S);
	return searchBegin(Begin, MaxEnds, Size, s,m);
}

int searchEnd(const int End, int Type, vector<AlignmentSigs> &ASs, const int * SigIndexes, int s=0, int S=-1)
{
	if (s==S) return s;
	if (S==-1) S=ASs.size();
	int m=(s+S)/2;
	if (ASs[SigIndexes[m]].TypeBeginMost[Type]>End) return searchEnd(End, Type, ASs, SigIndexes, s, m);
	return searchEnd(End, Type, ASs, SigIndexes,m==s?s+1:m, S);
}

void getRelevants(int Type, AlignmentSigs *pAlignmentSigs, vector<AlignmentSigs> * pAlignmentsSigs, const int * SigIndexes, int &RelevantBeginI, int &RelevantEndI, const int * MaxEnds)
{
	RelevantBeginI=searchBegin(pAlignmentSigs->BeginMost,MaxEnds, pAlignmentsSigs->size());
	RelevantEndI=searchEnd(pAlignmentSigs->BeginMost,Type,*pAlignmentsSigs, SigIndexes);
}

void omniBMergeType(int Type, vector<Signature> &Signatures, int Index, vector<AlignmentSigs> * pAlignmentsSigs, const int *SigIndexes, const int *MaxEnds, Arguments *pArgs)
{
	if (pAlignmentsSigs->at(Index).TypeSignatures[Type].size()==0) return;
	int RelevantBeginIndexI, RelevantEndIndexI;
	getRelevants(Type, &(*pAlignmentsSigs)[Index], pAlignmentsSigs, SigIndexes, RelevantBeginIndexI, RelevantEndIndexI, MaxEnds);
	vector<int> Edges=getAllProperEdges(pAlignmentsSigs->at(Index).TypeSignatures[Type]);
	vector<Segment> BestSegs;
	int BestBestScore=INT_MIN;
	forAllCandidateLinks(Type, Index, pAlignmentsSigs, SigIndexes, RelevantBeginIndexI, RelevantEndIndexI,BestSegs,BestBestScore, Edges, *pArgs, pArgs->OmniBMaxEdges);
	getSigsFromSegs(Type,pAlignmentsSigs->at(Index).HP,BestSegs,Signatures,(*pAlignmentsSigs)[Index].TypeSignatures[Type][0].Tech,(*pAlignmentsSigs)[Index].TypeSignatures[Type][0].TemplateName.c_str(),(*pAlignmentsSigs)[Index].TypeSignatures[Type][0].Quality);
}

void omniBMerge(vector<vector<Signature>> * pTypeSignatures, int Index, vector<AlignmentSigs> * pAlignmentsSigs, const int **TypeSigIndexes, const int **TypeMaxEnds, Arguments *pArgs)
{
    if (Index % 10000==0) pArgs->Log.verbose("Merged %d",Index);
	omniBMergeType(0, (*pTypeSignatures)[0], Index, pAlignmentsSigs, TypeSigIndexes[0], TypeMaxEnds[0], pArgs);
	omniBMergeType(1, (*pTypeSignatures)[1], Index, pAlignmentsSigs, TypeSigIndexes[0], TypeMaxEnds[0], pArgs);
}

//for multithreading

void *omniBHandler(void * Args)
{
	OmniBMergeArgs * A=(OmniBMergeArgs*)Args;
	omniBMerge(&(A->pTypeSignatures->at(pthread_self())),A->Index,A->pAlignmentsSigs,A->TypeSigIndexes,A->TypeMaxEnds,A->pArgs);
	delete A;
	return NULL;
}

AlignmentSigs::AlignmentSigs(unsigned long long ID, const char *TempName):AlignmentID(ID),TemplateName(TempName),TypeSignatures({std::vector<Signature>(),std::vector<Signature>()}),BeginMost(-1),EndMost(-1),TypeBeginMost({-1,-1}),TypeEndMost({-1,-1}),HP(0)
{
}

int AlignmentSigs::getBeginMost()
{
	if (BeginMost==-1)
	{
		if (TypeSignatures[0].size()>0) TypeBeginMost[0]=TypeSignatures[0][0].Begin;
		if (TypeSignatures[1].size()>0) TypeBeginMost[1]=TypeSignatures[1][0].Begin;
		BeginMost=TypeBeginMost[0]==-1?TypeBeginMost[1]:(TypeBeginMost[1]==-1?TypeBeginMost[0]:min(TypeBeginMost[0],TypeBeginMost[1]));
	}
	return BeginMost;
}

int AlignmentSigs::getEndMost()
{
	if (EndMost==-1)
	{
		for (int i=0;i<TypeSignatures.size();++i)
		for (int j=0;j<TypeSignatures[i].size();++j)
		{
			TypeEndMost[i]=max(TypeEndMost[i],TypeSignatures[i][j].End);
		}
		EndMost=max(TypeEndMost[0],TypeEndMost[1]);
	}
	return EndMost;
}

AlignmentSigs::operator std::string()
{
	std::string Result=std::to_string(AlignmentID)+": "+TemplateName+", ("+std::to_string(BeginMost)+", "+std::to_string(EndMost)+")";
	for (int t=0;t<TypeSignatures.size();++t)
	{
		Result+="\n\tt"+std::to_string(t)+"("+std::to_string(TypeBeginMost[t])+", "+std::to_string(TypeEndMost[t])+")";
		for (int i=0;i<TypeSignatures[t].size();++i)
		{
			Result+="\n\t\t"+std::string(TypeSignatures[t][i]);
		}
	}
	return Result;
}

//to improve performance
void divideASs(vector<AlignmentSigs> &ASs)
{
    int OriginalSize=ASs.size();
    for (int i=0;i<OriginalSize;++i)
    {
        for (int t=0;t<ASs[i].TypeSignatures.size();++t)
        {
            if (ASs[i].TypeSignatures[t].size()==0) continue;
            for (int j=0;j<ASs[i].TypeSignatures[t].size()-1;++j)
            {
                if (ASs[i].TypeSignatures[t][j].End+1000<ASs[i].TypeSignatures[t][j+1].Begin)
                {
                    ASs.push_back(AlignmentSigs(ASs[i].AlignmentID));
			        ASs[ASs.size()-1].TypeSignatures[t].insert(ASs[ASs.size()-1].TypeSignatures[t].end(),make_move_iterator(ASs[i].TypeSignatures[t].begin()),make_move_iterator(ASs[i].TypeSignatures[t].begin()+j+1));
                }
            }
        }
    }
}

///=====Omni Merge Ends======

void simpleMergeSigs(int Type, vector<Signature> &Signatures, AlignmentSigs & AlignmentSigs, Arguments & Args, bool Regional=false)
{
    int CurrentBegin=-1, CurrentLength=0;
    double MaxMergeDisPortion=Args.DelMaxMergePortion;
    int MinMaxMergeDis=Args.DelMinMaxMergeDis;
	if (Type==1) MinMaxMergeDis=Args.InsMinMaxMergeDis;
	string InsBases="";
    for (int i=0;i<AlignmentSigs.TypeSignatures[Type].size();++i)
    {
		int ThisLength=AlignmentSigs.TypeSignatures[Type][i].Length;
		if (Type==1)
		{
			ThisLength=AlignmentSigs.TypeSignatures[Type][i].InsBases.size();
		}
        if (CurrentBegin==-1)
        {
            CurrentBegin=AlignmentSigs.TypeSignatures[Type][i].Begin;
			CurrentLength=ThisLength;
			InsBases=AlignmentSigs.TypeSignatures[Type][i].InsBases;
            continue;
        }
        int Begin=AlignmentSigs.TypeSignatures[Type][i].Begin;
		int CalcLength=CurrentLength;
		if (Regional) CalcLength=0;//If regional, compare with the distance from start not last
        if (Begin-CurrentBegin-CalcLength>=(CurrentLength*MaxMergeDisPortion>MinMaxMergeDis?CurrentLength*MaxMergeDisPortion:MinMaxMergeDis))
        {
            if(CurrentLength>=Args.MinSVLen)
            {
                Signature Temp(0,AlignmentSigs.TypeSignatures[Type][0].Tech,Type,CurrentBegin,CurrentBegin+CurrentLength,AlignmentSigs.TemplateName.c_str(),AlignmentSigs.TypeSignatures[Type][0].Quality,InsBases.c_str());
                Signatures.push_back(Temp);
            }
            CurrentBegin=Begin;
			CurrentLength=ThisLength;
			InsBases=AlignmentSigs.TypeSignatures[Type][i].InsBases;
        }
        else
        {
			CurrentLength+=ThisLength;
			InsBases+=AlignmentSigs.TypeSignatures[Type][i].InsBases;
        }
    }
	if (CurrentBegin!=-1)
	{
		if(CurrentLength>=Args.MinSVLen)
		{
                Signature Temp(0,AlignmentSigs.TypeSignatures[Type][0].Tech,Type,CurrentBegin,CurrentBegin+CurrentLength,AlignmentSigs.TemplateName.c_str(),AlignmentSigs.TypeSignatures[Type][0].Quality,InsBases.c_str());
                Signatures.push_back(Temp);
		}
	}
}

void simpleMergeType(int Type, vector<Signature> &Signatures, AlignmentSigs & AlignmentSigs, Arguments & Args, bool Regional)
{
    simpleMergeSigs(Type, Signatures,AlignmentSigs,Args,Regional);
}

void simpleMerge(vector<vector<Signature>> * pTypeSignatures, int Index, vector<AlignmentSigs> * pAlignmentsSigs, Arguments *pArgs, bool Regional)
{
	simpleMergeType(0,pTypeSignatures->at(0),(*pAlignmentsSigs)[Index],*pArgs, Regional);
	simpleMergeType(1,pTypeSignatures->at(1),(*pAlignmentsSigs)[Index],*pArgs, Regional);
}

void *simpleMergeHandler(void * Args)
{
	SimpleMergeArgs * A=(SimpleMergeArgs*)Args;
	simpleMerge(&(A->pTypeSignatures->at(pthread_self())),A->Index,A->pAlignmentsSigs,A->pArgs,A->Regional);
	delete A;
	return NULL;
}
