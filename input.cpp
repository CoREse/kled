#include "htslib/htslib/sam.h"
#include "htslib/htslib/faidx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>
#include <string>
#include <mutex>
#include <unistd.h>
// #include <omp.h>
#include "htslib/htslib/bgzf.h"
#include <thread>
#include <stdexcept>
#include "input.h"
#include "signature.h"
#include <filesystem>
#include "contig.h"
#include "crelib/crelib.h"
#include "StatsManager.hpp"
#include <math.h>
#include <algorithm>
#include "htslib/htslib/thread_pool.h"
#include <set>
#include "merge.h"
#include <unordered_map>
#include "htslib/thread_pool_internal.h"

using namespace std;
using namespace cre;

const int ReadThreadN=1;//old mechanism that leaves this number of threads handling reads. not good. obsoleted. use htslib's thread mechanisms instead

float MedianInsertionSize=481.2;
float ISSD=36.4;
//bool DynamicDRP=false;
//int DynamicCount=1000;
int RDWindowSize=100;

mutex TopLock, WriteLock, StdinLock;

struct HandleBrMutex
{
	pthread_mutex_t m_AllPrimarySeg;
	pthread_mutex_t m_Cov;
	pthread_mutex_t m_Sig[NumberOfSVTypes];
	pthread_mutex_t m_AlignmentsSigs;
	HandleBrMutex()
	{
		pthread_mutex_init(&m_AllPrimarySeg,NULL);
		pthread_mutex_init(&m_Cov,NULL);
		for (int i=0;i<NumberOfSVTypes;++i)
		{
			pthread_mutex_init(m_Sig+i,NULL);
		}
		pthread_mutex_init(&m_AlignmentsSigs,NULL);
	}
};

struct CoverageWindowMutex
{
	pthread_mutex_t m_Cov;
	CoverageWindowMutex()
	{
		pthread_mutex_init(&m_Cov,NULL);
	}
};

int getTechFromReads(bam_hdr_t *Header, htsFile* SamFile)
{
	int CheckN=10;
	bam1_t *br=bam_init1();
	for (int i=0;i<CheckN;++i)
	{
		sam_read1(SamFile,Header,br);
		if (bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br))>1000) return 0;
	}
	return 1;
}

int getTech(const char * ReferenceFileName, const char * BamFileName)//return 0: SMRT, 1:NGS, -1: others
{
	htsFile* SamFile;//the file of BAM/CRAM
	bam_hdr_t* Header;//header for BAM/CRAM file
	hts_idx_t* BamIndex;
	vector<Stats> AllStats;
	const char * SampleFileName=BamFileName;
	SamFile = hts_open(SampleFileName, "rb");
	BamIndex= sam_index_load(SamFile,SampleFileName);
	if (BamIndex==NULL)
	{
		hts_close(SamFile);
		die("Please index the samfile(%s) first!",SampleFileName);
	}
	//set reference file
	if (NULL != ReferenceFileName)
	{
		char referenceFilenameIndex[128];
		strcpy(referenceFilenameIndex, ReferenceFileName);
		strcat(referenceFilenameIndex, ".fai");
		int ret = hts_set_fai_filename(SamFile, referenceFilenameIndex);
	}
	Header = sam_hdr_read(SamFile);

	int Tech;//0: SMRT, 1:NGS
	kstring_t ks=KS_INITIALIZE;
	ks_resize(&ks,1000);
	char PL[1000];
	sam_hdr_find_line_pos(Header,"RG",0,&ks);
	int splitn;
	int * splitoffsets=ksplit(&ks,'\t',&splitn);
	for (int i=0;i<splitn;++i)
	{
		if (memcmp(ks.s+splitoffsets[i],"PL",2)==0)
		{
			int End=ks_len(&ks);
			if (i<splitn-1) End=splitoffsets[i+1];
			memcpy(PL,ks.s+splitoffsets[i]+3,End-splitoffsets[i]-3);
			PL[End-splitoffsets[i]-3]='\0';
		}
	}
	ks_free(&ks);
	if (strcmp(PL,"ILLUMINA")==0)
	{
		Tech=1;
	}
	else if (strcmp(PL,"PACBIO")==0 || strcmp(PL,"ONT")==0)
	{
		Tech=0;
	}
	else
	{
		Tech=getTechFromReads(Header,SamFile);
	}
	hts_close(SamFile);
	return Tech;
}

vector<int> getAllTechs(Arguments & Args)
{
	vector<int> AllTechs;
	for (int i=0;i<Args.BamFileNames.size();++i)
	{
		const char * SampleFileName=Args.BamFileNames[i];
		AllTechs.push_back(getTech(Args.ReferenceFileName,Args.BamFileNames[i]));
	}
	return AllTechs;
}

//get soft and hard clip length sum
float getClipLength(bam1_t * br)
{
	uint32_t* cigars=bam_get_cigar(br);
	float ClipLength=0;
	for (int i=0;i<br->core.n_cigar;++i)
	{
		if ((cigars[i]&0xf)==4 || (cigars[i]&0xf)==5) ClipLength+=cigars[i]>>4;
	}
	return ClipLength;
}

//get query length without head and tail clips
#define brGetClippedQlen(br) (getClippedQLen(br->core.n_cigar,bam_get_cigar(br)))
//get query length without head and tail clips
int getClippedQLen(uint32_t CIGARN, uint32_t* CIGARD)
{
	int Start=0,N=CIGARN;
	if (bam_cigar_op(CIGARD[0])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[0])==BAM_CHARD_CLIP)
	{
		Start=1;
		--N;
	}
	if (bam_cigar_op(CIGARD[CIGARN-1])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[CIGARN-1])==BAM_CHARD_CLIP)
	{
		--N;
	}
	return bam_cigar2qlen(N,CIGARD+Start);
}

//get read length(with all clipped length)
#define brGetReadLength(br) (getReadLength(br->core.n_cigar,bam_get_cigar(br)))
//get read length(with all clipped length)
int getReadLength(uint32_t CIGARN, uint32_t* CIGARD)//included hard clip
{
	int Start=0,N=CIGARN;
	int Clipped=0;
	if (bam_cigar_op(CIGARD[0])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[0])==BAM_CHARD_CLIP)
	{
		Start=1;
		--N;
		Clipped+=bam_cigar_oplen(CIGARD[0]);
	}
	if (bam_cigar_op(CIGARD[CIGARN-1])==BAM_CSOFT_CLIP || bam_cigar_op(CIGARD[CIGARN-1])==BAM_CHARD_CLIP)
	{
		--N;
		Clipped+=bam_cigar_oplen(CIGARD[CIGARN-1]);
	}
	return bam_cigar2qlen(N,CIGARD+Start)+Clipped;
}

void swap(int &a, int &b)
{
	int c=a;
	a=b;
	b=c;
}

struct Alignment
{
	int Cid;//Contig id of mapping position
	int Pos;
	int Length;//Length in Reference
	int End;
	int Strand;//0: -, 1: +
	int MapQ;
	int NM;
	vector<uint32_t> CIGAR;
	int InnerPos;//for same orient segments, just use innerpos, forward pos if for sorting of mixed orientations. Now is the same as ForwardPos
	int InnerLength;
	int InnerEnd;
	int ForwardPos;//for inner sort
	int ForwardEnd;
	// int ReadLength;
	// int FirstCIGAR;
	Alignment(int Cid, int Pos, int Strand, uint32_t* CIGARD, uint32_t CIGARN, int mapQ, int NM)
	:Cid(Cid), Pos(Pos), Strand(Strand), MapQ(mapQ), NM(NM)
	{
		construct(CIGARD,CIGARN);
	}
	void construct(uint32_t* CIGARD, uint32_t CIGARN)
	{
		Length=bam_cigar2rlen(CIGARN, CIGARD);
		End=Pos+Length;
		// FirstCIGAR=CIGARD[0];
		if (bam_cigar_opchr(*CIGARD)=='H' || bam_cigar_opchr(*CIGARD)=='S') InnerPos=bam_cigar_oplen(*CIGARD);
		else InnerPos=0;
		InnerLength=getClippedQLen(CIGARN,CIGARD);
		InnerEnd=InnerPos+InnerLength;
		// ReadLength=getReadLength(CIGARN,CIGARD);
		if (Strand==0)
		{
			int ReadLength=getReadLength(CIGARN,CIGARD);
			ForwardPos=ReadLength-InnerEnd;
			ForwardEnd=ReadLength-InnerPos;
			InnerPos=ForwardPos;
			InnerEnd=ForwardEnd;
		}
		else
		{
			ForwardPos=InnerPos;
			ForwardEnd=InnerEnd;
		}
		for (int i =0;i<CIGARN;++i) CIGAR.push_back(CIGARD[i]);
	}
	Alignment(int Cid, int Pos, int Strand, const char * CIGARS, int MapQ, int NM)
	:Cid(Cid), Pos(Pos), Strand(Strand), MapQ(MapQ), NM(NM)
	{
		size_t CIGARN=size_t(strlen(CIGARS));//sam_parse_cigar will reallocate if not enough
		uint32_t * ca=(uint32_t*) malloc(CIGARN*sizeof(uint32_t));
		int Processed=sam_parse_cigar(CIGARS,NULL,&ca,&CIGARN);
		if (Processed==-1) throw -1;
		CIGARN=Processed;
		construct(ca,CIGARN);
		free(ca);
	}
	bool operator<(const Alignment &other) const
	{
		return this->ForwardPos<other.ForwardPos;
		// if (this->Cid==other.Cid) return this->ForwardPos<other.ForwardPos;//TODO: a mistake, should not take cid into consideration?
		// return this->Cid<other.Cid;
		// return this->InnerPos<other.InnerPos;
	}
	double getQuality()
	{
		return 1.0-((double)(NM))/((double)InnerLength);//InnerLength-NM to estimate AS
	}
};

string cigar2string(uint32_t* CIGARD, uint32_t CIGARN)
{
	string result="";
	for (int i=0;i<CIGARN;++i)
	{
		result+=to_string(bam_cigar_oplen(CIGARD[i]));
		result+=bam_cigar_opchr(CIGARD[i]);
	}
	return result;
}

void searchDelFromAligns(bam1_t *br, int HP, Contig& TheContig, vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, Arguments & Args)
{
	// vector<Signature> &Signatures=TypeSignatures[TheContig.ID][0];
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<Aligns.size();++j)
		{
			if (Aligns[i].Strand!=Aligns[j].Strand || Aligns[i].Cid != Aligns[j].Cid)
			{
				continue;
				// if (i+1<Aligns.size())
				// {
				// 	if (Aligns[i-1].Strand!=Aligns[i+1].Strand) continue;
				// 	i2=i+1;
				// }
				// else continue;
			}
			if (abs(Aligns[i].Pos-Aligns[j].Pos)>1000000) continue;
			int FormerI=i, LatterI=j;
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			// if (Aligns[i-1].End<Aligns[i].Pos && Aligns[i].Pos-Aligns[i-1].End-(Aligns[i].InnerPos-Aligns[i-1].InnerEnd)>=Args.MinSVLen) Signatures.push_back(Signature(2,Tech,0,Aligns[i-1].End,Aligns[i].Pos,bam_get_qname(br),br->core.qual));
			int PosGap=Aligns[LatterI].Pos-Aligns[FormerI].End, InnerGap=Aligns[j].InnerPos-Aligns[i].InnerEnd;
			int GapDiff=PosGap-InnerGap;
			if (Aligns[FormerI].End<Aligns[LatterI].Pos && GapDiff>=Args.MinSVLen)
			{
				if (PosGap>2*(GapDiff)) break;
				// int End=Aligns[FormerI].End+Aligns[LatterI].Pos-Aligns[FormerI].End-(Aligns[LatterI].InnerPos-Aligns[FormerI].InnerEnd);
				int End=Aligns[LatterI].Pos;
				double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
				Signature TempSignature(2,Tech,0,Aligns[FormerI].End,End,bam_get_qname(br),Quality);
				TempSignature.HP=HP;
				// pthread_mutex_lock(&Mut->m_Sig[0]);
				TypeSignatures[Aligns[i].Cid][0].push_back(TempSignature);
				// pthread_mutex_unlock(&Mut->m_Sig[0]);
			}
			break;//Get the first one only. Prevent duplicate
		}
	}
}

void searchInsFromAligns(bam1_t *br, int HP, Contig& TheContig,vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, Arguments & Args)
{
	// vector<Signature> &Signatures=TypeSignatures[TheContig.ID][1];
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<Aligns.size();++j)
		{
			if (Aligns[i].Cid!=Aligns[j].Cid) continue;
			if (Aligns[i].Strand!=Aligns[j].Strand) continue;
			if (abs(Aligns[i].Pos-Aligns[j].Pos)>1000000) continue;
			int FormerI=i, LatterI=j;
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			// int Tolerance=0;
			// if (Aligns[LatterI].Pos+Tolerance<Aligns[FormerI].End) continue;
			// // if (Aligns[LatterI].Pos>Aligns[FormerI].End+Tolerance) continue;
			int Gap=(Aligns[j].InnerPos-Aligns[i].InnerEnd)-(Aligns[LatterI].Pos-Aligns[FormerI].End);
			if (Aligns[FormerI].End<Aligns[LatterI].Pos+Args.InsClipTolerance && Gap>=Args.MinSVLen)
			{
				double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
				Signature TempSignature(2,Tech,1,(Aligns[FormerI].End+Aligns[LatterI].Pos)/2,MIN(TheContig.Size-1,(Aligns[FormerI].End+Aligns[LatterI].Pos)/2+Gap),bam_get_qname(br),Quality);
				TempSignature.HP=HP;
				// pthread_mutex_lock(&Mut->m_Sig[1]);
				TypeSignatures[Aligns[FormerI].Cid][1].push_back(TempSignature);
				// Signatures.push_back(TempSignature);
				// pthread_mutex_unlock(&Mut->m_Sig[1]);
			}
			break;
		}
	}
}

bool continuous(const Alignment& Former, const Alignment& Latter, unsigned Endurance)
{
	// return true;
	// if(abs(Latter.Pos-Former.End)<Endurance && abs(Latter.InnerPos-Former.InnerEnd)<Endurance) return true;
	if (abs(Latter.Pos-Former.Pos)>1000000) return false;
	if (Latter.Pos+Endurance>=Former.End && Latter.Pos>Former.Pos) return true;
	return false;
}

void searchInvFromAligns(bam1_t *br, int HP, Contig& TheContig,vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, Arguments & Args)
{
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<=i+1;++j)
		{
			if (Aligns[i].Strand==Aligns[j].Strand) continue;
			if (Aligns[i].Cid!=Aligns[j].Cid) continue;
			int InvClipEndurance=1000;
			int FormerI=i, LatterI=j;//k is always the middle one, former and latter is for ref position
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			if (continuous(Aligns[FormerI],Aligns[LatterI],InvClipEndurance))
			{
				int k=0;
				for (k=j+1;k<Aligns.size();++k)
				{
					if (Aligns[i].Cid!=Aligns[k].Cid) continue;
					if (Aligns[i].Strand==0)
					{
						if (Aligns[j].Strand!=Aligns[k].Strand && continuous(Aligns[k],Aligns[j],InvClipEndurance)) break;
					}
					else
					{
						if (Aligns[j].Strand!=Aligns[k].Strand && continuous(Aligns[j],Aligns[k],InvClipEndurance)) break;
					}
				}
				if (k!=Aligns.size())
				{
					double Quality=Aligns[i].getQuality()+Aligns[j].getQuality()+Aligns[k].getQuality();
					Quality/=3.0;
					Signature TempSignature(2,Tech,3,Aligns[j].Pos,Aligns[j].End,bam_get_qname(br),Quality);
					TempSignature.HP=HP;
					TempSignature.setInvLeft(true);
					TempSignature.setInvRight(true);
					// pthread_mutex_lock(&Mut->m_Sig[3]);
					TypeSignatures[Aligns[i].Cid][3].push_back(TempSignature);
					// pthread_mutex_unlock(&Mut->m_Sig[3]);
					++i;
				}
				else
				{
					double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
					Signature Temp1Signature(2,Tech,3,Aligns[FormerI].Pos,Aligns[FormerI].End,bam_get_qname(br),Quality);
					Temp1Signature.HP=HP;
					Temp1Signature.setInvRight(true);
					Signature Temp2Signature(2,Tech,3,Aligns[LatterI].Pos,Aligns[LatterI].End,bam_get_qname(br),Quality);
					Temp2Signature.HP=HP;
					Temp2Signature.setInvLeft(true);
					// pthread_mutex_lock(&Mut->m_Sig[3]);
					TypeSignatures[Aligns[FormerI].Cid][3].push_back(Temp1Signature);
					TypeSignatures[Aligns[FormerI].Cid][3].push_back(Temp2Signature);
					// pthread_mutex_unlock(&Mut->m_Sig[3]);
				}
				break;
			}
		}
	}
}

inline void statCoverage(int Begin, int End, float *CoverageWindows, unsigned long ContigSize, CoverageWindowMutex *Mut, Arguments &Args, float Value=1.0)
{
	if (CoverageWindows==NULL) return;
	if (End> ContigSize) End=ContigSize-1;
	if (End<=Begin) return;
	int WBegin=Begin/Args.CoverageWindowSize;
	int WEnd=End/Args.CoverageWindowSize+1;
	double FirstPortion=((double)(((WBegin+1)*Args.CoverageWindowSize)-Begin))/((double)Args.CoverageWindowSize);
	if (Value>=0)
	{
		pthread_mutex_lock(&Mut->m_Cov);
		for (int i=WBegin+1;i<WEnd-1;++i) CoverageWindows[i]+=Value;
		CoverageWindows[WBegin]+=FirstPortion*Value;
		pthread_mutex_unlock(&Mut->m_Cov);
	}
	else
	{
		pthread_mutex_lock(&Mut->m_Cov);
		for (int i=WBegin+1;i<WEnd-1;++i) {CoverageWindows[i]+=Value;if (CoverageWindows[i]<0) CoverageWindows[i]=0;}
		CoverageWindows[WBegin]+=FirstPortion*Value;
		if (CoverageWindows[WBegin]<0) CoverageWindows[WBegin]=0;
		pthread_mutex_unlock(&Mut->m_Cov);
	}
	if (WEnd>WBegin+1)
	{
		double LastPortion=((double)(End+1-(WEnd-1)*Args.CoverageWindowSize))/((double)Args.CoverageWindowSize);
		pthread_mutex_lock(&Mut->m_Cov);
		CoverageWindows[WEnd-1]+=LastPortion*Value;
		if (CoverageWindows[WEnd-1]<0) CoverageWindows[WEnd-1]=0;
		pthread_mutex_unlock(&Mut->m_Cov);
	}
}
inline void statCoverageCigar(bam1_t * br, float *CoverageWindows, Contig & TheContig, CoverageWindowMutex *Mut, Arguments &Args)
{
	// printf("%s %d %d\n", bam_get_qname(br),br->core.pos, br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br)));
	if (br->core.qual<Args.MinMappingQuality) return;
	int Begin=br->core.pos;
	int End=Begin+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
	statCoverage(Begin,End, CoverageWindows, TheContig.Size, Mut, Args);
	// if (End>=TheContig.Size) End=TheContig.Size-1;
	// if (End<=Begin) return;
	// int WBegin=Begin/Args.CoverageWindowSize;
	// int WEnd=End/Args.CoverageWindowSize+1;
	// double FirstPortion=((double)(((WBegin+1)*Args.CoverageWindowSize)-Begin))/((double)Args.CoverageWindowSize);
	// pthread_mutex_lock(&Mut->m_Cov);
	// for (int i=WBegin+1;i<WEnd-1;++i) CoverageWindows[i]+=1;
	// CoverageWindows[WBegin]+=FirstPortion;
	// pthread_mutex_unlock(&Mut->m_Cov);
	// if (WEnd>WBegin+1)
	// {
	// 	double LastPortion=((double)(End+1-(WEnd-1)*Args.CoverageWindowSize))/((double)Args.CoverageWindowSize);
	// 	pthread_mutex_lock(&Mut->m_Cov);
	// 	CoverageWindows[WEnd-1]+=LastPortion;
	// 	pthread_mutex_unlock(&Mut->m_Cov);
	// }
}
//May need Skip the first one(primary alignment)
inline void statCoverageAligns(vector<Alignment> &Aligns, vector<float *> &CoverageWindowsPs, vector<unsigned long> &CoverageWindowsNs, int Skip, CoverageWindowMutex *Mut, Arguments &Args)
{
	if (Aligns.size()==0) return;
	sort(Aligns.data()+Skip,Aligns.data()+Aligns.size(),
	[](Alignment &a, Alignment &b){
		if (a.Cid==b.Cid) return a.Pos<b.Pos;
		return a.Cid<b.Cid;});
	int Cid=Aligns[Skip].Cid,Begin=Aligns[Skip].Pos,End=Aligns[Skip].End;
	for (int i=Skip+1;i<Aligns.size();++i)
	{
		if (Aligns[i].Cid==Cid && Aligns[i].Pos<=End)
		{
			End=Aligns[i].End;
		}
		else
		{
			statCoverage(Begin,End, CoverageWindowsPs[Cid], CoverageWindowsNs[Cid], Mut, Args);
			Cid=Aligns[i].Cid;
			Begin=Aligns[i].Pos;
			End=Aligns[i].End;
		}
	}
	statCoverage(Begin,End, CoverageWindowsPs[Cid], CoverageWindowsNs[Cid], Mut, Args);
}
bool conformDup(int S1, int E1, int S2, int E2, double Threshold=0.3)//from svim
{
	// return false;
	int L1=E1-S1, L2=E2-S2;
	double LengthRatio=double(abs(L1-L2))/(double(min(200,max(L1,L2))));
	double PositionRatio=abs(double(E1+S1)/2.0-double(E2+S2)/2.0)/(double(max(L1,L2)));
	if (LengthRatio + PositionRatio<Threshold) return true;
	// if (LengthRatio<Threshold && PositionRatio<Threshold) return true;
	return false;
}

void searchDupFromAligns(bam1_t *br,int HP,Contig& TheContig,vector<Alignment> &Aligns, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, vector<float*> &CoverageWindowsPs, vector<unsigned long> &CoverageWindowsNs, CoverageWindowMutex *Mut, Arguments & Args)
{
	// #define OLD
	#ifndef OLD
	int CurrentStart=-1, CurrentEnd=-1, CurrentN=0, CurrentCid=-1, CurrentStrand=-1;
	double CurrentQuality=0;
	bool Covered=false;
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<Aligns.size();++j)
		{
			if (j>i+1) break;
			if (Aligns[i].Cid!=Aligns[j].Cid) continue;
			if (Aligns[i].Strand!=Aligns[j].Strand) continue;
			if (abs(Aligns[i].Pos-Aligns[j].Pos)>1000000) continue;
			// vector<Segment> v;
			int FormerI=i, LatterI=j;
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			// int InnerGap=Aligns[j].InnerPos-Aligns[i].InnerEnd;
			if (Aligns[LatterI].Pos<Aligns[FormerI].End)
			{
				// int DupEnd=Aligns[FormerI].End+InnerGap;//doesn't make any sense of any perspect
				int DupEnd=Aligns[FormerI].End;
				int DupLength=DupEnd-Aligns[LatterI].Pos;
				if (DupLength>=Args.MinSVLen)// && DupLength<100000)
				{
					if (CurrentStart==-1)
					{
						CurrentStart=Aligns[LatterI].Pos;
						CurrentEnd=DupEnd;
						CurrentN=2;
						Covered=false;
						CurrentQuality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
						if (Aligns[LatterI].End> Aligns[FormerI].Pos) Covered=true;
						CurrentCid=Aligns[FormerI].Cid;
						CurrentStrand=Aligns[FormerI].Strand;
					}
					else
					{
						if (CurrentCid==Aligns[FormerI].Cid && CurrentStrand==Aligns[FormerI].Strand && conformDup(CurrentStart,CurrentEnd, Aligns[LatterI].Pos, DupEnd))
						{
							CurrentStart=((double)CurrentStart)*(double(CurrentN)/(double(CurrentN+1)))+((double)Aligns[LatterI].Pos)*(1.0/(double(CurrentN+1)));
							CurrentEnd=((double)CurrentEnd)*(double(CurrentN)/(double(CurrentN+1)))+((double)DupEnd)*(1.0/(double(CurrentN+1)));
							CurrentQuality=((double)CurrentQuality)*(double(CurrentN)/(double(CurrentN+1)))+0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality())*(1.0/(double(CurrentN+1)));
							if (Aligns[LatterI].End> Aligns[FormerI].Pos) Covered=true;
							++CurrentN;
						}
						else
						{
							if (CurrentEnd-CurrentStart>Args.MinSVLen)
							{
								Signature TempSignature(2,Tech,2,CurrentStart,CurrentEnd,bam_get_qname(br),CurrentQuality);
								TempSignature.HP=HP;
								// TempSignature.setCN(CurrentN+1);//assume the other strand's CN is 1
								TempSignature.setCN(CurrentN);//As defined in VCFv4.4, the info copy number shall be the CN for the allele
								TempSignature.Covered=Covered;
								// statCoverage(CurrentStart,CurrentEnd,CoverageWindowsPs[CurrentCid],CoverageWindowsNs[CurrentCid],Mut,Args,-1);
								// pthread_mutex_lock(&Mut->m_Sig[2]);
								TypeSignatures[CurrentCid][2].push_back(TempSignature);
								// pthread_mutex_unlock(&Mut->m_Sig[2]);
							}
							CurrentStart=Aligns[LatterI].Pos;
							CurrentEnd=DupEnd;
							CurrentN=2;
							Covered=false;
							CurrentQuality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
							if (Aligns[LatterI].End> Aligns[FormerI].Pos) Covered=true;
							CurrentCid=Aligns[FormerI].Cid;
							CurrentStrand=Aligns[FormerI].Strand;
						}
					}
				}
				break;
				// //   |=========>|
				// //|===>
				// if (DupLength>Aligns[LatterI].Length)
				// {
				// 	Dup=1;
				// }
				// //   |====>/
				// // /============>
				// else if (DupLength>Aligns[FormerI].Length)
				// {
				// 	Dup=1;
				// }
				// //   |====>/
				// // /=====>
				// else
				// {
				// 	if (Aligns[FormerI].End-Aligns[LatterI].Pos>=Args.MinSVLen)
				// 	{
				// 		Dup=1;
				// 	}
				// }
				// if (Dup)
				// {
				// 	double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
				// 	// v.push_back(Segment(Aligns[i-1].Pos,Aligns[i-1].End));
				// 	// v.push_back(Segment(Aligns[i].Pos,Aligns[i].End));
				// 	Signature TempSignature(2,Tech,2,Aligns[LatterI].Pos,Aligns[FormerI].End+InnerGap,bam_get_qname(br),Quality);
				// 	pthread_mutex_lock(&Mut->m_Sig[2]);
				// 	TypeSignatures[Aligns[FormerI].Cid][2].push_back(TempSignature);
				// 	pthread_mutex_unlock(&Mut->m_Sig[2]);
				// }
			}
			// break;
			//if (Aligns[i-1].End<Aligns[i].Pos && Aligns[i].Pos-Aligns[i-1].End-(Aligns[i].InnerPos-Aligns[i-1].InnerEnd)>=50) Signatures.push_back(Signature(2,Tech,2,Aligns[i-1].End,Aligns[i].Pos,bam_get_qname(br),br->core.qual));
		}
		// break;
	}
	if (CurrentStart!=-1)
	{
		if (CurrentEnd-CurrentStart>=Args.MinSVLen)
		{
			Signature TempSignature(2,Tech,2,CurrentStart,CurrentEnd,bam_get_qname(br),CurrentQuality);
			TempSignature.HP=HP;
			TempSignature.setCN(CurrentN);
			TempSignature.Covered=Covered;
			// statCoverage(CurrentStart,CurrentEnd,CoverageWindowsPs[CurrentCid],CoverageWindowsNs[CurrentCid],Mut,Args,-1);
			// pthread_mutex_lock(&Mut->m_Sig[2]);
			TypeSignatures[CurrentCid][2].push_back(TempSignature);
			// pthread_mutex_unlock(&Mut->m_Sig[2]);
		}
	}
	#endif
	#ifdef OLD
	//Break at first sig.
	for (int i=0;i<Aligns.size()-1;++i)
	{
		for (int j=i+1;j<Aligns.size();++j)
		{
			if (Aligns[i].Cid!=Aligns[j].Cid) continue;
			if (Aligns[i].Strand!=Aligns[j].Strand) continue;
			if (abs(Aligns[i].Pos-Aligns[j].Pos)>1000000) continue;
			// vector<Segment> v;
			int FormerI=i, LatterI=j;
			if (Aligns[i].Strand==0)
			{
				FormerI=j;
				LatterI=i;
			}
			int InnerGap=Aligns[j].InnerPos-Aligns[i].InnerEnd;
			if (Aligns[LatterI].Pos<Aligns[FormerI].End)
			{
				int DupLength=Aligns[FormerI].End+InnerGap-Aligns[LatterI].Pos;
				int Dup=0;
				if (DupLength>=Args.MinSVLen) Dup=1;
				// //   |=========>|
				// //|===>
				// if (DupLength>Aligns[LatterI].Length)
				// {
				// 	Dup=1;
				// }
				// //   |====>/
				// // /============>
				// else if (DupLength>Aligns[FormerI].Length)
				// {
				// 	Dup=1;
				// }
				// //   |====>/
				// // /=====>
				// else
				// {
				// 	if (Aligns[FormerI].End-Aligns[LatterI].Pos>=Args.MinSVLen)
				// 	{
				// 		Dup=1;
				// 	}
				// }
				if (Dup)
				{
					double Quality=0.5*(Aligns[FormerI].getQuality()+Aligns[LatterI].getQuality());
					// v.push_back(Segment(Aligns[i-1].Pos,Aligns[i-1].End));
					// v.push_back(Segment(Aligns[i].Pos,Aligns[i].End));
					Signature TempSignature(2,Tech,2,Aligns[LatterI].Pos,Aligns[FormerI].End+InnerGap,bam_get_qname(br),Quality);
					pthread_mutex_lock(&Mut->m_Sig[2]);
					TypeSignatures[Aligns[FormerI].Cid][2].push_back(TempSignature);
					pthread_mutex_unlock(&Mut->m_Sig[2]);
				}
			}
			//if (Aligns[i-1].End<Aligns[i].Pos && Aligns[i].Pos-Aligns[i-1].End-(Aligns[i].InnerPos-Aligns[i-1].InnerEnd)>=50) Signatures.push_back(Signature(2,Tech,2,Aligns[i-1].End,Aligns[i].Pos,bam_get_qname(br),br->core.qual));
		}
		break;
	}
	#endif
}

void dealClipConflicts(vector<Alignment> &Aligns, Arguments & Args)
{
	int Endurance=100;
	for (int i=1;i<Aligns.size();++i)
	{
		if (Aligns[i].Strand!=Aligns[i-1].Strand) continue;
		int FormerI=i-1, LatterI=i;
		if (Aligns[i-1].Strand==0)
		{
			FormerI=i;
			LatterI=i-1;
		}
		if (Aligns[FormerI].InnerEnd>Aligns[LatterI].InnerPos+Endurance)
		{
			int Diff=int(Aligns[FormerI].InnerEnd)-int(Aligns[LatterI].InnerPos);
			Diff/=2;
			Aligns[FormerI].InnerEnd-=Diff;
			Aligns[FormerI].End-=Diff;
			Aligns[FormerI].InnerLength-=Diff;
			Aligns[FormerI].Length-=Diff;
			Aligns[LatterI].InnerPos+=Diff;
			Aligns[LatterI].Pos+=Diff;
			Aligns[LatterI].InnerLength-=Diff;
			Aligns[LatterI].Length-=Diff;
		}
	}
}

void searchForClipSignatures(bam1_t *br, int HP, Contig & TheContig, Sam &SamFile, int Tech, vector<vector<vector<Signature>>> &TypeSignatures, vector<float *>&CoverageWindowsPs, vector<unsigned long> &CoverageWindowsNs, CoverageWindowMutex *Mut, Arguments & Args)
{
	if (!align_is_primary(br)) return;
	// int AS=bam_aux2i(bam_aux_get(br,"AS"));
	// int NM=bam_aux2i(bam_aux_get(br,"NM"));
	// int qlen=bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br));
	// int qlenc=brGetClippedQlen(br);//no clip qlen
	// if (bam_cigar_op(bam_get_cigar(br)[0])==BAM_CSOFT_CLIP) qlenc-=bam_cigar_oplen(bam_get_cigar(br)[0]);
	// if (bam_cigar_op(bam_get_cigar(br)[br->core.n_cigar-1])==BAM_CSOFT_CLIP) qlenc-=bam_cigar_oplen(bam_get_cigar(br)[br->core.n_cigar-1]);
	// printf("%d %d %d %d %d\n",AS,NM,qlenc,qlen,bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br)));
	// return;
	if (Tech==0) statCoverageCigar(br,CoverageWindowsPs[TheContig.ID],TheContig, Mut, Args);
	vector<Alignment> Aligns;
	char* SA_tag_char = bam_get_string_tag(br, "SA");
	if(SA_tag_char == NULL)
	{
		// statCoverageCigar(br,CoverageWindows,TheContig,Args);
		return;
	}
	// printf("name:%s; SA:%s;\n",bam_get_qname(br),SA_tag_char);
	int SACharLen=strlen(SA_tag_char);
	char SATag[SACharLen+1];
	strcpy(SATag,SA_tag_char);
	for (int i=0;i<SACharLen;++i) if (SATag[i]==',' || SATag[i]==';') SATag[i]='\0';
	int k=0;
	int pos,strand,cid,mapQ,NM;
	const char * cigars;
	pos=br->core.pos;
	strand= read_is_forward(br)?1:0;
	NM=0;
	// if (string("a03d97b4-835f-4770-a4d1-6ce8af10d71e")==bam_get_qname(br))
	// {
	// 	fprintf(stderr,"here\n");
	// }
	uint8_t * NMP=bam_aux_get(br,"NM");
	if (NMP!=NULL) NM=bam_aux2i(NMP);
	Aligns.push_back(Alignment(br->core.tid,pos,strand,bam_get_cigar(br),br->core.n_cigar,br->core.qual,NM));
	for (int i=0;i<SACharLen;i+=strlen(SATag+i)+1)
	{
		if (k%6==0)//rname
		{
			if (Args.CallByContig && strcmp(SATag+i,SamFile.Header->target_name[br->core.tid])!=0)//Pass other contig. Should be altered if want to do multi-chromosome sv.
			{
				++k;
				// i+=strlen(SATag+i)+1;
				for (;k%6!=0;++k) i+=strlen(SATag+i)+1;
				continue;
			}
			cid=Args.ContigNameID[SATag+i];
		}
		else if (k%6==1) pos=atoi(SATag+i);//pos
		else if (k%6==2) strand=SATag[i]=='+'?1:0;//strand
		else if (k%6==3) cigars=SATag+i;//CIGAR
		else if (k%6==4) mapQ=atoi(SATag+i);//mapQ
		else if (k%6==5)
		{
			NM=atoi(SATag+i);
			// printf(" %s",cigars);
			// if (abs(pos-br->core.pos)<=1000000)
			Aligns.push_back(Alignment(cid,pos,strand,cigars,mapQ,NM));
		}
		k+=1;
	}
	for (int i=1;i<Aligns.size();++i)
	{
		// statCoverage(Aligns[i].Pos, Aligns[i].End, CoverageWindows, TheContig, Mut, Args);//Still no other contigs.
		// statCoverage(Aligns[i].Pos, Aligns[i].End, CoverageWindowsPs[TheContig.ID], TheContig.Size, Mut, Args);//Still no other contigs.
		if (Aligns[i].Cid==TheContig.ID) statCoverage(Aligns[i].Pos, Aligns[i].End, CoverageWindowsPs[TheContig.ID], TheContig.Size, Mut, Args);//Still no other contigs.
		// statCoverage(Aligns[i].Pos, Aligns[i].End, CoverageWindowsPs[Aligns[i].Cid], CoverageWindowsNs[Aligns[i].Cid], Mut, Args);//Still no other contigs.
		// getDelFromCigar(Aligns[i].CIGAR.data(), Aligns[i].CIGAR.size(),Aligns[i].Pos, bam_get_qname(br), Tech, TypeSignatures[0], Args);
	}
	// statCoverageAligns(Aligns, CoverageWindowsPs, TheContig, 1, Mut, Args);
	sort(Aligns.data(),Aligns.data()+Aligns.size());
	// dealClipConflicts(Aligns,Args);//Careful use. Good for ont but not for sims.
	// pthread_mutex_lock(&Mut->m_Sig[1]);
	// printf("name:%s; size:%lu; SA:%s;",bam_get_qname(br),Aligns.size(),bam_get_string_tag(br, "SA"));
	// // printf("name:%s; %d,%d,%d;",bam_get_qname(br),Aligns[0].Strand,Aligns[0].Pos,Aligns[0].Length);
	// for (int i=0;i<Aligns.size();++i){
	// 	printf(" %d,%d,%d,%d,%d,%d,%d,CIGAR:",Aligns[i].Strand,Aligns[i].Pos,Aligns[i].Length,Aligns[i].InnerPos,Aligns[i].InnerEnd,Aligns[i].ForwardPos,Aligns[i].ForwardEnd);
	// 	if (i!=0) for (int j=0;j<Aligns[i].CIGAR.size();++j) printf("%lu%c,",bam_cigar_oplen(Aligns[i].CIGAR[j]),bam_cigar_opchr(Aligns[i].CIGAR[j]));
	// 	}
	// // for (int i=0;i<Aligns.size();++i) printf(" %d,%d,%d",Aligns[i].Strand,Aligns[i].ForwardPos,Aligns[i].ForwardEnd);
	// printf("\n");
	// pthread_mutex_unlock(&Mut->m_Sig[1]);
	searchDelFromAligns(br,HP,TheContig,Aligns,Tech,TypeSignatures, Args);
	searchInsFromAligns(br,HP,TheContig,Aligns,Tech,TypeSignatures, Args);
	searchDupFromAligns(br,HP,TheContig,Aligns,Tech,TypeSignatures, CoverageWindowsPs, CoverageWindowsNs, Mut, Args);
	searchInvFromAligns(br,HP,TheContig,Aligns,Tech,TypeSignatures, Args);
}

//This kind of signature should - some normal isize when calc svlen
void getDRPSignature(bam1_t * br, const Stats& SampleStats, Contig& TheContig, vector<vector<vector<Signature>>> &TypeSignatures)
{
	// const char * qname=bam_get_qname(br);
	// if (strcmp(qname, "HISEQ1:93:H2YHMBCXX:2:1112:15493:51859")==0)
	// {
	// 	string CIGAR=cigar2string(bam_get_cigar(br),br->core.n_cigar);
	// 	int a=1;
	// }
	if ((!read_is_unmapped(br)) && (!read_mate_is_unmapped(br)) && br->core.mtid==br->core.tid && br->core.isize>=0)
	{
		int ForwardBase=-1, ReverseBase, ForwardEnd, ReverseEnd, DupStart, DupEnd, DupSize;// ForwardBase|-->ForwardEnd  ReverseEnd<--|ReverseBase
		// if (read_is_read1(br) && read_is_forward(br) && (!read_mate_is_forward(br)) || !(read_is_read1(br)) && (!read_is_forward(br)) && read_mate_is_forward(br))
		if (read_is_forward(br) && (!read_mate_is_forward(br)))
		{
			ForwardBase=br->core.pos;
			ForwardEnd=br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
			ReverseBase=ForwardBase+br->core.isize;
			ReverseEnd=br->core.mpos;
			DupStart=ForwardBase;
			DupEnd=ReverseBase;
			DupSize=SampleStats.Mean-(DupEnd-DupStart);
			if (DupSize<0) DupSize=0;
			// uint32_t * cigars=bam_get_cigar(br);
			// int isize=br->core.mpos+bam_cigar2qlen(br->core.n_cigar,cigars)-int(br->core.pos);
			// if (br->core.isize>SampleStats.Mean+3*SampleStats.SD) TypeSignatures[0].push_back(Signature(1,1,0,br->core.pos,br->core.pos+br->core.isize,bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize-SampleStats.Mean));
			// else if (br->core.isize<SampleStats.Mean-3*SampleStats.SD && br->core.isize>-1000) TypeSignatures[2].push_back(Signature(1,1,2,br->core.mpos,bam_endpos(br),bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize>0?SampleStats.Mean-br->core.isize:abs(br->core.isize)+brGetClippedQlen(br)));
			// if (isize>SampleStats.Mean+3*SampleStats.SD) TypeSignatures[0].push_back(Signature(1,1,0,br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br)),br->core.mpos,bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize-SampleStats.Mean));
			// else if (isize<SampleStats.Mean-3*SampleStats.SD && br->core.isize>-1000) TypeSignatures[2].push_back(Signature(1,1,2,br->core.mpos,bam_endpos(br),bam_get_qname(br),br->core.qual,Segment(br->core.pos,bam_endpos(br)),Segment(br->core.mpos,br->core.pos+br->core.isize),br->core.isize>0?SampleStats.Mean-br->core.isize:abs(br->core.isize)+brGetClippedQlen(br)));
		}
		else if ((!read_is_forward(br)) && read_mate_is_forward(br))
		{
			ForwardBase=br->core.mpos;
			ForwardEnd=br->core.pos+br->core.isize;
			ReverseBase=br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
			ReverseEnd=br->core.pos;
			DupStart=ReverseEnd;
			DupEnd=ForwardEnd;
			DupSize=DupEnd-DupStart;
		}
		if (ForwardBase!=-1)
		{
			int BaseGap=ReverseBase-ForwardBase;
			// if (BaseGap>100000)
			// {
			// 	string CIGAR=cigar2string(bam_get_cigar(br),br->core.n_cigar);
			// 	char * qname=bam_get_qname(br);
			// 	int a=1;
			// }
			int DelLength=BaseGap-SampleStats.Mean;
			if (BaseGap>SampleStats.Mean+5*SampleStats.SD)
			{
				// pthread_mutex_lock(&Mut->m_Sig[0]);
				TypeSignatures[TheContig.ID][0].push_back(Signature(1,1,0,ForwardBase+int(SampleStats.Mean/2),ForwardBase+int(SampleStats.Mean/2)+DelLength,bam_get_qname(br),br->core.qual,Segment(ForwardBase,ForwardEnd),Segment(ReverseEnd,ReverseBase),DelLength));
				// pthread_mutex_unlock(&Mut->m_Sig[0]);
			}
			else if (BaseGap<SampleStats.Mean-3*SampleStats.SD)
			{
				// pthread_mutex_lock(&Mut->m_Sig[2]);
				TypeSignatures[TheContig.ID][2].push_back(Signature(1,1,2,DupStart,DupEnd,bam_get_qname(br),br->core.qual,Segment(ForwardBase,ForwardEnd),Segment(ReverseEnd,ReverseBase),DupSize));
				// pthread_mutex_unlock(&Mut->m_Sig[2]);
			}
		}
	}
}
const char * BamBases="NACNGNNNTNNNNNNN";
inline void getInsFromCigar(bam1_t * br, int Tech, vector<Signature>& Signatures, Arguments & Args)
{
	if (br->core.qual<Args.MinMappingQuality) return;
	uint32_t * cigars=bam_get_cigar(br);
	unsigned n_cigar=br->core.n_cigar, pos=br->core.pos;
	const char * qname=bam_get_qname(br);
	int qual=br->core.qual;
// 	return getDelFromCigar(bam_get_cigar(br), br->core.n_cigar, br->core.pos, bam_get_qname(br), br->core.qual, Tech, Signatures, Mut, Args);
	int TLength= bam_cigar2qlen(n_cigar,cigars);
	if (TLength<Args.MinTemplateLength) return;
	int Begin=pos;
	int QueryBegin=0;
	//int MergeDis=500;
	int MinMaxMergeDis=Args.InsMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	float MaxMergeDisPortion=Args.DelMaxMergePortion;
	double MergeScore=0;
	int AS=0;
	uint8_t * ASP=bam_aux_get(br,"AS");
	if (ASP!=NULL) AS=bam_aux2i(ASP);
	int ClippedQLen=brGetClippedQlen(br);
	double Quality=((double)AS)/((double)ClippedQLen);
	int MinSigLen=Args.MinSVLen;
	string Allele="";
	// if (Args.OmniMerge) MinSigLen=10;
	#ifdef DEBUG
	vector<string> MergeStrings;
	int Cumulated=0;
	map<int,int> ItoD;
	int PassedD=0;
	for (int i=0;i<n_cigar;++i)
	{
		if (bam_cigar_op(cigars[i])==BAM_CDEL && bam_cigar_oplen(cigars[i])>=Args.MinSVLen)
		{
			MergeStrings.push_back(to_string(Cumulated)+"M"+to_string(bam_cigar_oplen(cigars[i]))+"D");
			Cumulated=0;
			ItoD[i]=PassedD;
			++PassedD;
		}
		else
		{
			// if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Cumulated+=bam_cigar_oplen(cigars[i]);
			Cumulated+=bam_cigar_oplen(cigars[i]);
		}
	}
	#endif
	for (int i=0;i<n_cigar;++i)
	{
		if (bam_cigar_op(cigars[i])==BAM_CINS && bam_cigar_oplen(cigars[i])>=MinSigLen)
		{
			// int rlen=bam_cigar2rlen(1,cigars+i);
			// int rlen=bam_cigar_oplen(cigars[i]);
			int qlen=bam_cigar_oplen(cigars[i]);
			if(qlen>=MinSigLen)
			{
				Allele="";
				for (int j=0;j<qlen;++j)
				{
					Allele+=BamBases[bam_seqi(bam_get_seq(br),QueryBegin+j)];
				}
				assert(Allele.size()==qlen);
				Signature Temp(0,Tech,1,Begin,Begin+qlen,qname,Quality,Allele.c_str());//for clustering and later processing, the end shall become Begin+qlen(Allele.size())
				// #ifdef DEBUG
				// 	string MergeString="";
				// 	for (int k=0;k<ItoD[BeginI];++k) MergeString+=MergeStrings[k];
				// 	MergeString+="[";
				// 	for (int k=ItoD[BeginI];k<ItoD[i];++k) MergeString+=MergeStrings[k];
				// 	if (MergeString.find("5004M33D10682M59D")!=string::npos)
				// 	{
				// 		MergeString=MergeString;
				// 	}
				// 	MergeString+="]";
				// 	for (int k=ItoD[i];k<MergeStrings.size();++k) MergeString+=MergeStrings[k];
				// 	Temp.setMergeString(MergeString);
				// 	// fprintf(stderr, Temp.MergeString.c_str());
				// #endif
				Signatures.push_back(Temp);
			}
		}
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==1||bam_cigar_op(cigars[i])==4||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) QueryBegin+=bam_cigar_oplen(cigars[i]);
	}
}
// void getInsFromCigar(bam1_t *br, int Tech, vector<Signature>& Signatures, HandleBrMutex *Mut, Arguments & Args)
// {
// 	if (br->core.qual<Args.MinMappingQuality) return;
// 	int TLength= bam_cigar2qlen(br->core.n_cigar,bam_get_cigar(br));
// 	if (TLength<Args.MinTemplateLength) return;
// 	uint32_t * cigars=bam_get_cigar(br);
// 	int CurrentStart=-1, CurrentLength=0;
// 	int Begin=br->core.pos;
// 	int QueryBegin=0;
// 	int CurrentQueryStart=QueryBegin;
// 	string Allele;
// 	//int MergeDis=500;
// 	int MinMaxMergeDis=Args.InsMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
// 	float MaxMergeDisPortion=Args.DelMaxMergePortion;
// 	int AS=0;
// 	uint8_t * ASP=bam_aux_get(br,"AS");
// 	if (ASP!=NULL) AS=bam_aux2i(ASP);
// 	int ClippedQLen=brGetClippedQlen(br);
// 	double Quality=((double)AS)/((double)ClippedQLen);
// 	set<int> IEnds;//recording i of merge ends
// 	set<int> SingleIs;
// 	int Begins[br->core.n_cigar];
// 	int QueryBegins[br->core.n_cigar];
// 	int BeginI=-1;
// 	for (int i=0;i<br->core.n_cigar;++i)
// 	{
// 		Begins[i]=-1;
// 		QueryBegins[i]=-1;
// 	}
// 	for (int i=0;i<br->core.n_cigar;++i)
// 	{
// 		if (bam_cigar_op(cigars[i])==BAM_CINS && bam_cigar_oplen(cigars[i])>=Args.MinSVLen)
// 		{
// 			if (Begins[i]!=-1)
// 			{
// 				Begin=Begins[i];//not the first time
// 				QueryBegin=QueryBegins[i];
// 			}
// 			Begins[i]=Begin;
// 			QueryBegins[i]=QueryBegin;
// 			// int rlen=bam_cigar2rlen(1,cigars+i);
// 			int qlen=bam_cigar_oplen(cigars[i]);
// 			// printf("%d %d %s\n",Begin,rlen,bam_get_qname(br));
// 			string ThisAllele="";
// 			for (int j=0;j<qlen;++j)
// 			{
// 				ThisAllele+=BamBases[bam_seqi(bam_get_seq(br),QueryBegin+j)];
// 			}
// 			if (Args.IndependantMerge && SingleIs.count(i)!=1)
// 			{
// 				Signature Temp(0,Tech,1,Begin,Begin+qlen,bam_get_qname(br),Quality,ThisAllele.c_str());
// 				pthread_mutex_lock(&Mut->m_Sig[1]);
// 				Signatures.push_back(Temp);
// 				pthread_mutex_unlock(&Mut->m_Sig[1]);
// 				SingleIs.insert(i);
// 			}
// 			if (CurrentStart==-1)
// 			{
// 				BeginI=i;
// 				CurrentStart=Begin;
// 				CurrentLength=qlen;
// 				CurrentQueryStart=QueryBegin;
// 				Allele="";
// 			}
// 			else
// 			{
// 				if (Begin-CurrentStart-CurrentLength>=(CurrentLength*MaxMergeDisPortion>MinMaxMergeDis?CurrentLength*MaxMergeDisPortion:MinMaxMergeDis))
// 				{
// 					if(CurrentLength>=Args.MinSVLen)
// 					{
// 						if (!Args.IndependantMerge || ( Args.IndependantMerge && IEnds.count(i)==0))
// 						{
// 							Signature TempSignature(0,Tech,1,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br),Quality,Allele.c_str());
// 							pthread_mutex_lock(&Mut->m_Sig[1]);
// 							Signatures.push_back(TempSignature);
// 							pthread_mutex_unlock(&Mut->m_Sig[1]);
// 						}
// 					}
// 					if (Args.IndependantMerge)
// 					{
// 						i=BeginI;
// 						Begin=Begins[i];
// 						CurrentStart=-1;
// 					}
// 					else
// 					{
// 						CurrentStart=Begin;
// 						CurrentLength=qlen;
// 						CurrentQueryStart=QueryBegin;
// 						Allele="";
// 					}
// 				}
// 				else
// 				{
// 					CurrentLength+=qlen;
// 				}
// 			}
// 			Allele+=ThisAllele;
// 		}
// 		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
// 		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==1||bam_cigar_op(cigars[i])==4||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) QueryBegin+=bam_cigar_oplen(cigars[i]);
// 		//Begin+=bam_cigar2rlen(1,cigars+i);
// 	}
// 	if (CurrentStart!=-1)
// 	{
// 		if(CurrentLength>=Args.MinSVLen)
// 		{
// 			int AS=bam_aux2i(bam_aux_get(br,"AS"));
// 			int ClippedQLen=brGetClippedQlen(br);
// 			double Quality=((double)AS)/((double)ClippedQLen);
// 			Signature TempSignature(0,Tech,1,CurrentStart,CurrentStart+CurrentLength,bam_get_qname(br),Quality,Allele.c_str());
// 			pthread_mutex_lock(&Mut->m_Sig[1]);
// 			Signatures.push_back(TempSignature);
// 			pthread_mutex_unlock(&Mut->m_Sig[1]);
// 		}
// 	}
// }

inline void getDelFromCigar(bam1_t * br, int Tech, vector<Signature>& Signatures, Arguments & Args)
{
	if (br->core.qual<Args.MinMappingQuality) return;
	uint32_t * cigars=bam_get_cigar(br);
	unsigned n_cigar=br->core.n_cigar, pos=br->core.pos;
	const char * qname=bam_get_qname(br);
	int qual=br->core.qual;
// 	return getDelFromCigar(bam_get_cigar(br), br->core.n_cigar, br->core.pos, bam_get_qname(br), br->core.qual, Tech, Signatures, Mut, Args);
	int TLength= bam_cigar2qlen(n_cigar,cigars);
	if (TLength<Args.MinTemplateLength) return;
	int CurrentStart=-1, CurrentLength=0;
	int Begin=pos;
	int LastBegin=Begin;
	//int MergeDis=500;
	int MinMaxMergeDis=Args.DelMinMaxMergeDis;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	float MaxMergeDisPortion=Args.DelMaxMergePortion;
	double MergeScore=0;
	vector<Segment> ShortSigs;
	int AS=0;
	uint8_t * ASP=bam_aux_get(br,"AS");
	if (ASP!=NULL) AS=bam_aux2i(ASP);
	int ClippedQLen=brGetClippedQlen(br);
	double Quality=((double)AS)/((double)ClippedQLen);
	set<int> IEnds;//recording i of merge ends
	set<int> SingleIs;
	//for retrospection
	int Begins[n_cigar];
	int BeginI=-1;
	for (int i=0;i<n_cigar;++i)
	{
		Begins[i]=-1;
	}
	int MinSigLen=Args.MinSVLen;
	// if (Args.OmniMerge) MinSigLen=10;
	#ifdef DEBUG
	vector<string> MergeStrings;
	int Cumulated=0;
	map<int,int> ItoD;
	int PassedD=0;
	for (int i=0;i<n_cigar;++i)
	{
		if (bam_cigar_op(cigars[i])==BAM_CDEL && bam_cigar_oplen(cigars[i])>=Args.MinSVLen)
		{
			MergeStrings.push_back(to_string(Cumulated)+"M"+to_string(bam_cigar_oplen(cigars[i]))+"D");
			Cumulated=0;
			ItoD[i]=PassedD;
			++PassedD;
		}
		else
		{
			// if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Cumulated+=bam_cigar_oplen(cigars[i]);
			Cumulated+=bam_cigar_oplen(cigars[i]);
		}
	}
	#endif
	for (int i=0;i<n_cigar;++i)
	{
		if (bam_cigar_op(cigars[i])==BAM_CDEL && bam_cigar_oplen(cigars[i])>=MinSigLen)
		{
			if (Begins[i]!=-1) Begin=Begins[i];//not the first time
			Begins[i]=Begin;
			// int rlen=bam_cigar2rlen(1,cigars+i);
			int rlen=bam_cigar_oplen(cigars[i]);
			CurrentLength=rlen;
			if(CurrentLength>=MinSigLen)
			{
				Signature Temp(0,Tech,0,Begin,Begin+rlen,qname,Quality);
				#ifdef DEBUG
					string MergeString="";
					for (int k=0;k<ItoD[BeginI];++k) MergeString+=MergeStrings[k];
					MergeString+="[";
					for (int k=ItoD[BeginI];k<ItoD[i];++k) MergeString+=MergeStrings[k];
					if (MergeString.find("5004M33D10682M59D")!=string::npos)
					{
						MergeString=MergeString;
					}
					MergeString+="]";
					for (int k=ItoD[i];k<MergeStrings.size();++k) MergeString+=MergeStrings[k];
					Temp.setMergeString(MergeString);
					// fprintf(stderr, Temp.MergeString.c_str());
				#endif
				Signatures.push_back(Temp);
			}
		}
		if (bam_cigar_op(cigars[i])==0 ||bam_cigar_op(cigars[i])==2||bam_cigar_op(cigars[i])==7||bam_cigar_op(cigars[i])==8) Begin+=bam_cigar_oplen(cigars[i]);
	}
}

struct HandleBrArgs
{
	bam1_t *br;
	Contig *pTheContig;
	Sam *pSamFile;
	int Tech;
	const Stats *pSampleStats;
	unordered_map<pthread_t,vector<AlignmentSigs>> * pAlignmentsSigs;
	unordered_map<pthread_t,vector<vector<vector<Signature>>>> *pTypeSignatures;
	unordered_map<pthread_t,SegmentSet> *pAllPrimarySegments;
	vector<float *> *pCoverageWindowsPs;
	vector<unsigned long>* pCoverageWindowsNs;
	CoverageWindowMutex *mut;
	// HandleBrMutex *mut;
	Arguments * pArgs;
};
// double meanlength=0,tcount=0;

void handlebr(bam1_t *br, Contig * pTheContig, Sam *pSamFile, int Tech, const Stats *pSampleStats, vector<AlignmentSigs> *pAlignmentsSigs, vector<vector<vector<Signature>>> *pTypeSignatures, SegmentSet * pAllPrimarySegments, vector<float*> &CoverageWindowsPs, vector<unsigned long> &CoverageWindowsNs, CoverageWindowMutex *Mut, Arguments * pArgs)
{
	Contig & TheContig=*pTheContig;
	Sam & SamFile=*pSamFile;
	const Stats & SampleStats=*pSampleStats;
	SegmentSet & AllPrimarySegments=*pAllPrimarySegments;
	Arguments & Args=*pArgs;
	if (align_is_primary(br))
	{
		hts_pos_t end=br->core.pos+bam_cigar2rlen(br->core.n_cigar,bam_get_cigar(br));
		// pthread_mutex_lock(&mut->m_AllPrimarySeg);
		AllPrimarySegments.add(br->core.pos,end);
		// pthread_mutex_unlock(&mut->m_AllPrimarySeg);
	}
	int HP=0;
	uint8_t * HPP=bam_aux_get(br,"HP");
	if (HPP!=NULL) HP=bam_aux2i(HPP);
	unsigned long long BrHash=0;
	// BrHash=br->core.tid;
	BrHash=br->l_data;
	BrHash<<=24;
	BrHash+=br->core.pos;
	BrHash<<=16;
	BrHash+=br->core.flag;
	BrHash<<=16;
	BrHash+=br->core.n_cigar&0xffff;
	AlignmentSigs TheAlignment(BrHash, bam_get_qname(br));
	// AlignmentSigs TheAlignment(br->id, bam_get_qname(br));
	getDelFromCigar(br,Tech,TheAlignment.TypeSignatures[0], Args);
	getInsFromCigar(br,Tech,TheAlignment.TypeSignatures[1], Args);
	if (TheAlignment.TypeSignatures[0].size()!=0 || TheAlignment.TypeSignatures[1].size()!=0)
	{
		// pthread_mutex_lock(&mut->m_AlignmentsSigs);
		TheAlignment.HP=HP;
		pAlignmentsSigs->push_back(TheAlignment);
		// pthread_mutex_unlock(&mut->m_AlignmentsSigs);
	}
	if (Tech==1)
	{
		if (read_is_paired(br))
		{
			getDRPSignature(br, SampleStats, TheContig, *pTypeSignatures);
		}
	}
	searchForClipSignatures(br, HP, TheContig, SamFile, Tech, *pTypeSignatures, CoverageWindowsPs, CoverageWindowsNs, Mut, Args);
}

void * handlebrWrapper(void * args)
{
	pthread_t tid=pthread_self();
	HandleBrArgs * A=(HandleBrArgs*)args;
	handlebr(A->br, A->pTheContig, A->pSamFile, A->Tech, A->pSampleStats, &(A->pAlignmentsSigs->at(tid)), &(A->pTypeSignatures->at(tid)), &(A->pAllPrimarySegments->at(tid)), *(A->pCoverageWindowsPs), *(A->pCoverageWindowsNs), A->mut, A->pArgs);
	bam_destroy1(A->br);
	delete A;
	return NULL;
}

void calcMeanSD(StatsManager & SM, const char * SampleFileName, float &Mean, float &SD)
{
	Mean=0,SD=0;
	const int Size=10000;
	for (int i=0;i<=Size;++i)
	{
		Mean+=SM.getInsertLen(SampleFileName,float(i)/float(Size));
	}
	Mean/=float(Size+1);
	for (int i=0;i<=Size;++i)
	{
		float V=SM.getInsertLen(SampleFileName,float(i)/float(Size));
		//fprintf(stderr,"%f\n",V);
		V-=Mean;
		V*=V;
		SD+=V;
	}
	SD/=float(Size+1);
	SD=pow(SD,0.5);
}

Stats getSampleStats(const char * ReferenceFileName, const char * SampleFileName,int Tech=1)
{
		Stats SampleStats={0,0,0,0,0};
		if (Tech!=1) return SampleStats;
		StatsManager SM(ReferenceFileName,"");
		SM.handleBamCramStats(SampleFileName);
		SampleStats.MedianIS=SM.getInsertLen(SampleFileName,0.5f);
		SampleStats.UpIS=SM.getInsertLen(SampleFileName,0.99f);
		SampleStats.BelowIS=SM.getInsertLen(SampleFileName,0.01f);
		calcMeanSD(SM, SampleFileName, SampleStats.Mean, SampleStats.SD);
		return SampleStats;
}

vector<Stats> getAllStats(const char * ReferenceFileName, const vector<const char *> & BamFileNames, vector<int> AllTechs)
{
	vector<Stats> AllStats;
	for (int k=0;k<BamFileNames.size();++k)
	{
		const char * SampleFileName=BamFileNames[k];
		AllStats.push_back(getSampleStats(ReferenceFileName,SampleFileName,AllTechs[k]));
	}
	return AllStats;
}

string getSampleName(bam_hdr_t* Header)
{
	char SampleName[1024];
	kstring_t ks=KS_INITIALIZE;
	ks_resize(&ks,1000);
	sam_hdr_find_line_pos(Header,"RG",0,&ks);
	int splitn;
	int * splitoffsets=ksplit(&ks,'\t',&splitn);
	for (int i=0;i<splitn;++i)
	{
		if (memcmp(ks.s+splitoffsets[i],"SM",2)==0)
		{
			int End=ks_len(&ks);
			if (i<splitn-1) End=splitoffsets[i+1];
			memcpy(SampleName,ks.s+splitoffsets[i]+3,End-splitoffsets[i]-3);
			SampleName[End-splitoffsets[i]-3]='\0';
		}
	}
	ks_free(&ks);
	return string(SampleName);
}
Sam::Sam():SamFile(NULL), Header(NULL), BamIndex(NULL){}
void Sam::close()
{
	if (SamFile!=NULL) hts_close(SamFile);
}
htsThreadPool p = {NULL, 0};
vector<Sam> initSam(Arguments & Args)
{
	vector<Sam> SamFiles;
	if (Args.ThreadN > 1) {
		if (!(p.pool = hts_tpool_init(Args.ThreadN))) {
			die("Error creating thread pool\n");
		}
		p.qsize=Args.ThreadN*2;
	}
	for (int i=0;i<Args.BamFileNames.size();++i)
	{
		SamFiles.push_back(Sam());
		const char * SampleFileName=Args.BamFileNames[i];
		SamFiles[i].SamFile = hts_open(SampleFileName, "rb");
		SamFiles[i].BamIndex= sam_index_load(SamFiles[i].SamFile,SampleFileName);
		if (SamFiles[i].BamIndex==NULL)
		{
			die("Please index the samfile(%s) first!",SampleFileName);
		}
		//set reference file
		const char * ReferenceFileName=Args.ReferenceFileName;
		if (NULL != ReferenceFileName)
		{
			char referenceFilenameIndex[128];
			strcpy(referenceFilenameIndex, ReferenceFileName);
			strcat(referenceFilenameIndex, ".fai");
			int ret = hts_set_fai_filename(SamFiles[i].SamFile, referenceFilenameIndex);
		}
		SamFiles[i].Header = sam_hdr_read(SamFiles[i].SamFile);
		if (Args.SampleName=="*")
		{
			string SampleName=getSampleName(SamFiles[i].Header);
			if (SampleName!="") Args.SampleName=SampleName;
		}
		if (p.pool)
		{
			hts_set_opt(SamFiles[i].SamFile,  HTS_OPT_THREAD_POOL, &p);
			// if (settings.out) hts_set_opt(settings.out, HTS_OPT_THREAD_POOL, &p);
		}
	}
	return SamFiles;
}

void closeSam(vector<Sam> &SamFiles)
{
	for (int i=0;i<SamFiles.size();++i) SamFiles[i].close();
	if (p.pool) hts_tpool_destroy(p.pool);
}

struct ReduceCWArgs
{
	double * CoverageWindows;
	int Begin;
	int End;
	unordered_map<pthread_t,double *> &ProcessCovW;
};
void * reduceCW(void * Args)
{
	ReduceCWArgs *A=(ReduceCWArgs*)Args;
	for (unsigned long i=A->Begin;i<A->End;++i)
	{
		for (int pi=0;pi<A->ProcessCovW.size();++pi)
		A->CoverageWindows[i]+=A->ProcessCovW[(p.pool->t+pi)->tid][i];
	}
	delete A;
	return NULL;
}

// unsigned long long TotalReduceTime=0;
void collectSignatures(Contig &TheContig, vector<vector<vector<Signature>>> &TypeSignatures, SegmentSet & AllPrimarySegments, Arguments & Args, vector<Sam>& SamFiles, const vector<Stats> &AllStats, const vector<int> &AllTechs, vector<float*>& CoverageWindowsPs, vector<unsigned long> &CoverageWindowsNs, const char * DataSource)
{
	// unsigned long long ReduceTime=0;
	const char * ReferenceFileName=Args.ReferenceFileName;
	const vector<const char *> & BamFileNames=Args.BamFileNames;
	vector<AlignmentSigs> AlignmentsSigs;//All CIGAR sigs for this contig is here, so we can proceed them here.
	CoverageWindowMutex CWMutex;
	for (int k=0;k<BamFileNames.size();++k)
	{
		const char * SampleFileName=BamFileNames[k];
		Args.Log.info("Reading %s for region %s...",SampleFileName,TheContig.Name.c_str());
		const Stats &SampleStats=AllStats[k];

		string Region=TheContig.Name;

		int Tech=AllTechs[k];

		int ReadCount=0, UnmappedCount=0;
		if (DataSource!=0)
		{
			Args.Log.error("Deprecated function!");
			exit(100);
		}
		else
		{
			if (Args.ThreadN==1)
			{
				bam1_t *br=bam_init1();
				hts_itr_t* RegionIter=sam_itr_querys(SamFiles[k].BamIndex,SamFiles[k].Header,Region.c_str());
				while(sam_itr_next(SamFiles[k].SamFile, RegionIter, br) >=0)//read record
				{
					handlebr(br,&TheContig, &SamFiles[k], Tech, &SampleStats, &AlignmentsSigs, &TypeSignatures, &AllPrimarySegments, CoverageWindowsPs, CoverageWindowsNs, &CWMutex, &Args);
				}
				bam_destroy1(br);
			}
			else
			{
				unordered_map<pthread_t,vector<AlignmentSigs>> ProcessASVec;
				unordered_map<pthread_t,vector<vector<vector<Signature>>> > ProcessSigVec;
				unordered_map<pthread_t,SegmentSet> ProcessSegS;
				time_t StartR=getTimeInMuse();
				for (int pi=0;pi<p.pool->tsize;++pi)
				{
					ProcessASVec[(p.pool->t+pi)->tid]=vector<AlignmentSigs>();
					ProcessSigVec[(p.pool->t+pi)->tid]=vector<vector<vector<Signature>>>();
					ProcessSegS[(p.pool->t+pi)->tid]=SegmentSet();
					// double * PCW=(double *)malloc(sizeof(double)*CoverageWindowsN);
					// assert(ProcessCovW[(p.pool->t+pi)->tid]==NULL);
					// for (unsigned long i=0;i<CoverageWindowsN;++i) PCW[i]=0;
					// ProcessCovW[(p.pool->t+pi)->tid]=PCW;
					for (int i=0;i<TypeSignatures.size();++i)
					{
						ProcessSigVec[(p.pool->t+pi)->tid].push_back(vector<vector<Signature>>());
						for (int j=0;j<TypeSignatures[i].size();++j)
						{
							ProcessSigVec[(p.pool->t+pi)->tid][i].push_back(vector<Signature>());
						}
					}
				}
				// ReduceTime+=getTimeInMuse()-StartR;
				hts_tpool_process *HandlebrProcess=hts_tpool_process_init(p.pool,p.qsize,1);
				bam1_t *br=bam_init1();
				hts_itr_t* RegionIter=sam_itr_querys(SamFiles[k].BamIndex,SamFiles[k].Header,Region.c_str());
				while(sam_itr_next(SamFiles[k].SamFile, RegionIter, br) >=0)//read record
				{
					bam1_t *cbr=bam_dup1(br);
					HandleBrArgs *A=new HandleBrArgs{cbr,&TheContig, &SamFiles[k], Tech, &SampleStats, &ProcessASVec, &ProcessSigVec, &ProcessSegS, &CoverageWindowsPs, &CoverageWindowsNs, &CWMutex, &Args};
					hts_tpool_dispatch(p.pool,HandlebrProcess,handlebrWrapper,(void *)A);
				}
				bam_destroy1(br);
				hts_tpool_process_flush(HandlebrProcess);
				// StartR=getTimeInMuse();
				for (int pi=0;pi<p.pool->tsize;++pi)
				{
					AlignmentsSigs.insert(AlignmentsSigs.end(),make_move_iterator(ProcessASVec[(p.pool->t+pi)->tid].begin()),make_move_iterator(ProcessASVec[(p.pool->t+pi)->tid].end()));
					for (int i=0;i<TypeSignatures.size();++i)
					{
						for (int j=0;j<TypeSignatures[i].size();++j)
						{
							//make move iterator works even with different life intervals.
							TypeSignatures[i][j].insert(TypeSignatures[i][j].end(),make_move_iterator(ProcessSigVec[(p.pool->t+pi)->tid][i][j].begin()),make_move_iterator(ProcessSigVec[(p.pool->t+pi)->tid][i][j].end()));
						}
					}
					AllPrimarySegments.merge(ProcessSegS[(p.pool->t+pi)->tid]);
					// // #pragma omp parallel for
					// for (unsigned long i=0;i<CoverageWindowsN;++i) CoverageWindows[i]+=ProcessCovW[(p.pool->t+pi)->tid][i];
				}
				hts_tpool_process_destroy(HandlebrProcess);
				// ReduceTime+=getTimeInMuse()-StartR;
			}
		}
		// TotalReduceTime+=ReduceTime;
		// Args.Log.debug("Reduce time for %s is %lfs, total: %lfs.",TheContig.Name.c_str(),double(ReduceTime)/1000000.0,double(TotalReduceTime)/1000000.0);
	}

	if (AlignmentsSigs.size()>0)
	{
		Args.Log.verbose("Number of CIGAR sig alignments:%lu",AlignmentsSigs.size());
		if (Args.Log.Verbosity>1)
		{
			unsigned long DELN=0,INSN=0;
			for (int i=0;i<AlignmentsSigs.size();++i)
			{
				DELN+=AlignmentsSigs[i].TypeSignatures[0].size();
				INSN+=AlignmentsSigs[i].TypeSignatures[1].size();
			}
			Args.Log.debug("Size:%lu, DELN:%lu, INSN:%lu",AlignmentsSigs.size(),DELN,INSN);
		}
		if (Args.CigarMerge==0)
		{
			int *TypeSigIndexes[2];
			TypeSigIndexes[0]=(int *)malloc(AlignmentsSigs.size()*sizeof(int));
			TypeSigIndexes[1]=(int *)malloc(AlignmentsSigs.size()*sizeof(int));
			for (int i=0;i<AlignmentsSigs.size();++i)
			{
				TypeSigIndexes[0][i]=i;
				TypeSigIndexes[1][i]=i;
				AlignmentsSigs[i].getBeginMost();
				AlignmentsSigs[i].getEndMost();
			}
			sort(TypeSigIndexes[0],TypeSigIndexes[0]+AlignmentsSigs.size(),[&AlignmentsSigs](int a, int b) {
				if (AlignmentsSigs[a].TypeBeginMost[0]==AlignmentsSigs[b].TypeBeginMost[0]) if (AlignmentsSigs[a].AlignmentID==AlignmentsSigs[b].AlignmentID) return AlignmentsSigs[a].TemplateName<AlignmentsSigs[b].TemplateName; else return AlignmentsSigs[a].AlignmentID<AlignmentsSigs[b].AlignmentID;
				return AlignmentsSigs[a].TypeBeginMost[0]<AlignmentsSigs[b].TypeBeginMost[0];
			});
			sort(TypeSigIndexes[1],TypeSigIndexes[1]+AlignmentsSigs.size(),[&AlignmentsSigs](int a, int b) {
				if (AlignmentsSigs[a].TypeBeginMost[1]==AlignmentsSigs[b].TypeBeginMost[1]) if (AlignmentsSigs[a].AlignmentID==AlignmentsSigs[b].AlignmentID) return AlignmentsSigs[a].TemplateName<AlignmentsSigs[b].TemplateName; else return AlignmentsSigs[a].AlignmentID<AlignmentsSigs[b].AlignmentID;
				return AlignmentsSigs[a].TypeBeginMost[1]<AlignmentsSigs[b].TypeBeginMost[1];
			});
			int *TypeMaxEnds[2];
			TypeMaxEnds[0]=(int*)malloc(AlignmentsSigs.size()*sizeof(int));
			TypeMaxEnds[1]=(int*)malloc(AlignmentsSigs.size()*sizeof(int));
			for (int T=0;T<2;++T)
			{
				TypeMaxEnds[T][0]=AlignmentsSigs[TypeSigIndexes[T][0]].TypeEndMost[T];
				for (int i=1;i<AlignmentsSigs.size();++i)
				{
					TypeMaxEnds[T][i]=max(TypeMaxEnds[T][i-1],AlignmentsSigs[TypeSigIndexes[T][i]].TypeEndMost[T]);
				}
			}
			// for (int i=0;i<AlignmentsSigs.size();++i) printf("%s\n",string(AlignmentsSigs[TypeSigIndexes[0][i]]).c_str());
			// return;
			Args.Log.info("Start merging for %s...",TheContig.Name.c_str());
			if (Args.ThreadN==1)
			{
				for (int i=0;i<AlignmentsSigs.size();++i)
				{
					omniBMerge(&(TypeSignatures[TheContig.ID]),i,&AlignmentsSigs, (const int **)TypeSigIndexes, (const int **)TypeMaxEnds, &Args);
				}
			}
			else
			{
				unordered_map<pthread_t,vector<vector<Signature>> > ProcessSigVec;
				for (int pi=0;pi<p.pool->tsize;++pi)
				{
					ProcessSigVec[(p.pool->t+pi)->tid]=vector<vector<Signature>>();
					for (int i=0;i<2;++i)
					{
						ProcessSigVec[(p.pool->t+pi)->tid].push_back(vector<Signature>());
					}
				}
				hts_tpool_process *MergingSigProcess=hts_tpool_process_init(p.pool,p.qsize,1);
				for (int i=0;i<AlignmentsSigs.size();++i)
				{
					// omniBMerge(&TypeSignatures[TheContig.ID], &AlignmentsSigs[i], &AlignmentsSigs, MaxEnds, &mut);
					OmniBMergeArgs *A=new OmniBMergeArgs{&ProcessSigVec, i, &AlignmentsSigs, (const int **)TypeSigIndexes, (const int **)TypeMaxEnds, &Args};
					hts_tpool_dispatch(p.pool,MergingSigProcess,omniBHandler,A);
				}
				hts_tpool_process_flush(MergingSigProcess);
				hts_tpool_process_destroy(MergingSigProcess);
				for (int pi=0;pi<p.pool->tsize;++pi)
				{
					for (int i=0;i<2;++i)
					{
							//make move iterator works even with different life intervals.
						TypeSignatures[TheContig.ID][i].insert(TypeSignatures[TheContig.ID][i].end(),make_move_iterator(ProcessSigVec[(p.pool->t+pi)->tid][i].begin()),make_move_iterator(ProcessSigVec[(p.pool->t+pi)->tid][i].end()));
					}
				}
			}
			free(TypeMaxEnds[0]);
			free(TypeMaxEnds[1]);
			free(TypeSigIndexes[0]);
			free(TypeSigIndexes[1]);
		}
		else
		{
			if (Args.ThreadN==1)
			{
				for (int i=0;i<AlignmentsSigs.size();++i)
				{
					simpleMerge(&TypeSignatures[TheContig.ID], i, &AlignmentsSigs, &Args, Args.CigarMerge==1?false:true);
				}
			}
			else
			{
				unordered_map<pthread_t,vector<vector<Signature>> > ProcessSigVec;
				for (int pi=0;pi<p.pool->tsize;++pi)
				{
					ProcessSigVec[(p.pool->t+pi)->tid]=vector<vector<Signature>>();
					for (int i=0;i<2;++i)
					{
						ProcessSigVec[(p.pool->t+pi)->tid].push_back(vector<Signature>());
					}
				}
				hts_tpool_process *MergingSigProcess=hts_tpool_process_init(p.pool,p.qsize,1);
				for (int i=0;i<AlignmentsSigs.size();++i)
				{
					// omniBMerge(&TypeSignatures[TheContig.ID], &AlignmentsSigs[i], &AlignmentsSigs, MaxEnds, &mut);
					SimpleMergeArgs *A=new SimpleMergeArgs{&ProcessSigVec, i, &AlignmentsSigs, &Args, Args.CigarMerge==1?false:true};
					hts_tpool_dispatch(p.pool,MergingSigProcess,simpleMergeHandler,A);
				}
				hts_tpool_process_flush(MergingSigProcess);
				hts_tpool_process_destroy(MergingSigProcess);
				for (int pi=0;pi<p.pool->tsize;++pi)
				{
					for (int i=0;i<2;++i)
					{
							//make move iterator works even with different life intervals.
						TypeSignatures[TheContig.ID][i].insert(TypeSignatures[TheContig.ID][i].end(),make_move_iterator(ProcessSigVec[(p.pool->t+pi)->tid][i].begin()),make_move_iterator(ProcessSigVec[(p.pool->t+pi)->tid][i].end()));
					}
				}
			}
		}
		Args.Log.info("Done merging for %s...",TheContig.Name.c_str());
	}
}

Contig * getContigs(Arguments & Args, int& NSeq, int RDWindowSize)
{
	faidx_t * Ref=fai_load(Args.ReferenceFileName);
	NSeq=faidx_nseq(Ref);
	Contig * Contigs=(Contig*) malloc(sizeof(Contig)*(NSeq));
	for (int i=0;i<(NSeq);++i)
	{
		const char * ContigName=faidx_iseq(Ref,i);
		int SeqLen=faidx_seq_len(Ref,ContigName);
		new (Contigs+i) Contig(i,ContigName, SeqLen);
		Args.ContigNameID[ContigName]=i;
	}
	fai_destroy(Ref);
	return Contigs;
}