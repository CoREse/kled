#pragma once
#ifndef KLED_KLED_H
#define KLED_KLED_H

// #define DEBUG

#define COMMA ,

#ifdef DEBUG
#define DEBUG_CODE(Code) Code
#else
#define DEBUG_CODE(Code) 
#endif

#include <vector>
#include <map>
#include <string>
#include "crelib/crelib.h"
#include "optutils/OptHelper.h"

const int NumberOfSVTypes=4;//Default is static so is fine.
const char * const SVTypeNames[NumberOfSVTypes]={"DEL","INS","DUP","INV"};//Default is static so is fine.

struct Arguments {
	const char * Version="1.2.9H11";
	bool ShowVersion=false;
	bool ShowHelp=false;
	cre::Logger Log;
	std::string CommandLine;
	std::string RunHash;//Storing std::hash<std::string>(Version+Arguments) in hex.
	OptHelper * OH=NULL;
	int TestN=0;
	const char * ReferenceFileName=0;
	std::vector<const char *> BamFileNames;
	std::vector<const char *> CallingContigs;
	std::map<std::string,int> ContigNameID;
	std::string SampleName="*";
	int ThreadN=8;
	bool NoFilter=false;
	bool Filter2ST=false;
	int MinSVLen=30;
	int MinMappingQuality=20;
	int MinTemplateLength=500;
	int DelMinMaxMergeDis=500;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	int DelMaxMaxMergeDis=50000;
	int InsMinMaxMergeDis=500;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	double DelMaxMergePortion=0.2;
	int CoverageWindowSize=100;
	int InsClipTolerance=10;
	int BreakerLength=50;//It doesn't make any sense if there is a insertion detected within a deletion, there for merge breaks when meet a long insertion
	int InsMaxGapSize=50000;
	int MaxClusterSize=1000;//will limit (sampling) Cluster size to [MaxClusterSize, 2*MaxClusterSize)
	int ClusteringMaxMergeRange=20000;
	int ClusteringBatchSize=10000;
    double BrotherhoodTypeRatios[NumberOfSVTypes]={1.5,0.5,0.2,0.1};//For SV Types
    int BrotherhoodTypeForceBrothers[NumberOfSVTypes]={50,150,0,100};
    double BrotherhoodTypeLengthRatios[NumberOfSVTypes]={0.3,0.1,0.1,0.1};//For SV Types
    int BrotherhoodTypeLengthMinEndurance[NumberOfSVTypes]={10,20,0,0};
    int BrotherhoodNearRanges[NumberOfSVTypes]={-1,-1,-1,-1};//If -1, use brotherhoodcluster1
    int BrotherhoodTypeForceBrothers2[NumberOfSVTypes]={10,50,500,100};
    double BrotherhoodTypeLengthRatios2[NumberOfSVTypes]={0.5,0.3,0.1,0.5};//For SV Types
    double BrotherhoodCLRTypeRatios[NumberOfSVTypes]={1.5,1.2,0.2,0.1};//For SV Types
       int BrotherhoodCLRTypeForceBrothers[NumberOfSVTypes]={50,10,0,100};
    double BrotherhoodCLRTypeLengthRatios[NumberOfSVTypes]={0.3,0.3,0.1,0.1};//For SV Types
       int BrotherhoodCLRTypeLengthMinEndurance[NumberOfSVTypes]={60,0,0,0};
       int BrotherhoodCLRNearRanges[NumberOfSVTypes]={-1,-1,-1,-1};//If -1, use brotherhoodcluster1
       int BrotherhoodCLRTypeForceBrothers2[NumberOfSVTypes]={10,50,500,100};
    double BrotherhoodCLRTypeLengthRatios2[NumberOfSVTypes]={0.5,0.3,0.1,0.5};//For SV Types
    double BrotherhoodCCSTypeRatios[NumberOfSVTypes]={2.0,2.5,0.2,0.1};//For SV Types
    int BrotherhoodCCSTypeForceBrothers[NumberOfSVTypes]={100,200,0,100};
    double BrotherhoodCCSTypeLengthRatios[NumberOfSVTypes]={0.5,0.5,0.1,0.1};//For SV Types
    int BrotherhoodCCSTypeLengthMinEndurance[NumberOfSVTypes]={10,50,0,0};
       int BrotherhoodCCSNearRanges[NumberOfSVTypes]={-1,-1,-1,-1};//If -1, use brotherhoodcluster1
       int BrotherhoodCCSTypeForceBrothers2[NumberOfSVTypes]={10,50,500,100};
    double BrotherhoodCCSTypeLengthRatios2[NumberOfSVTypes]={0.5,0.3,0.1,0.5};//For SV Types
	bool AllCLR=false;
	bool AllCCS=false;
	double ASSBases[NumberOfSVTypes][2]=//Layers of base filter of the addition of supporting segmentations and templates
	{{1,0}//DEL
	,{1,0}//INS
	,{0,1}//DUP
	,{2,1}};//INV
	double ASSCoverageMulti[NumberOfSVTypes][2]={{0.28,0.16},
	{0.19,0.14},
	{0.23,0.12},
	{0,0.18}};
	double LSDRSs[NumberOfSVTypes][2]={{70, 95},//For Legnth Standard Deviation Ratio Scores(100-ratio*100)
	{60,90},
	{90,96},
	{98,0}};
	double CLRASSBases[NumberOfSVTypes][2]=//Layers of base filter of the addition of supporting segmentations and templates
	{{1,0}//DEL
	,{2,0}//INS
	,{4,0}//DUP
	,{2,1}};//INV
	double CLRASSCoverageMulti[NumberOfSVTypes][2]={{0.13,0.12},
	{0.12,0.5},
	{0.13,0.1},
	{0,0.18}};
	double CLRLSDRSs[NumberOfSVTypes][2]={{65, 95},//For Legnth Standard Deviation Ratio Scores(100-ratio*100)
	{50,95},
	{85,99},
	{98,0}};
	double CCSASSBases[NumberOfSVTypes][2]={{0,0},
	{0,0},
	{4,0},
	{2,1}};
	double CCSASSCoverageMulti[NumberOfSVTypes][2]={{0.31,0.1},
	{0.33,0.13},
	{0.13,0.1},
	{0,0.18}};
	double CCSLSDRSs[NumberOfSVTypes][2]={{55, 95},
	{55,95},
	{85,99},
	{98,0}};
    double PreciseStandard=3;
    int MinimumPreciseTemplates=5;

	double TotalCoverage=0;

	std::string CustomFilterParas[NumberOfSVTypes]={"","","",""};
	std::string CustomClusterParas[NumberOfSVTypes]={"","","",""};

	bool CallByContig=false;
	bool WeightPosLength=false;
	bool WeightFilter=false;
	bool IndependantMerge=false;

	bool CigarMerge=0;//0: Omni.B, 1:Simple, 2:Simple Regional
	int OmniA[2]={800,200};
	int OmniB[2]={800,200};
	int OmniBMaxEdges=8;
	int OmniBCountLimit=20;
	double OmniBScoreBRatio=0.1;
	int MinPosSTD[NumberOfSVTypes]={-1,-1,-1,-1};//Filter out clusters that have position stds > MinPosSTD, -1: don't filter.
	bool CalcPosSTD=false;
	bool FID=true;

	double HPSameCompareRatio[NumberOfSVTypes]={0.8,0.8,0.4,1.0};
	double HPDiffCompareRatio[NumberOfSVTypes]={1.4,1.5,1.4,1.3};

	double HPRatio[NumberOfSVTypes]={0.8,0.7,0.2,0.5};
	double HomoRatio[NumberOfSVTypes]={0.8,0.75,0.9,0.65};
	double HomoMinus[NumberOfSVTypes]={0,0.05,0,0.7};
	double HomoMinusRatio[NumberOfSVTypes]={0.2,0.25,0.6,0.2};
	double NonHomoMinus[NumberOfSVTypes]={0.0,0.0,0.0,0.0};
	double NonHomoMinusRatio[NumberOfSVTypes]={-0.05,0.3,0,0.4};
	double LowHPPlus[NumberOfSVTypes]={0,0,0,0};
	double LowHPPlusRatio[NumberOfSVTypes]={0.25,0.25,0.9,0};
};

#endif
