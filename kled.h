#pragma once
#ifndef KLED_KLED_H
#define KLED_KLED_H

#include <vector>

const int NumberOfSVTypes=4;//Default is static so is fine.

struct Arguments {
	const char * Version="1.0";
	bool ShowVersion=false;
	bool ShowHelp=false;
	std::string CommandLine;
	std::string RunHash;//Storing std::hash<std::string>(Version+Arguments) in hex.
	int TestN=0;
	const char * ReferenceFileName=0;
	std::vector<const char *> BamFileNames;
	std::vector<const char *> CallingContigs;
	std::string SampleName="*";
	int ThreadN=8;
	bool NoFilter=false;
	int MinSVLen=30;
	int MinMappingQuality=20;
	int MinTemplateLength=500;
	int DelMinMaxMergeDis=500;//min maxmergedis, if CurrentLength*MaxMergeDisPortion>MinMaxMergeDis, MaxMergeDiss=CurrentLength*MaxMergeDisPortion
	double DelMaxMergePortion=0.2;
	int CoverageWindowSize=100;
	int InsClipTolerance=10;
	int InsMaxGapSize=50000;
	int MaxClusterSize=1000;//will limit (sampling) Cluster size to [MaxClusterSize, 2*MaxClusterSize)
	int ClusteringMaxMergeRange=20000;
	int ClusteringBatchSize=10000;
    float BrotherhoodTypeRatios[NumberOfSVTypes]={0.5,0.5,0.1,0.5};//For SV Types
    int BrotherhoodTypeForceBrothers[NumberOfSVTypes]={10,50,500,100};
    float BrotherhoodCCSTypeRatios[NumberOfSVTypes]={2.0,2.5,2.0,2.0};//For SV Types
    int BrotherhoodCCSTypeForceBrothers[NumberOfSVTypes]={100,200,100,100};
	bool AllCLR=false;
	bool AllCCS=false;
	double ASSBases[NumberOfSVTypes][2]=//Layers of base filter of the addition of supporting segmentations and templates
	{{3,1}//DEL
	,{1,1}//INS
	,{0,2}//DUP
	,{1,1}};//INV, only ST no SS
	double ASSCoverageMulti[NumberOfSVTypes][2]={{0.5,0.3},
	{0.4,0.3},
	{0.4,0.4},
	{0.5,0.5}};
	double LSDRSs[NumberOfSVTypes][2]={{80, 95},//For Legnth Standard Deviation Ratio Scores(100-ratio*100)
	{65,90},
	{95,90},
	{55,90}};
	double CLRASSBases[NumberOfSVTypes][2]=//Layers of base filter of the addition of supporting segmentations and templates
	{{3,1}//DEL
	,{2,4}//INS
	,{0,2}//DUP
	,{1,1}};//INV, only ST no SS
	double CLRASSCoverageMulti[NumberOfSVTypes][2]={{0.2,0.2},
	{0.4,0.1},
	{0.4,0.4},
	{0.5,0.5}};
	double CLRLSDRSs[NumberOfSVTypes][2]={{30, 95},//For Legnth Standard Deviation Ratio Scores(100-ratio*100)
	{55,85},
	{95,90},
	{55,90}};
	double CCSASSBases[NumberOfSVTypes][2]={{10,3},
	{5,1},
	{10,3},//Need further polish
	{1,1}};//Need further polish
	double CCSASSCoverageMulti[NumberOfSVTypes][2]={{0.5,0.1},
	{0.3,0.2},
	{0.5,0.1},
	{0.5,0.5}};
	double CCSLSDRSs[NumberOfSVTypes][2]={{0, 80},
	{55,80},
	{0,80},
	{55,90}};
    double PreciseStandard=3;
    int MinimumPreciseTemplates=5;

	std::string CustomFilterParas[NumberOfSVTypes]={"","","",""};
	std::string CustomClusterParas[NumberOfSVTypes]={"","","",""};
};

#endif