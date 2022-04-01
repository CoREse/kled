#pragma once
#ifndef KLED_KLED_H
#define KLED_KLED_H

#include <vector>

const int NumberOfSVTypes=3;//Default is static so is fine.

struct Arguments {
	int TestN=0;
	const char * ReferenceFileName=0;
	std::vector<const char *> BamFileNames;
	std::vector<const char *> CallingContigs;
	std::string SampleName="*";
	int ThreadN=8;
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
	bool AllCCS=false;
	double ASSBases[NumberOfSVTypes][2]=//Layers of base filter of the addition of supporting segmentations and templates
	{{10,3}//DEL
	,{10,1}//INS
	,{20,2}};//DUP
	double ASSCoverageMulti[NumberOfSVTypes][2]={{0.5,0.3},
	{0.5,0.3},
	{0.5,0.3}};
	double LSDRSs[NumberOfSVTypes][2]={{0, 80},//For Legnth Standard Deviation Ratio Scores(100-ratio*100)
	{55,80},
	{60,90}};
	// double InsASSBases[2]={10,1};
	// double InsASSCoverageMulti[2]={0.5,0.3};
	// double InsLSDRSs[2]={55, 80};
	double CCSASSBases[NumberOfSVTypes][2]={{10,3},
	{5,1},
	{10,3}};//Need further polish
	double CCSASSCoverageMulti[NumberOfSVTypes][2]={{0.5,0.1},
	{0.3,0.2},
	{0.5,0.1}};
	double CCSLSDRSs[NumberOfSVTypes][2]={{0, 80},
	{55,80},
	{0,80}};
	// double CCSInsASSBases[2]={5,1};
	// double CCSInsASSCoverageMulti[2]={0.3,0.2};
	// double CCSInsLSDRSs[2]={55, 80};
    double PreciseStandard=3;
    int MinimumPreciseTemplates=5;
};

#endif