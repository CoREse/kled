#include <vector>
// #include "variant.h"
#include "signature.h"
#include "optutils/OptHelper.h"
#include "contig.h"
#include "input.h"
#include "kled.h"
#include "clustering.h"
#include "report.h"
#include "htslib/htslib/faidx.h"
#include <algorithm>
#include "crelib/crelib.h"
#include <omp.h>
#include <iterator>
#include <functional>
#include <sstream>
#include <iostream>
#ifdef DEBUG
#include <fstream>

// include headers that implement a archive in simple text format
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif
using namespace std;
using namespace cre;

void sortAndDeDup(vector<Signature> &V)
{
	if (V.size()==0) return;
	sort(V.begin(),V.begin()+V.size());
	for (int i=1;i<V.size();++i)
	{
		if (V[i]==V[i-1]) V[i].Type=-1;
	}
}

void showVersion(Arguments & Args)
{
	printf("Kled version %s.\n",Args.Version);
}

vector<string> split(string line, string delimiter=" ")
{
    vector<string> items;
    while (line!="")
    {
        int pos=line.find(delimiter);
        items.push_back(line.substr(0, pos));
        if (pos!=-1)
            line.erase(0, pos + delimiter.length());
        else
            line.erase(0,pos);
    }
    return items;
}

bool analyzeCustomParas(Arguments & Args)
{
	for (int i=0;i<NumberOfSVTypes;++i)
	{
		if (Args.CustomClusterParas[i]!="")
		{
			vector<string> sl=split(Args.CustomClusterParas[i], ",");
			if (sl.size()!=4 && sl.size()!=7 && sl.size()!=5)
			{
				fprintf(stderr,"Error: wrong parameter format %s.\n", Args.CustomClusterParas[i].c_str());
				return false;
			}
			if (sl[0][0]!='*')
			Args.BrotherhoodTypeForceBrothers[i]=atoi(sl[0].c_str());
			if (sl[1][0]!='*')
			Args.BrotherhoodTypeRatios[i]=atof(sl[1].c_str());
			if (sl[2][0]!='*')
			Args.BrotherhoodTypeLengthMinEndurance[i]=atoi(sl[2].c_str());
			if (sl[3][0]!='*')
			Args.BrotherhoodTypeLengthRatios[i]=atof(sl[3].c_str());
			if (sl.size()==4 || sl.size()==5)
			{
				Args.BrotherhoodNearRanges[i]=-1;
			}
			else
			{
				if (sl[4][0]!='*')
				Args.BrotherhoodNearRanges[i]=atoi(sl[4].c_str());
				if (sl[5][0]!='*')
				Args.BrotherhoodTypeForceBrothers2[i]=atoi(sl[5].c_str());
				if (sl[6][0]!='*')
				Args.BrotherhoodTypeLengthRatios2[i]=atof(sl[6].c_str());
			}
		}
		if (Args.CustomFilterParas[i]!="")
		{
			vector<string> sl=split(Args.CustomFilterParas[i], ",");
			if (sl.size()!=6)
			{
				fprintf(stderr,"Error: wrong parameter format %s.\n",Args.CustomFilterParas[i].c_str());
				return false;
			}
			if (sl[0][0]!='*')
			Args.ASSBases[i][0]=atof(sl[0].c_str());
			if (sl[3][0]!='*')
			Args.ASSBases[i][1]=atof(sl[3].c_str());
			if (sl[1][0]!='*')
			Args.ASSCoverageMulti[i][0]=atof(sl[1].c_str());
			if (sl[4][0]!='*')
			Args.ASSCoverageMulti[i][1]=atof(sl[4].c_str());
			if (sl[2][0]!='*')
			Args.LSDRSs[i][0]=atof(sl[2].c_str());
			if (sl[5][0]!='*')
			Args.LSDRSs[i][1]=atof(sl[5].c_str());
		}
	}
	return true;
}

bool toCall(const Contig & C, const Arguments &Args)
{
	if (Args.CallingContigs.size()!=0)
	{
		bool ToCall=false;
		for (int k=0;k<Args.CallingContigs.size();++k)
		{
			if (C.Name==Args.CallingContigs[k])
			{
				ToCall=true;
				break;
			}
		}
		return ToCall;
	}
	return true;
}

#pragma omp declare reduction(RecordVectorConc: vector<VCFRecord>: omp_out.insert(omp_out.end(),make_move_iterator(omp_in.begin()),make_move_iterator(omp_in.end())))
//#pragma omp declare reduction(RecordListConc: list<VCFRecord>: omp_out.splice(omp_out.end(),omp_in))

Arguments Args;
int main(int argc, const char* argv[])
{
	string RunString=Args.Version;
	for (int i=1;i<argc;++i) RunString+=string(" ")+argv[i];
	Args.CommandLine=argv[0];
	for (int i=1;i<argc;++i) Args.CommandLine+=string(" ")+argv[i];
	size_t Hash=hash<string>()(RunString);
	stringstream ss;
	ss<<std::hex<<Hash;
	ss>>Args.RunHash;
	bool NoHeader=false;
	OptHelper OH=OptHelper("kled [Options] Bam1 [Bam2] [Bam3] ...");
    // OH.addOpt('N', 0, 1, "TestNumber", "for test notation",'i',&(Args.TestN));
    OH.addOpt('R', "Ref", 1, "FileName", "Indicate Reference Fasta File(required)",'s',&(Args.ReferenceFileName));
    OH.addOpt('C', 0, 1, "ContigName", "Only call variants in Contig(s), can occur multiple times",'s',&(Args.CallingContigs),true);
    OH.addOpt('S', 0, 1, "SampleName", "Sample name, if not given, kled will try to get it from the first bam file",'S',&(Args.SampleName));
    OH.addOpt('t', "threads", 1, "Number", "Number of threads.",'i',&(Args.ThreadN));
    OH.addOpt('h', "help", 0, "", "Show this help and exit.",'b',&(Args.ShowHelp));
    OH.addOpt('v', "version", 0, "", "Show version and exit.",'b',&(Args.ShowVersion));
    OH.addOpt(0, "BC", 0, "", "Calling contig by contig, cost less memory.",'b',&(Args.CallByContig));
    OH.addOpt(0, "CCS", 0, "", "Use default parameters for CCS data (will overwrite cluster and filter parameters).",'b',&(Args.AllCCS));
    OH.addOpt(0, "CLR", 0, "", "Use default parameters for CLR data (will overwrite cluster and filter parameters).",'b',&(Args.AllCLR));
    OH.addOpt(0, "DelClusterParas", 1, "Fixed,Ratio,MinLengthEndurance,LengthRatio[,NearRange,LengthDiff,LengthRatio2]", "Custom clustering parameters for deletions, if later 3 not given or NearRange=-1 use single layer clustering. Value * for defaults.",'S',&(Args.CustomClusterParas[0]));
		OH.addOpt(0, "DelClusterFixed", 1, "FixedDistance", "Fixed distance clustering parameter for deletions.",'i',&(Args.BrotherhoodTypeForceBrothers[0]));
		OH.addOpt(0, "DelClusterRatio", 1, "DistanceRatio", "Distance ratio clustering parameter for deletions.",'F',&(Args.BrotherhoodTypeRatios[0]));
		OH.addOpt(0, "DelClusterMinLengthEdurance", 1, "MinLengthEndurance", "Min length endurance clustering parameter for deletions.",'i',&(Args.BrotherhoodTypeLengthMinEndurance[0]));
		OH.addOpt(0, "DelClusterLengthRatio", 1, "LenthRatio", "Length ratio clustering parameter for deletions.",'F',&(Args.BrotherhoodTypeLengthRatios[0]));
		OH.addOpt(0, "DelClusterNearRange", 1, "NearRange", "Near range clustering parameter for deletions.",'i',&(Args.BrotherhoodNearRanges[0]));
		OH.addOpt(0, "DelClusterFixed2", 1, "FixedDistance", "Fixed distance 2 clustering parameter for deletions.",'i',&(Args.BrotherhoodTypeForceBrothers2[0]));
		OH.addOpt(0, "DelClusterLengthRatio2", 1, "LenthRatio", "Length ratio 2 clustering parameter for deletions.",'F',&(Args.BrotherhoodTypeLengthRatios2[0]));
    OH.addOpt(0, "InsClusterParas", 1, "Fixed,Ratio,MinLengthEndurance,LengthRatio[,NearRange,LengthDiff,LengthRatio2]", "Custom clustering parameters for insertions, if later 3 not given or NearRange=-1 use single layer clustering. Value * for defaults.",'S',&(Args.CustomClusterParas[1]));
    	OH.addOpt(0, "InsClusterFixed", 1, "FixedDistance", "Fixed distance clustering parameter for insertions.",'i',&(Args.BrotherhoodTypeForceBrothers[1]));
		OH.addOpt(0, "InsClusterRatio", 1, "DistanceRatio", "Distance ratio clustering parameter for insertions.",'F',&(Args.BrotherhoodTypeRatios[1]));
		OH.addOpt(0, "InsClusterMinLengthEdurance", 1, "MinLengthEndurance", "Min length endurance clustering parameter for insertions.",'i',&(Args.BrotherhoodTypeLengthMinEndurance[1]));
		OH.addOpt(0, "InsClusterLengthRatio", 1, "LenthRatio", "Length ratio clustering parameter for insertions.",'F',&(Args.BrotherhoodTypeLengthRatios[1]));
		OH.addOpt(0, "InsClusterNearRange", 1, "NearRange", "Near range clustering parameter for insertions.",'i',&(Args.BrotherhoodNearRanges[1]));
		OH.addOpt(0, "InsClusterFixed2", 1, "FixedDistance", "Fixed distance 2 clustering parameter for insertions.",'i',&(Args.BrotherhoodTypeForceBrothers2[1]));
		OH.addOpt(0, "InsClusterLengthRatio2", 1, "LenthRatio", "Length ratio 2 clustering parameter for insertions.",'F',&(Args.BrotherhoodTypeLengthRatios2[1]));
	OH.addOpt(0, "DupClusterParas", 1, "Fixed,Ratio,MinLengthEndurance,LengthRatio[,NearRange,LengthDiff,LengthRatio2]", "Custom clustering parameters for duplications, if later 3 not given or NearRange=-1 use single layer clustering. Value * for defaults.",'S',&(Args.CustomClusterParas[2]));
    	OH.addOpt(0, "DupClusterFixed", 1, "FixedDistance", "Fixed distance clustering parameter for duplications.",'i',&(Args.BrotherhoodTypeForceBrothers[2]));
		OH.addOpt(0, "DupClusterRatio", 1, "DistanceRatio", "Distance ratio clustering parameter for duplications.",'F',&(Args.BrotherhoodTypeRatios[2]));
		OH.addOpt(0, "DupClusterMinLengthEdurance", 1, "MinLengthEndurance", "Min length endurance clustering parameter for duplications.",'i',&(Args.BrotherhoodTypeLengthMinEndurance[2]));
		OH.addOpt(0, "DupClusterLengthRatio", 1, "LenthRatio", "Length ratio clustering parameter for duplications.",'F',&(Args.BrotherhoodTypeLengthRatios[2]));
		OH.addOpt(0, "DupClusterNearRange", 1, "NearRange", "Near range clustering parameter for duplications.",'i',&(Args.BrotherhoodNearRanges[2]));
		OH.addOpt(0, "DupClusterFixed2", 1, "FixedDistance", "Fixed distance 2 clustering parameter for duplications.",'i',&(Args.BrotherhoodTypeForceBrothers2[2]));
		OH.addOpt(0, "DupClusterLengthRatio2", 1, "LenthRatio", "Length ratio 2 clustering parameter for duplications.",'F',&(Args.BrotherhoodTypeLengthRatios2[2]));
	OH.addOpt(0, "InvClusterParas", 1, "Fixed,Ratio,MinLengthEndurance,LengthRatio[,NearRange,LengthDiff,LengthRatio2]", "Custom clustering parameters for inversions, if later 3 not given or NearRange=-1 use single layer clustering. Value * for defaults.",'S',&(Args.CustomClusterParas[3]));
    	OH.addOpt(0, "InvClusterFixed", 1, "FixedDistance", "Fixed distance clustering parameter for inversions.",'i',&(Args.BrotherhoodTypeForceBrothers[3]));
		OH.addOpt(0, "InvClusterRatio", 1, "DistanceRatio", "Distance ratio clustering parameter for inversions.",'F',&(Args.BrotherhoodTypeRatios[3]));
		OH.addOpt(0, "InvClusterMinLengthEdurance", 1, "MinLengthEndurance", "Min length endurance clustering parameter for inversions.",'i',&(Args.BrotherhoodTypeLengthMinEndurance[3]));
		OH.addOpt(0, "InvClusterLengthRatio", 1, "LenthRatio", "Length ratio clustering parameter for inversions.",'F',&(Args.BrotherhoodTypeLengthRatios[3]));
		OH.addOpt(0, "InvClusterNearRange", 1, "NearRange", "Near range clustering parameter for inversions.",'i',&(Args.BrotherhoodNearRanges[3]));
		OH.addOpt(0, "InvClusterFixed2", 1, "FixedDistance", "Fixed distance 2 clustering parameter for inversions.",'i',&(Args.BrotherhoodTypeForceBrothers2[3]));
		OH.addOpt(0, "InvClusterLengthRatio2", 1, "LenthRatio", "Length ratio 2 clustering parameter for inversions.",'F',&(Args.BrotherhoodTypeLengthRatios2[3]));
	OH.addOpt(0, "DelFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for deletions. Value * for defaults.",'S',&(Args.CustomFilterParas[0]));
    OH.addOpt(0, "InsFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for insertions. Value * for defaults.",'S',&(Args.CustomFilterParas[1]));
    OH.addOpt(0, "DupFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for duplications. Value * for defaults.",'S',&(Args.CustomFilterParas[2]));
    OH.addOpt(0, "InvFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for inversions. Value * for defaults.",'S',&(Args.CustomFilterParas[3]));
    OH.addOpt(0, "NOF", 0, "", "No filter, output all results.",'b',&(Args.NoFilter));
    OH.addOpt(0, "F2", 0, "", "Output all results with ST>=2.",'b',&(Args.Filter2ST));
    OH.addOpt('m', 0, 1, "SVLEN", "Minimum SV length.",'i',&(Args.MinSVLen));
    OH.addOpt('q', 0, 1, "Quality", "Minimum mapping quality.",'i',&(Args.MinMappingQuality));
    OH.addOpt('l', 0, 1, "Length", "Minimum template length.",'i',&(Args.MinTemplateLength));
    OH.addOpt('d', 0, 1, "Distance", "Minimum max merge distance of signature merging during CIGAR signature collection for simple merge.",'i',&(Args.DelMinMaxMergeDis));
    OH.addOpt('D', 0, 1, "Distance", "Maximum max merge distance of signature merging during CIGAR signature collection for simple merge.",'i',&(Args.DelMaxMaxMergeDis));
    OH.addOpt('p', 0, 1, "Portion", "Max merge portion of signature merging during CIGAR signature collection for simple merge.",'F',&(Args.DelMaxMergePortion));
    OH.addOpt('c', 0, 1, "Size", "Coverage window size.",'i',&(Args.CoverageWindowSize));
    OH.addOpt('M', 0, 1, "Size", "Max cluster size, will resize to this value if a cluster is larger than this.",'i',&(Args.MaxClusterSize));
    OH.addOpt(0, "InsClipTolerance", 1, "Size", "Insertion clip signature distance tolerance.",'i',&(Args.InsClipTolerance));
    OH.addOpt(0, "InsMaxGapSize", 1, "Size", "Insertion clip signature max gap size.",'i',&(Args.InsMaxGapSize));
    OH.addOpt(0, "ClusteringBatchSize", 1, "Size", "Batch size of multihreading when clustering.",'i',&(Args.ClusteringBatchSize));
    OH.addOpt(0, "CigarMerge", 1, "MergeType", "CigarMergeType, 0 for Omni.B, 1 for simple, 2 for simple regional.",'i',&(Args.CigarMerge));
    OH.addOpt(0, "OMaxE", 1, "Size", "Max edges(depth) for omni.b merge. This will grow complexity exponentially.",'i',&(Args.OmniBMaxEdges));
    OH.addOpt(0, "OCountLimit", 1, "Size", "Relevant alignments count limit for omni.b merge.",'i',&(Args.OmniBCountLimit));
    OH.addOpt(0, "OScoreBRatio", 1, "Ratio", "ScoreB ratio for omni.b merge.",'F',&(Args.OmniBScoreBRatio));
    OH.addOpt(0, "OmniA", 1, "ValueA", "A for omni.b merge.",'i',&(Args.OmniA));
    OH.addOpt(0, "OmniB", 1, "ValueB", "B for omni.b merge.",'i',&(Args.OmniB));
	#ifdef DEBUG
    OH.addOpt(0, "NOH", 0, "", "No header, for test",'b',&(NoHeader));
	string WriteSigDataFileName="";
	string ReadSigDataFileName="";
    OH.addOpt(0, "WS", 1, "Data file name", "File name to write signature data",'S',&(WriteSigDataFileName));
    OH.addOpt(0, "RS", 1, "Data file name", "File name to read signature data",'S',&(ReadSigDataFileName));
	#endif
    OH.getOpts(argc,argv);

	if (Args.ShowHelp)
	{
		OH.showhelp();
		exit(0);
	}

	if (Args.ShowVersion)
	{
		showVersion(Args);
		exit(0);
	}

	Args.BamFileNames=OH.Args;

	if (Args.ReferenceFileName==0 || Args.BamFileNames.size()==0)
	{
		OH.showhelp();
		exit(1);
	}
	
	if (!analyzeCustomParas(Args))
	{
		OH.showhelp();
		exit(1);
	}

	if (Args.AllCLR)
	{
		for (int i=0;i<NumberOfSVTypes;++i)
		{
			Args.BrotherhoodTypeRatios[i]=				Args.BrotherhoodCLRTypeRatios[i];
			Args.BrotherhoodTypeForceBrothers[i]=		Args.BrotherhoodCLRTypeForceBrothers[i];
			Args.BrotherhoodTypeLengthRatios[i]=		Args.BrotherhoodCLRTypeLengthRatios[i];
			Args.BrotherhoodTypeLengthMinEndurance[i]=	Args.BrotherhoodCLRTypeLengthMinEndurance[i];
			Args.BrotherhoodNearRanges[i]=				Args.BrotherhoodCLRNearRanges[i];
			Args.BrotherhoodTypeForceBrothers2[i]=		Args.BrotherhoodCLRTypeForceBrothers2[i];
			Args.BrotherhoodTypeLengthRatios2[i]=		Args.BrotherhoodCLRTypeLengthRatios2[i];
		}
	}
	else if (Args.AllCCS)
	{
		for (int i=0;i<NumberOfSVTypes;++i)
		{
			Args.BrotherhoodTypeRatios[i]=				Args.BrotherhoodCCSTypeRatios[i];
			Args.BrotherhoodTypeForceBrothers[i]=		Args.BrotherhoodCCSTypeForceBrothers[i];
			Args.BrotherhoodTypeLengthRatios[i]=		Args.BrotherhoodCCSTypeLengthRatios[i];
			Args.BrotherhoodTypeLengthMinEndurance[i]=	Args.BrotherhoodCCSTypeLengthMinEndurance[i];
			Args.BrotherhoodNearRanges[i]=				Args.BrotherhoodCCSNearRanges[i];
			Args.BrotherhoodTypeForceBrothers2[i]=		Args.BrotherhoodCCSTypeForceBrothers2[i];
			Args.BrotherhoodTypeLengthRatios2[i]=		Args.BrotherhoodCCSTypeLengthRatios2[i];
		}
	}

	omp_set_num_threads(Args.ThreadN);
	// ThreadPool ThePool(8);

	updateTime("Starting kled, reading reference...");
	int NSeq;
	Contig * Contigs=getContigs(Args,NSeq);//,RDWindowSize);
    
	updateTime("Reading reference","Getting stats...");
	vector<int> AllTechs=getAllTechs(Args);
	// fprintf(stderr,"All techs:");for (int i=0;i<AllTechs.size();++i) fprintf(stderr," %d",AllTechs[i]);fprintf(stderr,"\n");

	vector<Stats> AllStats=getAllStats(Args.ReferenceFileName,Args.BamFileNames,AllTechs);

	for (int i=0;i<AllStats.size();++i) fprintf(stderr,"%f %f %f %f %f\n",AllStats[i].BelowIS,AllStats[i].MedianIS,AllStats[i].UpIS,AllStats[i].Mean,AllStats[i].SD);
	//exit(0);

	VCFHeader Header(Args.ReferenceFileName);
	addKledEntries(Header);
	for (int i=0;i<NSeq;++i)
	{
		if (toCall(Contigs[i],Args)) Header.addContig(Contigs[i]);
	}

	faidx_t * Ref=fai_load(Args.ReferenceFileName);
	// vector<vector<Variant>> VariantsByContig;
	// bool FirstBam=true;
	vector<Sam> SamFiles=initSam(Args);
	double TotalCoverage=0;//Accumulative
	unsigned ProcessedLength=0;
	// FILE * WindowsFile=fopen("/home/cre/workspace/kled/data/wins.txt","wb");
	updateTime("Getting stats", "Starting calling...");
	unsigned int SVCounts[NumberOfSVTypes];for (int i=0;i<NumberOfSVTypes;++i) SVCounts[i]=0;

	vector<vector<vector<Signature>>> TypeSignatures;//Contig-Type-Signatures
	vector<SegmentSet> ContigsAllPrimarySegments;
	vector<double *> CoverageWindowsPs;
	for (int i=0;i<NSeq;++i)
	{
		TypeSignatures.push_back(vector<vector<Signature>>());
		for (int j=0;j<NumberOfSVTypes;++j) TypeSignatures[i].push_back(vector<Signature>());
	}
	#ifdef DEBUG
	ifstream ifs;
	ofstream ofs;
	#endif
	if (!Args.CallByContig)
	{
		for (int i=0;i<NSeq;++i)
		{
			if (! toCall(Contigs[i],Args)) continue;
			updateTime("","Calling...");
			SegmentSet AllPrimarySegments;
			// vector<Signature> ContigTypeSignatures[NumberOfSVTypes];//For supported SV type
			unsigned int CoverageWindowSize=Args.CoverageWindowSize;
			unsigned int NumberOfCoverageWindows=Contigs[i].Size/CoverageWindowSize+1;
			double *CoverageWindows=new double[NumberOfCoverageWindows];
			#ifdef DEBUG
			if (ReadSigDataFileName!="")
			{
				if (!ifs.is_open())
					ifs.open(ReadSigDataFileName.c_str(),ios::binary);
				boost::archive::binary_iarchive ia(ifs);
				ia >> NumberOfCoverageWindows;
				for (int j=0;j<NumberOfCoverageWindows;++j) ia>>(CoverageWindows[j]);
				ia >> AllPrimarySegments;
				ia >> (TypeSignatures[i]);
			}
			else
			{
			#endif
			for (int k=0;k<Contigs[i].Size/CoverageWindowSize+1;++k) CoverageWindows[k]=0;
			collectSignatures(Contigs[i],TypeSignatures,AllPrimarySegments,Args,SamFiles,AllStats,AllTechs,CoverageWindows,0);
			AllPrimarySegments.sortNStat();
			#ifdef DEBUG
			}
			#endif
			ContigsAllPrimarySegments.push_back(AllPrimarySegments);
			CoverageWindowsPs.push_back(CoverageWindows);
		}
	}
	// int Skipped=0;
	vector<vector<VCFRecord>> ContigOutputs;
	vector<double> ContigWholeCoverage;
	vector<int> ContigIndex;
	int Called=0;
	for (int i=0;i<NSeq;++i)
	{
		if (toCall(Contigs[i],Args))
		{
			ContigIndex.push_back(Called++);
		}
		else ContigIndex.push_back(-1);
		ContigOutputs.push_back(vector<VCFRecord>());
		ContigWholeCoverage.push_back(0.0);
	}
	#pragma omp parallel for
	for (int i=0;i<NSeq;++i)
	{
		if (! toCall(Contigs[i],Args))
		{
			// ++Skipped;
			continue;
		}
		unsigned int CoverageWindowSize=Args.CoverageWindowSize;
		unsigned int NumberOfCoverageWindows=Contigs[i].Size/CoverageWindowSize+1;
		if (Args.CallByContig)
		{
			updateTime("","Calling...");
			SegmentSet AllPrimarySegments;
			// vector<Signature> ContigTypeSignatures[NumberOfSVTypes];//For supported SV type
			double *CoverageWindows=new double[NumberOfCoverageWindows];
			for (int k=0;k<Contigs[i].Size/CoverageWindowSize+1;++k) CoverageWindows[k]=0;
			collectSignatures(Contigs[i],TypeSignatures,AllPrimarySegments,Args,SamFiles,AllStats,AllTechs,CoverageWindows,0);
			AllPrimarySegments.sortNStat();
			ContigsAllPrimarySegments.push_back(AllPrimarySegments);
			CoverageWindowsPs.push_back(CoverageWindows);
		}
		double *CoverageWindows=CoverageWindowsPs[ContigIndex[i]];
		SegmentSet &AllPrimarySegments=ContigsAllPrimarySegments[ContigIndex[i]];
		vector<vector<Signature>> &ContigTypeSignatures=TypeSignatures[i];
		#ifdef DEBUG
		if (WriteSigDataFileName!="")
		{
			if (!ofs.is_open())
				ofs.open(WriteSigDataFileName.c_str(),ios::binary);
			boost::archive::binary_oarchive oa(ofs);
			oa << NumberOfCoverageWindows;
			for (int j=0;j<NumberOfCoverageWindows;++j) oa<<(CoverageWindows[j]);
			oa << AllPrimarySegments;
			oa << ContigTypeSignatures;
			continue;
		}
		else if (ReadSigDataFileName=="")
		{
			if (!ofs.is_open())
				ofs.open("data/SigData.dat",ios::binary);
			boost::archive::binary_oarchive oa(ofs);
			oa << NumberOfCoverageWindows;
			for (int j=0;j<NumberOfCoverageWindows;++j) oa<<(CoverageWindows[j]);
			oa << AllPrimarySegments;
			oa << ContigTypeSignatures;
		}
		#endif
		// fprintf(stderr,"%u\n",Contigs[i].Size-1);
		double *CoverageWindowsSums=NULL;//=(double*) malloc(sizeof(double)*(int)(NumberOfCoverageWindows+1));
		CoverageWindows[0]=0;
		// CoverageWindowsSums[0]=0;
		int CheckPointInterval=10000;
		double *CheckPoints=NULL;//=(double *)malloc(sizeof(double)*(int)(NumberOfCoverageWindows/CheckPointInterval+1));
		// CheckPoints[0]=0;
		// for (int i=1;i<NumberOfCoverageWindows+1;++i)
		// {
		// 	CoverageWindowsSums[i]=CoverageWindowsSums[i-1]+CoverageWindows[i];
		// 	if (i%CheckPointInterval==0)
		// 	{
		// 		CheckPoints[(int)i/CheckPointInterval]=CoverageWindowsSums[i];
		// 		CoverageWindowsSums[i]=0;
		// 	}
		// }
		double &WholeCoverage=ContigWholeCoverage[i];
		WholeCoverage=getAverageCoverage(0,Contigs[i].Size-1,CoverageWindows,Args, CoverageWindowsSums, CheckPoints, CheckPointInterval);
		// int NameLength=Contigs[i].Name.length();
		// fwrite(&(NameLength),sizeof(int),1,WindowsFile);
		// fwrite(Contigs[i].Name.c_str(),1,Contigs[i].Name.length(),WindowsFile);
		// fwrite(&(Contigs[i].Size),sizeof(unsigned),1,WindowsFile);
		// fwrite(&(NumberOfCoverageWindows),sizeof(unsigned),1,WindowsFile);
		// fwrite(CoverageWindows,sizeof(double),NumberOfCoverageWindows,WindowsFile);
		TotalCoverage=TotalCoverage*((double)(ProcessedLength)/(double)(ProcessedLength+Contigs[i].Size));
		TotalCoverage+=WholeCoverage*((double)(Contigs[i].Size)/(double)(ProcessedLength+Contigs[i].Size));
		Args.TotalCoverage=TotalCoverage;
		ProcessedLength+=Contigs[i].Size;
		// double WholeCoverage=CoverageWindowsSums[(int)(Contigs[i].Size/CoverageWindowSize+1)]/(Contigs[i].Size/CoverageWindowSize+1);
		// continue;
		int totalsig=0,cigardel=0, cigarins=0, cigardup=0, drpdel=0, drpdup=0, clipdel=0, clipins=0, clipdup=0, clipinv=0, NGS=0, SMRT=0, NGSCigar=0, NGSClip=0;
		for (int m=0;m<NumberOfSVTypes;++m)
		{
			vector<Signature>& ContigSignatures=ContigTypeSignatures[m];
			totalsig+=ContigSignatures.size();
			for (int j=0;j<ContigSignatures.size();++j)
			{
				if (ContigSignatures[j].Type==0)
				{
					if (ContigSignatures[j].SupportedSV==0) ++cigardel;
					if (ContigSignatures[j].SupportedSV==1) ++cigarins;
					if (ContigSignatures[j].SupportedSV==2) ++cigardup;
					if (ContigSignatures[j].Tech==1) ++NGSCigar;
				}
				else if (ContigSignatures[j].Type==1)
				{
					if (ContigSignatures[j].SupportedSV==0) ++drpdel;
					if (ContigSignatures[j].SupportedSV==2) ++drpdup;
				}
				else
				{
					if (ContigSignatures[j].SupportedSV==0) ++clipdel;
					if (ContigSignatures[j].SupportedSV==1) ++clipins;
					if (ContigSignatures[j].SupportedSV==2) ++clipdup;
					if (ContigSignatures[j].SupportedSV==3) ++clipinv;
					if (ContigSignatures[j].Tech==1) ++NGSClip;
				}
				if (ContigSignatures[j].Tech==0) ++SMRT;
				else ++NGS;
			}
		}
		fprintf(stderr,"%s: %d\n, cigardel: %d, cigarins: %d, cigardup: %d, drpdel: %d, drpdup: %d, clipdel: %d, clipins: %d, clipdup: %d, clipinv: %d, NGS: %d(Cigar: %d, Clip: %d), SMRT: %d. Contig Size:%ld, Average Coverage: %lf, Total Average Coverage: %lf\n",Contigs[i].Name.c_str(),totalsig,cigardel, cigarins, cigardup, drpdel, drpdup, clipdel, clipins, clipdup, clipinv, NGS, NGSCigar, NGSClip, SMRT, Contigs[i].Size, WholeCoverage, TotalCoverage);
		// extern int forwarded, reversed, nread1, nread2;
		// fprintf(stderr,"%forwarded:%d,reversed:%d, nread1:%d, nread2:%d",forwarded,reversed,nread1,nread2);
		updateTime("Getting signatures","Clustering...");
		vector<vector<Signature>> SignatureTypeClusters[NumberOfSVTypes];
		vector<ClusterCore> SignatureTypeClusterCores[NumberOfSVTypes];
		for (int k=0;k<NumberOfSVTypes;++k)
		{
			sortAndDeDup(ContigTypeSignatures[k]);
			for (unsigned d=0;d<ContigTypeSignatures[k].size();++d) ContigTypeSignatures[k][d].setID(d);
			clustering(ContigTypeSignatures[k],SignatureTypeClusters[k],SignatureTypeClusterCores[k],AllStats[i],Args);
		}
		updateTime("Clustering","Generating results...");
		vector<vector<Signature>> SignatureClusters;
		vector<ClusterCore> SignatureClusterCores;
		for (int k=0;k<NumberOfSVTypes;++k)
		{
			SignatureClusters.insert(SignatureClusters.end(),make_move_iterator(SignatureTypeClusters[k].begin()),make_move_iterator(SignatureTypeClusters[k].end()));
			SignatureClusterCores.insert(SignatureClusterCores.end(),make_move_iterator(SignatureTypeClusterCores[k].begin()),make_move_iterator(SignatureTypeClusterCores[k].end()));
		}
		
		totalsig=0,cigardel=0, cigarins=0, cigardup=0, drpdel=0, drpdup=0, clipdel=0, clipins=0, clipdup=0, clipinv=0, NGS=0, SMRT=0;
		for (int m=0;m<SignatureClusters.size();++m)
		{
			// vector<Signature>& ContigSignatures=ContigTypeSignatures[m];
			// totalsig+=ContigSignatures.size();
			for (int j=0;j<SignatureClusters[m].size();++j)
			{
				if (SignatureClusters[m][j].Type==0)
				{
					if (SignatureClusters[m][j].SupportedSV==0) ++cigardel;
					if (SignatureClusters[m][j].SupportedSV==1) ++cigarins;
					if (SignatureClusters[m][j].SupportedSV==2) ++cigardup;
				}
				else if (SignatureClusters[m][j].Type==1)
				{
					if (SignatureClusters[m][j].SupportedSV==0) ++drpdel;
					if (SignatureClusters[m][j].SupportedSV==2) ++drpdup;
				}
				else
				{
					if (SignatureClusters[m][j].SupportedSV==0) ++clipdel;
					if (SignatureClusters[m][j].SupportedSV==1) ++clipins;
					if (SignatureClusters[m][j].SupportedSV==2) ++clipdup;
					if (SignatureClusters[m][j].SupportedSV==3) ++clipinv;
				}
				if (SignatureClusters[m][j].Tech==0) ++SMRT;
				else ++NGS;
			}
		}
		fprintf(stderr,"%s: %d\n, cigardel: %d, cigarins: %d, cigardup: %d, drpdel: %d, drpdup: %d, clipdel: %d, clipins: %d, clipdup: %d, clipinv: %d, NGS: %d, SMRT: %d. Contig Size:%ld, Average Coverage: %lf, Total Average Coverage: %lf\n",Contigs[i].Name.c_str(),totalsig,cigardel, cigarins, cigardup, drpdel, drpdup, clipdel, clipins, clipdup, clipinv, NGS, SMRT, Contigs[i].Size, WholeCoverage, TotalCoverage);

		vector<VCFRecord> Records;
		// #pragma omp parallel for reduction(RecordVectorConc:Records)
		for (int j=0;j<SignatureClusters.size();++j)
		{
			// ++Times[omp_get_thread_num()];
			if (SignatureClusters[j].size()==0) continue;
			ClusterCore Core;
			if (SignatureClusterCores.size()!=0) Core=SignatureClusterCores[j];
			Records.push_back(VCFRecord(Contigs[i],Ref,SignatureClusters[j], Core, AllPrimarySegments,CoverageWindows, TotalCoverage, Args, CoverageWindowsSums, CheckPoints, CheckPointInterval));
		}
		updateTime("Results generation","Sorting results...");
		sort(Records.data(),Records.data()+Records.size());
		// Records.sort();
		updateTime("Results sorting","");

		for (auto r: Records)
		{
			if (!r.Keep) continue;
			ContigOutputs[i].push_back(r);
			// r.genotype(Contigs[i],AllPrimarySegments,CoverageWindows,CoverageWindowsSums,CheckPoints,CheckPointInterval,Args);
			// r.resolveRef(Contigs[i],Ref,SVCounts[r.getSVTypeI()], WholeCoverage,Args);
			// ++SVCounts[r.getSVTypeI()];
			// // printf("\n%s",string(r).c_str());
		}
		delete CoverageWindows;
		CoverageWindowsPs[i]=nullptr;
		TypeSignatures[i].clear();
		ContigsAllPrimarySegments[i]=SegmentSet();
		// free(CoverageWindowsSums);
		// free(CheckPoints);
	}
	#ifdef DEBUG
	if (ifs.is_open()) ifs.close();
	if (ofs.is_open()) ofs.close();
	#endif
	updateTime("All contigs finished.","Outputing results...");
	if (!NoHeader)
	{
		Header.addSample(Args.SampleName.c_str());
		printf(Header.genHeader(Args).c_str());
	}
	for (int i=0;i<NSeq;++i) for (int j=0;j<ContigOutputs[i].size();++j)
	{
		VCFRecord &r=ContigOutputs[i][j];
		if (!r.Keep) continue;
		// r.genotype(Contigs[i],AllPrimarySegments,CoverageWindows,CoverageWindowsSums,CheckPoints,CheckPointInterval,Args);
		r.resolveRef(Contigs[i],Ref,SVCounts[r.getSVTypeI()], ContigWholeCoverage[i],Args);
		++SVCounts[r.getSVTypeI()];
		printf("\n%s",string(r).c_str());
	}

	//report(VariantsByContig);
	fai_destroy(Ref);
	free(Contigs);
	closeSam(SamFiles);
	// fclose(WindowsFile);
    return 0;
}