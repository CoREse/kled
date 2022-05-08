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
			if (sl.size()!=2)
			{
				fprintf(stderr,"Error: wrong parameter format %s.\n", Args.CustomClusterParas[i].c_str());
				return false;
			}
			Args.BrotherhoodTypeForceBrothers[i]=atoi(sl[0].c_str());
			Args.BrotherhoodTypeRatios[i]=atof(sl[1].c_str());
		}
		if (Args.CustomFilterParas[i]!="")
		{
			vector<string> sl=split(Args.CustomFilterParas[i], ",");
			if (sl.size()!=6)
			{
				fprintf(stderr,"Error: wrong parameter format %s.\n",Args.CustomFilterParas[i].c_str());
				return false;
			}
			Args.ASSBases[i][0]=atof(sl[0].c_str());
			Args.ASSBases[i][1]=atof(sl[3].c_str());
			Args.ASSCoverageMulti[i][0]=atof(sl[1].c_str());
			Args.ASSCoverageMulti[i][1]=atof(sl[4].c_str());
			Args.LSDRSs[i][0]=atof(sl[2].c_str());
			Args.LSDRSs[i][1]=atof(sl[5].c_str());
		}
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
    OH.addOpt('t', "threads", 1, "Number", "Number of threads. (8)",'i',&(Args.ThreadN));
    OH.addOpt('h', "help", 0, "", "Show this help and exit.",'b',&(Args.ShowHelp));
    OH.addOpt('v', "version", 0, "", "Show version and exit.",'b',&(Args.ShowVersion));
    // OH.addOpt(0, "NOH", 0, "", "No header, for test",'b',&(NoHeader));
    OH.addOpt(0, "CCS", 0, "", "Use default parameters for CCS data.",'b',&(Args.AllCCS));
    OH.addOpt(0, "CLR", 0, "", "Use default parameters for CLR data.",'b',&(Args.AllCLR));
    OH.addOpt(0, "DelClusterParas", 1, "Fixed,Ratio", "Custom clustering parameters for deletions.",'S',&(Args.CustomClusterParas[0]));
    OH.addOpt(0, "InsClusterParas", 1, "Fixed,Ratio", "Custom clustering parameters for insertions.",'S',&(Args.CustomClusterParas[1]));
    OH.addOpt(0, "DupClusterParas", 1, "Fixed,Ratio", "Custom clustering parameters for duplications.",'S',&(Args.CustomClusterParas[2]));
    OH.addOpt(0, "InvClusterParas", 1, "Fixed,Ratio", "Custom clustering parameters for inversions.",'S',&(Args.CustomClusterParas[3]));
    OH.addOpt(0, "DelFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for deletions.",'S',&(Args.CustomFilterParas[0]));
    OH.addOpt(0, "InsFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for insertions.",'S',&(Args.CustomFilterParas[1]));
    OH.addOpt(0, "DupFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for duplications.",'S',&(Args.CustomFilterParas[2]));
    OH.addOpt(0, "InvFilterParas", 1, "Base1,Ratio1,SDScore1,Base2,Ratio2,SDScore2", "Custom filter parameters for inversions.",'S',&(Args.CustomFilterParas[3]));
    OH.addOpt(0, "NOF", 0, "", "No filter, output all results.(default false)",'b',&(Args.NoFilter));
    OH.addOpt('m', 0, 1, "SVLEN", "Minimum SV length. (30)",'i',&(Args.MinSVLen));
    OH.addOpt('q', 0, 1, "Quality", "Minimum mapping quality. (20)",'i',&(Args.MinMappingQuality));
    OH.addOpt('l', 0, 1, "Length", "Minimum template length. (500)",'i',&(Args.MinTemplateLength));
    OH.addOpt('d', 0, 1, "Distance", "Minimum max merge distance of signature merging during CIGAR signature collection. (500)",'i',&(Args.DelMinMaxMergeDis));
    OH.addOpt('p', 0, 1, "Portion", "Max merge portion of signature merging during CIGAR signature collection. (0.2)",'F',&(Args.DelMaxMergePortion));
    OH.addOpt('c', 0, 1, "Size", "Coverage window size. (100)",'i',&(Args.CoverageWindowSize));
    OH.addOpt('M', 0, 1, "Size", "Max cluster size, will resize to this value if a cluster is larger than this. (1000)",'i',&(Args.MaxClusterSize));
    OH.addOpt(0, "InsClipTolerance", 1, "Size", "Insertion clip signature distance tolerance. (10)",'i',&(Args.InsClipTolerance));
    OH.addOpt(0, "InsMaxGapSize", 1, "Size", "Insertion clip signature max gap size. (50000)",'i',&(Args.InsMaxGapSize));
    OH.addOpt(0, "ClusteringBatchSize", 1, "Size", "Batch size of multihreading when clustering. (10000)",'i',&(Args.ClusteringBatchSize));
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

	omp_set_num_threads(8);
	// ThreadPool ThePool(8);

	updateTime("Starting kled, reading reference...");
	int NSeq;
	Contig * Contigs=getContigs(Args.ReferenceFileName,NSeq);//,RDWindowSize);
    
	updateTime("Reading reference","Getting stats...");
	vector<int> AllTechs=getAllTechs(Args);

	vector<Stats> AllStats=getAllStats(Args.ReferenceFileName,Args.BamFileNames,AllTechs);

	// for (int i=0;i<AllStats.size();++i) fprintf(stderr,"%f %f %f %f %f\n",AllStats[i].BelowIS,AllStats[i].MedianIS,AllStats[i].UpIS,AllStats[i].Mean,AllStats[i].SD);
	//exit(0);

	VCFHeader Header(Args.ReferenceFileName);
	addKledEntries(Header);
	for (int i=0;i<NSeq;++i)
	{
		if (Args.CallingContigs.size()!=0)
		{
			bool ToCall=false;
			for (int k=0;k<Args.CallingContigs.size();++k)
			{
				if (Contigs[i].Name==Args.CallingContigs[k])
				{
					ToCall=true;
					break;
				}
			}
			if (!ToCall) continue;
		}
		Header.addContig(Contigs[i]);
	}

	faidx_t * Ref=fai_load(Args.ReferenceFileName);
	// vector<vector<Variant>> VariantsByContig;
	bool FirstBam=true;
	vector<Sam> SamFiles=initSam(Args);
	double TotalCoverage=0;//Accumulative
	unsigned ProcessedLength=0;
	// FILE * WindowsFile=fopen("/home/cre/workspace/kled/data/wins.txt","wb");
	updateTime("Getting stats", "Starting calling...");
	for (int i=0;i<NSeq;++i)
	{
		if (Args.CallingContigs.size()!=0)
		{
			bool ToCall=false;
			for (int k=0;k<Args.CallingContigs.size();++k)
			{
				if (Contigs[i].Name==Args.CallingContigs[k])
				{
					ToCall=true;
					break;
				}
			}
			if (!ToCall) continue;
		}
		updateTime("","Calling...");
		SegmentSet AllPrimarySegments;
		vector<Signature> ContigTypeSignatures[NumberOfSVTypes];//For supported SV type
		unsigned int CoverageWindowSize=Args.CoverageWindowSize;
		unsigned int NumberOfCoverageWindows=Contigs[i].Size/CoverageWindowSize+1;
		double *CoverageWindows=new double[NumberOfCoverageWindows];
		for (int k=0;k<Contigs[i].Size/CoverageWindowSize+1;++k) CoverageWindows[k]=0;
		collectSignatures(Contigs[i],ContigTypeSignatures,AllPrimarySegments,Args,SamFiles,AllStats,AllTechs,CoverageWindows,0);
		AllPrimarySegments.sortNStat();
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
		double WholeCoverage=getAverageCoverage(0,Contigs[i].Size-1,CoverageWindows,Args, CoverageWindowsSums, CheckPoints, CheckPointInterval);
		// int NameLength=Contigs[i].Name.length();
		// fwrite(&(NameLength),sizeof(int),1,WindowsFile);
		// fwrite(Contigs[i].Name.c_str(),1,Contigs[i].Name.length(),WindowsFile);
		// fwrite(&(Contigs[i].Size),sizeof(unsigned),1,WindowsFile);
		// fwrite(&(NumberOfCoverageWindows),sizeof(unsigned),1,WindowsFile);
		// fwrite(CoverageWindows,sizeof(double),NumberOfCoverageWindows,WindowsFile);
		TotalCoverage=TotalCoverage*((double)(ProcessedLength)/(double)(ProcessedLength+Contigs[i].Size));
		TotalCoverage+=WholeCoverage*((double)(Contigs[i].Size)/(double)(ProcessedLength+Contigs[i].Size));
		ProcessedLength+=Contigs[i].Size;
		// double WholeCoverage=CoverageWindowsSums[(int)(Contigs[i].Size/CoverageWindowSize+1)]/(Contigs[i].Size/CoverageWindowSize+1);
		// continue;
		if (!NoHeader and FirstBam)
		{
			Header.addSample(Args.SampleName.c_str());
			printf(Header.genHeader(Args).c_str());
			FirstBam=false;
		}
		// int totalsig=0,cigardel=0, cigarins=0, cigardup=0, drpdel=0, drpdup=0, clipdel=0, clipins=0, clipdup=0, clipinv=0;
		// for (int m=0;m<NumberOfSVTypes;++m)
		// {
		// 	vector<Signature>& ContigSignatures=ContigTypeSignatures[m];
		// 	totalsig+=ContigSignatures.size();
		// 	for (int j=0;j<ContigSignatures.size();++j)
		// 	{
		// 		if (ContigSignatures[j].Type==0)
		// 		{
		// 			if (ContigSignatures[j].SupportedSV==0) ++cigardel;
		// 			if (ContigSignatures[j].SupportedSV==1) ++cigarins;
		// 			if (ContigSignatures[j].SupportedSV==2) ++cigardup;
		// 		}
		// 		else if (ContigSignatures[j].Type==1)
		// 		{
		// 			if (ContigSignatures[j].SupportedSV==0) ++drpdel;
		// 			if (ContigSignatures[j].SupportedSV==2) ++drpdup;
		// 		}
		// 		else
		// 		{
		// 			if (ContigSignatures[j].SupportedSV==0) ++clipdel;
		// 			if (ContigSignatures[j].SupportedSV==1) ++clipins;
		// 			if (ContigSignatures[j].SupportedSV==2) ++clipdup;
		// 			if (ContigSignatures[j].SupportedSV==3) ++clipinv;
		// 		}
		// 	}
		// }
		// fprintf(stderr,"%s: %d\n, cigardel: %d, cigarins: %d, cigardup: %d, drpdel: %d, drpdup: %d, clipdel: %d, clipins: %d, clipdup: %d, clipinv: %d. Contig Size:%ld, Average Coverage: %lf, Total Average Coverage: %lf\n",Contigs[i].Name.c_str(),totalsig,cigardel, cigarins, cigardup, drpdel, drpdup, clipdel, clipins, clipdup, clipinv, Contigs[i].Size, WholeCoverage, TotalCoverage);

		updateTime("Getting signatures","Clustering...");
		vector<vector<Signature>> SignatureTypeClusters[NumberOfSVTypes];
		for (int k=0;k<NumberOfSVTypes;++k)
		{
			sortAndDeDup(ContigTypeSignatures[k]);
			clustering(ContigTypeSignatures[k],SignatureTypeClusters[k],AllStats[i],Args);
		}
		updateTime("Clustering","Generating results...");
		vector<vector<Signature>> SignatureClusters;
		for (int k=0;k<NumberOfSVTypes;++k) SignatureClusters.insert(SignatureClusters.end(),make_move_iterator(SignatureTypeClusters[k].begin()),make_move_iterator(SignatureTypeClusters[k].end()));
		vector<VCFRecord> Records;
		#pragma omp parallel for reduction(RecordVectorConc:Records)
		for (int j=0;j<SignatureClusters.size();++j)
		{
			// ++Times[omp_get_thread_num()];
			Records.push_back(VCFRecord(Contigs[i],Ref,SignatureClusters[j], AllPrimarySegments,CoverageWindows, TotalCoverage, Args, CoverageWindowsSums, CheckPoints, CheckPointInterval));
		}
		updateTime("Results generation","Sorting results...");
		sort(Records.data(),Records.data()+Records.size());
		// Records.sort();
		updateTime("Results sorting","Outputing results...");
		unsigned int SVCounts[NumberOfSVTypes];for (int i=0;i<NumberOfSVTypes;++i) SVCounts[i]=0;
		for (auto r: Records)
		{
			if (!r.Keep) continue;
			// r.genotype(Contigs[i],AllPrimarySegments,CoverageWindows,CoverageWindowsSums,CheckPoints,CheckPointInterval,Args);
			r.resolveRef(Contigs[i],Ref,SVCounts[r.getSVTypeI()],Args);
			++SVCounts[r.getSVTypeI()];
			printf("\n%s",string(r).c_str());
		}
		delete CoverageWindows;
		// free(CoverageWindowsSums);
		// free(CheckPoints);
	}

	//report(VariantsByContig);
	fai_destroy(Ref);
	free(Contigs);
	closeSam(SamFiles);
	// fclose(WindowsFile);
    return 0;
}