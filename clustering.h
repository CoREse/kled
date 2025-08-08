#ifndef KLED_CLUSTERING_H
#define KLED_CLUSTERING_H

#include "signature.h"
#include <vector>
#include "input.h"
#include "kled.h"

struct ClusterCore
{
    int Begin,End;
    ClusterCore(int B=0,int E=0):Begin(B),End(E){}
};

float distance(const Signature &A, const Signature &B, bool Partial=true, float *PPD=NULL);
int precisionLevel(const Signature &A);
int bestPrecision(const Signature &A,const Signature &B);
int worstPrecision(const Signature &A,const Signature &B);
void clustering(int SVTypeI, std::string & ContigName, std::vector<Signature> & SortedSignatures, std::vector<std::vector<Signature>> &Clusters, std::vector<ClusterCore> &Cores, Stats BamStats, Arguments& Args);

#endif