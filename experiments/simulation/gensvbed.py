import sys
import random
import pysam

CHRPrefix="chr"

Alphabet="ATGC"
def genInsBases(SVLen):
    Base=""
    for i in range(SVLen):
        Base+=Alphabet[random.randrange(4)]
    return Base

def genAF(Min=0.01):
    return 0.4#Portion of 1/1=0.25

class SV:
    def __init__(self,Type="", Chr="",Start=0, End=0, SVLen=0, AF=0, ID=0):
        self.Type=Type
        self.Start=Start
        self.End=End
        self.Chr=Chr
        self.AF=AF
        self.ID=ID
        if SVLen==0:
            self.SVLen=End-Start
        else:
            self.SVLen=SVLen
        if (Type=="insertion"):
            self.genIns()
        elif (Type=="tandem duplication"):
            self.genCN()
    def genIns(self):
        self.InsBases=genInsBases(self.SVLen)
    def genCN(self):
        CN=random.gauss(2,2)
        CN=int(CN)
        if (CN<2):
            CN=2
        self.CN=CN
    def __lt__(self,other):
        if (self.Chr==other.Chr):
            return self.Start<other.Start
        try:
            ChrN=int(self.Chr)
            OChrN=int(other.Chr)
            return ChrN<OChrN
        except:
            return self.Chr<other.Chr
    def __eq__(self,other):
        return self.ID==other.ID
    def __str__(self):
        Str="%s\t"%(self.Chr if len(self.Chr)>2 else CHRPrefix+self.Chr)
        if (self.Type=="insertion"):
            if self.Start==1:
                self.Start=2
            Str+="%s\t%s\t"%(self.Start-1,self.Start)
        else:
            Str+="%s\t%s\t"%(self.Start,self.End-1)
        Str+="%s\t"%(self.Type)
        if (self.Type=="tandem duplication"):
            Str+="%s\t"%(self.CN)
        elif (self.Type=="insertion"):
            Str+=self.InsBases+"\t"
        else:
            Str+="None\t"
        Str+="0"
        return Str

SVIndelPosDist={}
TotalSVIndel=0
SVNearRange=5000

def getSVIndelPos():
    global SVIndelPosDist
    global TotalSVIndel
    Pos=random.randint(0,max(0,TotalSVIndel-1))
    Chr=""
    for i in SVIndelPosDist.keys():
        if Pos<len(SVIndelPosDist[i]):
            Chr=i
            break
        Pos-=len(SVIndelPosDist[i])
    if Chr=="":
        return -1,-1
    PosR=random.random()
    PosP=SVIndelPosDist[Chr][random.randint(0,len(SVIndelPosDist[Chr])-1)]
    Pos=float(ChromLengths[Chr])*PosP
    Pos=(1.0-PosR)*max(0,Pos-SVNearRange)+PosR*min(ChromLengths[Chr]-1,Pos+SVNearRange)
    Pos=int(Pos)
    return Chr,Pos

def getSVLenDist(VCFFileName,Bed=None):
    SVLenDist={"DEL":[],"INS":[]}
    global SVIndelPosDist
    global TotalSVIndel
    Included=None
    if Bed!=None:
        Included={}
        BedFile=open(Bed,"r")
        for line in BedFile:
            sl=line.split()
            if (sl[0] not in Included.keys()):
                Included[sl[0]]=[]
            Included[sl[0]].append((int(sl[1]),int(sl[2])))
        BedFile.close()
    vcf=pysam.VariantFile(VCFFileName,"r")
    for record in vcf.fetch():
        if Included!=None: 
            pos=record.pos
            chr=record.chrom
            go=False
            for interval in Included[chr]:
                if pos>=interval[0] and pos<=interval[1]:
                    go=True
                    break
            if not go:
                continue
        SVTYPE=record.info["SVTYPE"]
        chrom="chr"+record.chrom
        for SVLEN in record.info["SVLEN"]:
            if chrom not in SVIndelPosDist.keys():
                SVIndelPosDist[chrom]=[]
            # if abs(SVLEN)<50:
            #     continue
            if (SVTYPE=="DEL"):
                SVLenDist["DEL"].append(abs(int(SVLEN)))
                SVIndelPosDist[chrom].append(float(record.pos)/float(vcf.header.contigs.get(record.chrom).length))
            elif (SVTYPE=="INS"):
                SVLenDist["INS"].append(abs(int(SVLEN)))
                SVIndelPosDist[chrom].append(float(record.pos)/float(vcf.header.contigs.get(record.chrom).length))
            else:
                continue
    vcf.close()
    SVLenDist["DEL"]=sorted(SVLenDist["DEL"])
    SVLenDist["INS"]=sorted(SVLenDist["INS"])
    for chrom in SVIndelPosDist.keys():
        SVIndelPosDist[chrom]=sorted(SVIndelPosDist[chrom])
        TotalSVIndel+=len(SVIndelPosDist[chrom])
    return SVLenDist

def deOverlap(SVs):
    SVs.sort()
    i=0
    while i< len(SVs)-1:
        j=i+1
        while j<len(SVs):
            if SVs[i].Chr!=SVs[j].Chr:
                break
            if (SVs[i].End>=SVs[j].Start-1):
                del SVs[j]
                continue
            break
        i+=1

def genHaplotypes(SVs,H1SVs,H2SVs):
    for V in SVs:
        P=V.AF*V.AF/(1.0-(1.0-V.AF)**2)
        # P=0.25
        if random.uniform(0,1)<P:#homo
            H1SVs.append(V)
            H2SVs.append(V)
        else:
            if random.uniform(0,1)<0.5:
                H1SVs.append(V)
            else:
                H2SVs.append(V)

def removeDuplicants(SVs):
    SVs=SVs.copy()
    SVs.sort()
    i=0
    while i< len(SVs)-1:
        if (i+1<len(SVs) and SVs[i]==SVs[i+1]):
            del SVs[i+1]
            continue
        i+=1  
    return SVs

ChromLengths={}
def getChrs(FaFileName):
    Chroms=[]
    TotalLength=0
    fa=pysam.FastaFile(FaFileName)
    for chrom in fa.references:
        Chroms.append((chrom,fa.get_reference_length(chrom)))
        ChromLengths[chrom]=fa.get_reference_length(chrom)
        TotalLength+=fa.get_reference_length(chrom)
    return Chroms,TotalLength

def getPosChr(Pos,Chroms):
    for i in range(len(Chroms)):
        if Pos<Chroms[i][1]:
            return Chroms[i][0],Pos
        Pos-=Chroms[i][1]
    return -1,-1

def find(SortedArray, Value):
    for i in range(len(SortedArray)):
        if Value<=SortedArray[i]:
            return i
    return len(SortedArray)

#calc from HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz
IndelLengths=[0, 296893, 113890, 39432, 70475, 16385, 20240, 4573, 17304, 4094, 8254, 2170, 9189, 1635, 3513, 2091, 4017, 939, 2054, 774, 2353, 666, 1018, 502, 1362, 492, 623, 348, 658, 259, 480, 221, 407, 168, 231, 153, 274, 142, 175, 118, 198, 118, 136, 65, 93, 55, 88, 53, 99, 59]
IndelLengthAcc=[]
IndelN=0
for i in IndelLengths:
    IndelN+=i
    IndelLengthAcc.append(IndelN)
def getRandomIndelLength():
    r=random.randint(0,IndelN)
    for i in range(1,len(IndelLengthAcc)):
        if r<=IndelLengthAcc[i]:
            return i
    return 49

ComplexRegionRatio=0.1
ComplexRegionRange=250
ComplexRetionIndel=0.1
TypeNames=["DEL","INS","DUP","INV","INDEL"]
TypeLongNames=["deletion","insertion","tandem duplication","inversion","indel"]
def genSVs(Chroms,TotalLength, SVLenDist,Numbers=(70000,70000,1000,1000,700000)):
    SVs=[]
    Indels=[]
    ComplexRegions=[]
    for Typei in range(5):
        Type=TypeNames[Typei]
        if Type!="INDEL":
            for i in range(Numbers[Typei]):
                RandI=random.randint(0,len(SVLenDist[Type])-1)
                SVLEN=SVLenDist[Type][RandI]
                if RandI!=len(SVLenDist[Type])-1:
                    RandR=random.random()
                    SVLEN=SVLEN*RandR+SVLenDist[Type][RandI+1]*(1-RandR)
                SVLEN=SVLEN+((-1)**random.randint(0,1))*SVLEN*random.random()*0.1
                SVLEN=int(SVLEN)
                if SVLEN<50:
                    SVLEN=50
                AF=genAF()
                GPos=random.randint(0,TotalLength-1)
                Chr,Pos=getPosChr(GPos,Chroms)
                if Chr==-1:
                    raise -1
                if ChromLengths[Chr]<Pos+SVLEN:
                    Pos=ChromLengths[Chr]-1-SVLEN
                if Pos<1:
                    Pos=1
                SVs.append(SV(TypeLongNames[Typei],Chr,Pos,Pos+SVLEN,SVLEN,AF,"SIMULATED.%s.%s"%(Type,i)))
                if random.random()<ComplexRegionRatio:
                    ComplexRegions.append((Chr,max(0, Pos-ComplexRegionRange),min(ChromLengths[Chr],(Pos if Typei==1 else Pos+SVLEN)+ComplexRegionRange)))
        else:
            for i in range(Numbers[Typei]):
                SVLEN=getRandomIndelLength()
                AF=genAF()
                GPos=random.randint(0,TotalLength-1)
                Chr,Pos=getPosChr(GPos,Chroms)
                if Chr==-1:
                    raise -1
                if ChromLengths[Chr]<Pos+SVLEN:
                    Pos=ChromLengths[Chr]-1-SVLEN
                if Pos<1:
                    Pos=1
                InsDel=random.randint(0,1)
                Indels.append(SV(TypeLongNames[InsDel],Chr,Pos,Pos+SVLEN,SVLEN,AF,"SIMULATED.%s.%s.%s"%(Type,TypeNames[InsDel],i)))
            for r in ComplexRegions:
                AF=0.7
                Chr=r[0]
                Start=r[1]
                Stop=r[2]
                for i in range(int((Stop-Start)*ComplexRetionIndel)):
                    SVLEN=getRandomIndelLength()
                    Pos=random.randint(Start,Stop-1)
                    if Stop<Pos+SVLEN:
                        Pos=Stop-1-SVLEN
                    if Pos<1:
                        Pos=1
                    InsDel=random.randint(0,1)
                    Indels.append(SV(TypeLongNames[InsDel],Chr,Pos,Pos+SVLEN,SVLEN,AF,"SIMULATED.%s.%s.%s"%(Type,TypeNames[InsDel],i)))
    return SVs,Indels

import copy
def fixSVPos(Indels, SVs):
    All=SVs+Indels
    All=copy.deepcopy(All)
    All=sorted(All)
    Chr=""
    Acc=0
    ChrLength=0
    Result=[]
    for i in range(len(All)):
        V=All[i]
        if Chr!=V.Chr:
            if Chr!="":
                ChrLength+=Acc
                print("%s : Original: %s, After indel: %s"%(Chr,ChromLengths[Chr],ChrLength),file=sys.stderr)
                j=i-1
                while j>=0:
                    V2=All[j]
                    if V2.Chr!=Chr:
                        break
                    if V2.Type=="insertion":
                        if V2.Start>=ChrLength:
                            V2.Start=max(1,ChrLength-1)
                    elif V2.End>ChrLength:
                        V2.Start-=ChrLength-2-V2.End
                        if V2.Start<1:
                            V2.Start=1
                        V2.End=ChrLength-1
                    j-=1 
            Chr=V.Chr
            ChrLength=ChromLengths[Chr]
            Acc=0
        if V.SVLen<50:
            if V.Type=="insertion":
                Acc+=V.SVLen
            else:
                Acc-=V.SVLen
        else:
            V.Start+=Acc
            V.End+=Acc
            Result.append(V)
    return Result

import sys
def run():
    TotalSVN=25000
    Chroms, TotalLength=getChrs(sys.argv[1])
    Bed=None
    if len(sys.argv)>6:
        Bed=sys.argv[6]
    SVLenDist=getSVLenDist(sys.argv[2],Bed)
    H1BedName=sys.argv[3]
    H2BedName=sys.argv[4]
    Seed=int(sys.argv[5])
    random.seed(Seed)
    print("Seed: %s."%Seed, file=sys.stderr)
    MinDelLen=50
    MinInsLen=50
    MinDupLen=400
    MinInvLen=200
    MaxInvLen=10000
    SVLenDist["DUP"]=SVLenDist["INS"][find(SVLenDist["INS"],MinDupLen):]
    SVLenDist["INV"]=SVLenDist["DEL"][find(SVLenDist["DEL"],MinInvLen):find(SVLenDist["DEL"],MaxInvLen)]
    SVLenDist["DEL"]=SVLenDist["DEL"][find(SVLenDist["DEL"],MinDelLen):]
    SVLenDist["INS"]=SVLenDist["INS"][find(SVLenDist["INS"],MinInsLen):]
    SVs,Indels=genSVs(Chroms,TotalLength,SVLenDist)
    SVs=random.sample(SVs,TotalSVN)
    # SVs+=Indels

    H1SVs=[]
    H2SVs=[]
    genHaplotypes(SVs,H1SVs,H2SVs)
    deOverlap(H1SVs)
    deOverlap(H2SVs)
    AllSVs=removeDuplicants(H1SVs+H2SVs)

    
    H1Indels=[]
    H2Indels=[]
    genHaplotypes(Indels,H1Indels,H2Indels)
    deOverlap(H1Indels)
    deOverlap(H2Indels)
    H1IndelFile=open(H1BedName+".indel","w")
    for V in H1Indels:
        print(V,file=H1IndelFile)
    H1IndelFile.close()
    H2IndelFile=open(H2BedName+".indel","w")
    for V in H2Indels:
        print(V,file=H2IndelFile)
    H2IndelFile.close()

    NDEL=0
    NINS=0
    NDUP=0
    NINV=0
    NSINS=0
    NSDEL=0

    H1NSV=0
    H2NSV=0
    for V in AllSVs+H1Indels+H2Indels:
    # for V in AllSVs:
        if V.Type=="deletion":
            if V.SVLen>=50:
                NDEL+=1
            else:
                NSDEL+=1
        elif V.Type=="insertion":
            if V.SVLen>=50:
                NINS+=1
            else:
                NSINS+=1
        elif V.Type=="tandem duplication":
            NDUP+=1
        elif V.Type=="inversion":
            NINV+=1
    NSV=NDEL+NINS+NDUP+NINV
    for V in H1SVs:
        if V.SVLen>=50:
            H1NSV+=1
    for V in H2SVs:
        if V.SVLen>=50:
            H2NSV+=1

    print("All SVs: %d, SVs for H1: %d, SVs for H2: %d. DEL, INS, DUP and INV: %d %d %d %d. Indels: %s, %s DEL, %s INS."%(NSV,H1NSV,H2NSV,NDEL,NINS,NDUP,NINV,NSDEL+NSINS,NSDEL,NSINS),file=sys.stderr)

    H1SVIDs=set()
    H2SVIDs=set()
    for V in H1SVs:
        H1SVIDs.add(V.ID)
    for V in H2SVs:
        H2SVIDs.add(V.ID)
    
    for V in AllSVs:
        GT=""
        if (V.ID in H1SVIDs):
            GT+="1"
        else:
            GT+="0"
        GT+="|"
        if (V.ID in H2SVIDs):
            GT+="1"
        else:
            GT+="0"
        print(str(V)+"\t"+GT+"\t"+V.ID)

    H1SVs=fixSVPos(H1Indels,H1SVs)
    H2SVs=fixSVPos(H2Indels,H2SVs)

    H1F=open(H1BedName,"w")
    for V in H1SVs:
        print(V,file=H1F)
        # H1SVIDs.add(V.ID)
    H1F.close()
    H2F=open(H2BedName,"w")
    for V in H2SVs:
        print(V,file=H2F)
        # H2SVIDs.add(V.ID)
    H2F.close()


if __name__=="__main__":
    run()
