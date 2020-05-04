#!usr/bin/python
## Tcr Repertoire Utilities for Solid Tumors (TRUST)

##-------------------------------------------------------
# Repurposing non-selected RNA-sequencing data to analyze
# T cell receptor repertoire in solid tumors or tissues
# All rights preserved by Bo Li: bli@jimmy.harvard.edu
# Sep 30, 2015
# Updates: May 31, 2016, add single end component
# Updates: Nov 08, election day, 2016: add multiple-variable gene assignment
#               force TRUST to screen for all vgene assignment, slow
##-------------------------------------------------------

##-------------------------------------------------------
# Update: 2016-10-31 making this compatible with VIPER
# 1. fixing -o option
# 2. setting correct static file directory
# Update: 2016-11-14 updating the TRUST2.1 script w/ these changes
# TRUST _version# 2.1
##-------------------------------------------------------

import sys,re,os,resource,subprocess
try:
    import pysam
except:
    print '''Please install pysam first: https://code.google.com/p/pysam/'''
    sys.exit()
import numpy as np
import random,time
from itertools import chain
from optparse import OptionParser
from copy import deepcopy
from collections import defaultdict
from Bio import pairwise2
pairwise2.MAX_ALIGNMENTS=1
t0=time.time()
curDir=re.sub('TRUST2.py','',os.path.dirname(os.path.realpath(sys.argv[0])))+'/'
#curDir='./'

#ASSUME it's being run in the VIPER PROJECT directory
#MAYBE this should be a command line option in the future!!!
_static_path_dir=os.path.join(os.getcwd(), "viper", "static", "cdr3")

try:
    from gmpy import popcount as PopC
    print "Successfully import gmpy"
except:
    PopC = lambda x:bin(x).count('1')

## Prepare prerequisite data
BinCodeDict={'A':(0,0),'T':(1,1),'C':(1,0),'G':(0,1)}
DNACodeDict={'00':'A','11':'T','10':'C','01':'G'}
sys.setrecursionlimit(10000)
RCdict={'A':'T','T':'A','C':'G','G':'C','N':'N'}
AAcode={
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    }

#CDR3patV={"V0":'CASS',"V10_3":"CAIS","V12_5":"CASG","V15":"CATS","V20_1":"CASR",
#          "V24_1":"CATS","V29_1":"CSVE","V30":"CAWS"}
#CDR3patJ={"J1_1":"FGQG","J1_2":"FGSG","J1_3":"FGEG","J1_4":"FGSG",
#          "J1_5":"FGDG","J1_6":"FGNG","J2_1":"FGPG","J2_2":"FGEG",
#          "J2_2P":"LGGG","J2_3":"FGPG","J2_4":"FGAG","J2_5":"FGPG",
#          "J2_6":"FGAG","J2_7_01":"FGPG","J2_7_02":"VGPG"}
CDR3patV={"V0":'CASS',"V1":"CA.{1}S","V2":"CAS","V3":"CSVE",
          "V4":"CA[VGMLFYE]{1}","V5":"CVV","V6":"CI.{1}R",
          "V7":"Y[LF]{1}CA","V8":"CA[TL]{1}W","V9":"RASS",
          "V10":"CSA","V11":"SQTS",'V12':'CTSS','V13':'FCA',
          'V14':'YCA','V15':'RKSS','V16':'FSLQI','V17':'YFC',
          'V18':'YYCLL','V19':'DTTLN','V20':'SSLYL','V21':'SATYL',
          'V22':'GLEEK','V23':'KPSVQ','V24':'AMYF'}
CDR3patJ={"J0":'FG.{1}G',"J1":"LGGG","J2":"VGPG","J3":"WG.G","J4":"FA.G"}

CDR3patVB={'VB1':"YYCA",'VB2':"YYCV",'VB3':'YYCT','VB4':'DTA[TV]{1}YY','VB5':'RSDDT',
           'VB6':"LSSLR",'VB7':'DTSKN','VB8':'DSKNS','VB9':'LKAED','VB10':'KASDT',
           'VB11':'YHCA','VB12':'DYCA','VB13':'VYCC'}
CDR3patJB={'JB':"WG.G"}

patDict={'BJforward':re.compile('[T,C]{1}TT[T,C]{1}GG[A,T,G,C]{4}GG'),
         'BJreverse':re.compile('CC[A,T,G,C]{4}CC[A,G]{1}AA[A,G]{1}'),
         'BV1forward':re.compile('[C,T,G]{1}TG[C,T]{1}GCC'),
         'BV2forward':re.compile('CT[G,A]{1}CAG[C,T]{1}'),
         'BV3forward':re.compile('GCCAG[C,T,A]{1}A'),
         'BV1reverse':re.compile('GGC[G,A]{1}CA[G,A,C]{1}'),
         'BV2reverse':re.compile('[G,A]{1}CTG[C,T]{1}AG'),
         'BV3reverse':re.compile('T[G,A,T]{1}CTGGC'),
         'DJforward':re.compile('T[C,T]{1}TT[C,T]{1}GG[A,C]{1}A[A,C]{1}[A,G,T]{1}GG'),
         'DJreverse':re.compile('CC[T,C,A]{1}[T,G]{1}T[T,G]{1}CC[G,A]{1}AA[G,A]{1}A'),
         'DVforward':re.compile('TACT[A,C,T]{2}TGTGC'),
         'DJreverse':re.compile('GCACA[T,G,A]{2}AGTA')}

patDictAll={'AJforward':re.compile('T[TG]{1}.G[TC]{1}....GG.A'),
            'AJreverse':re.compile('T.CC....[AG]{1}C.[AC]{1}A'),
            'AV1forward':re.compile('A...AG..[AC]{1}.[CT]{1}TA'),
            'AV2forward':re.compile('T.TA[CT]{1}T[AT]{1}CTG'),
            'AV3forward':re.compile('C..C....TA[CT]{1}.[AT]{1}T'),
            'AV4forward':re.compile('TGTT.TG'),
            'AV1reverse':re.compile('TA[GA]{1}.[TG]{1}..CT...T'),
            'AV2reverse':re.compile('CAG[AT]{1}A[GA]{1}TA.A'),
            'AV3reverse':re.compile('A[AT]{1}.[GA]{1}TA....G..G'),
            'AV4reverse':re.compile('CA.AACA'),
         'BJforward':re.compile('[TC]{1}TT[TC]{1}GG....GG'),
         'BJreverse':re.compile('CC....CC[AG]{1}AA[AG]{1}'),
         'BV1forward':re.compile('[CTG]{1}TG[CT]{1}GCC'),
         'BV2forward':re.compile('CT[GA]{1}CAG[CT]{1}'),
         'BV3forward':re.compile('GCCAG[CTA]{1}A'),
         'BV1reverse':re.compile('GGC[GA]{1}CA[GAC]{1}'),
         'BV2reverse':re.compile('[GA]{1}CTG[CT]{1}AG'),
         'BV3reverse':re.compile('T[GAT]{1}CTGGC'),
         'GJforward':re.compile('AA[AG]{1}[CA]{1}[TAC]{1}[CAG]{1}T[CT]{1}[AT]{1}[AG]{1}[GC]{1}'),
         'GJreverse':re.compile('[GC]{1}[TC]{1}[A,T]{1}[G,A]{1}A[G,T,C]{1}[A,T,G]{1}[G,T]{1}[T,C]{1}TT'),
         'GVforward':re.compile('[AG]{1}[TGA]{1}G[CAG]{1}C[T,A,C]{1}[T,C]{1}[G,A,C]{1}[T,A,C]{1}G[G,A]{1}G'),
         'GVreverse':re.compile('C[CT]{1}C[GTA]{1}[G,C,T]{1}[A,G]{1}[A,T,G]{1}G[G,T,C]{1}C[A,C,T]{1}[T,C]{1}'),
         'DJforward':re.compile('T[CT]{1}TT[CT]{1}GG[A,C]{1}A[A,C]{1}[A,G,T]{1}GG'),
         'DJreverse':re.compile('CC[TCA]{1}[TG]{1}T[T,G]{1}CC[G,A]{1}AA[G,A]{1}A'),
         'DVforward':re.compile('TACT[ACT]{2}TGTGC'),
         'DJreverse':re.compile('GCACA[TGA]{2}AGTA')}

patDictB={'HJforward':re.compile('CTGGGG[CG]{1}[CA]{1}'),
        'HJreverse':re.compile('[GT]{1}[CG]{1}CCCCAG'),
        'HVforward':re.compile('TATTACTGT'),
        'HVreverse':re.compile('ACAGTAATA')}

def ParseFa(fname):
    InputStr=open(fname).readlines()
    FaDict={}
    seq=''
    for line in InputStr:
        if line.startswith('>'):
            if len(seq)>0:
                FaDict[seqHead]=seq
                seq=''
            seqHead=line.strip()
        else:
            seq+=line.strip()
    if seqHead not in FaDict:
        FaDict[seqHead]=seq
    return FaDict

AVFaDict=ParseFa(os.path.join(_static_path_dir,"TRAV-imgt-AA.fa"))
AJFaDict=ParseFa(os.path.join(_static_path_dir,"TRAJ-imgt-AA.fa"))
BVFaDict=ParseFa(os.path.join(_static_path_dir,"TRBV-imgt-AA.fa"))
BJFaDict=ParseFa(os.path.join(_static_path_dir,"TRBJ-imgt-AA.fa"))
DVFaDict=ParseFa(os.path.join(_static_path_dir,"TRDV-imgt-AA.fa"))
DJFaDict=ParseFa(os.path.join(_static_path_dir,"TRDJ-imgt-AA.fa"))
GVFaDict=ParseFa(os.path.join(_static_path_dir,"TRGV-imgt-AA.fa"))
GJFaDict=ParseFa(os.path.join(_static_path_dir,"TRGJ-imgt-AA.fa"))

AVFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRAV-imgt-DNA.fa"))
AJFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRAJ-imgt-DNA.fa"))
BVFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRBV-imgt-DNA.fa"))
BJFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRBJ-imgt-DNA.fa"))
DVFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRDV-imgt-DNA.fa"))
DJFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRDJ-imgt-DNA.fa"))
GVFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRGV-imgt-DNA.fa"))
GJFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"TRGJ-imgt-DNA.fa"))

bHVFaDict=ParseFa(os.path.join(_static_path_dir,"IGHV-imgt-AA.fa"))
bHJFaDict=ParseFa(os.path.join(_static_path_dir,"IGHJ-imgt-AA.fa"))
bHVFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"IGHV-imgt-DNA.fa"))
bHJFaDict_DNA=ParseFa(os.path.join(_static_path_dir,"IGHJ-imgt-DNA.fa"))

def GetJMotifs(FaDict,PAT=CDR3patJ):
        JMotifs={}
        for kk in FaDict:
            jj=FaDict[kk]
            for pp in PAT.values():
                mm=re.search(pp,jj)
                if mm is not None:
                    break
            mms=mm.span()
            kkJ=jj[mms[0]-3:mms[1]]
            if kkJ not in JMotifs:
                JMotifs[kkJ]=[kk]
            else:
                JMotifs[kkJ].append(kk)
        return JMotifs

def GetVMotifs(FaDict,PAT=CDR3patV):
        VMotifs={}
        VmotifList=[]
        for kk in FaDict:
            vv=FaDict[kk]
            for pp in PAT.values():
                mm=re.search(pp,vv)
                if mm is not None:
                    break
            mms=mm.span()
            kkV=vv[mms[0]-3:mms[1]] ## Keep at least 3 flanking amino acids
            VmotifList.append(kkV)
        for kk in FaDict:
            vv=FaDict[kk]
            for kkV in VmotifList:
                if kkV in vv:
                    if kkV not in VMotifs:
                        VMotifs[kkV]=[kk]
                    else:
                        VMotifs[kkV].append(kk)
        for kk in VMotifs:
            vv=VMotifs[kk]
            vv=list(set(vv))
            VMotifs[kk]=vv
        return VMotifs

MotifDict={'BJ':GetJMotifs(BJFaDict),'BV':GetVMotifs(BVFaDict),
           'DJ':GetJMotifs(DJFaDict),'DV':GetVMotifs(DVFaDict),
           'AJ':GetJMotifs(AJFaDict),'AV':GetVMotifs(AVFaDict),
           'GJ':GetJMotifs(GJFaDict),'GV':GetVMotifs(AVFaDict)}
DNAFaDict={'BJ':BJFaDict_DNA,'BV':BVFaDict_DNA,
           'DJ':DJFaDict_DNA,'DV':DVFaDict_DNA,
           'AJ':AJFaDict_DNA,'AV':AVFaDict_DNA,
           'GJ':GJFaDict_DNA,'GV':GVFaDict_DNA}
MotifDictB={'HJ':GetJMotifs(bHJFaDict,CDR3patJB),'HV':GetVMotifs(bHVFaDict,CDR3patVB)}
DNAFaDictB={'HJ':bHJFaDict_DNA,'HV':bHVFaDict_DNA}

def BuildPrefix(seq,MIN=1,direction=1):
	ns=len(seq)
	if ns<=MIN:
		print "sequence shorter than minimum length"
		raise
	if direction==1:
		PreSeq=[]
		for i in xrange(MIN,ns):
			PreSeq.append(seq[0:i])
	if direction==2:
		PreSeq=[]
		for i in xrange(0,ns-MIN):
			PreSeq.append(seq[i:ns])
	return PreSeq

def ReverseCompSeq(seq):
	seq_rc=''
	seqL=list(seq)
	for ss in seqL:
		seq_rc=RCdict[ss]+seq_rc
	return seq_rc

def FindOverlapSeq(Prefix1,Prefix2):
        tmp=list(set(Prefix1) & set(Prefix2))
        if len(tmp)>1:
            return max(tmp,key=len)
        if len(tmp)==1:
            return tmp[0]
        if len(tmp)==0:
            return ''
            
def CompareSuffixSeq(seq1,seq2,err=1):
        seq1c=ReverseCompSeq(seq1)
        x=ConvertDNAtoBinary(seq1)
        y=ConvertDNAtoBinary(seq2)
        xc=ConvertDNAtoBinary(seq1c)
        n1=len(seq1)
        n2=len(seq2)
	suffList=[]
	suf=CompareSuffixByBitSeq(x,y,n1,n2,err=err)
	suffList.append(suf)
	suf=CompareSuffixByBitSeq(y,x,n2,n1,err=err)	
	suffList.append(suf)
	suf=CompareSuffixByBitSeq(xc,y,n1,n2,err=err)
	suffList.append(suf)
	suf=CompareSuffixByBitSeq(y,xc,n2,n1,err=err)	
	suffList.append(suf)  
        MaxL=-1
        MaxOs=-1
        MaxErr=0
        for i in xrange(0,4):
            tmp=suffList[i]
            if tmp[0]>MaxL:
                MaxL=tmp[0]
                MaxOs=i
                MaxErr=tmp[1]
        return ((MaxL,MaxErr),MaxOs)

def ReverseComp(x,n=50):
	x0=2**n-1-x[0]
	x1=2**n-1-x[1]
	s0=bin(x0)[2:].zfill(n)
	x0c=int(s0[::-1],2)
	s1=bin(x1)[2:].zfill(n)
	x1c=int(s1[::-1],2)
#	s0=bin(x0)
#	d0='0'*(n+2-len(s0))
#	s1=bin(x1)
#	d1='0'*(n+2-len(s1))
#	x0c=int(s0[:1:-1]+d0,2)
#	x1c=int(s1[:1:-1]+d1,2)
	return (x0c,x1c)

def ConvertDNAtoBinary(seq):
	binArray=[]
	for s in list(seq):
		#if s=='N':
		#	s=BinCodeDict.keys()[random.randint(0,3)]
		binArray.append(BinCodeDict[s])
	a1=''
	b1=''
	for k in binArray:
		a1+=str(k[0])
		b1+=str(k[1])
	return (int(a1,2),int(b1,2))

def ConvertBitToDNA(x,n=50):
	a0=format(x[0],'0'+str(n)+'b')
	a1=format(x[1],'0'+str(n)+'b')
	seq=''
	for i in xrange(0,n):
		seq+=DNACodeDict[a0[i]+a1[i]]
	return seq

def CompareSuffixByBitSeq(x,y,nx,ny,err=1):
        if x[0]==y[0] and x[1]==y[1]:
            return (nx,0)
        err_min=[999,999]
        suffL=[-1,-1]
        bitIdx=[-1,-1]
        nm=min(nx,ny)
        for ss in xrange(0,2):
            xx=x[ss]
            yy=y[ss]
            OBJmin=9999
            for kk in xrange(0,nm-2):
                yyk=yy>>(ny-kk)
                xxk0=xx>>kk
                xxk=xx-(xxk0<<kk)
                tmp=xxk^yyk
                ee=PopC(tmp)
                obj=ee*3-kk
                if obj<=OBJmin:
                    OBJmin=obj
                    err_min[ss]=ee
                    suffL[ss]=kk
                    bitIdx[ss]=tmp
	if bitIdx[0]==bitIdx[1] and bitIdx[0]>0:
		err_t=err_min[0]
	else:
		err_t=err_min[0]+err_min[1]
	if err_t<= err and suffL[0]==suffL[1]:
		return (suffL[0],err_t)
	else:
		return (0,0)           

def CompareSuffixByBit(x,y,n=50,err=1):
	## find the overlap sequence between x and y, tolerate up to err mismatches, currently gap is not supported
	## input is ordered. assuming x is in front of y
	## x and y are two-bit coded, i.e. tuples with two integers
	## n is the read length
        if x[0]==y[0] and x[1]==y[1]:
            return (n,0)
	err_min=[999,999]
	suffL=[-1,-1]
	bitIdx=[-1,-1]
	for ss in xrange(0,2):
		xx=x[ss]
		yy=y[ss]
		OBJmin=9999
		for kk in xrange(0,n-2):
			yyk=yy>>kk
			xxk0=xx>>(n-kk)
			xxk=xx-(xxk0<<(n-kk))
			tmp=xxk^yyk
			ee=PopC(tmp)
			obj=ee*3+kk		
			if obj<=OBJmin:
				OBJmin=obj
				err_min[ss]=ee
				suffL[ss]=n-kk
				bitIdx[ss]=tmp
	if bitIdx[0]==bitIdx[1] and bitIdx[0]>0:
		err_t=err_min[0]
	else:
		err_t=err_min[0]+err_min[1]
	if err_t<= err and suffL[0]==suffL[1]:
		return (suffL[0],err_t)
	else:
		return (0,0)

def PairwiseReadComparison(x,y,n=50,err=1):
	## Compare a pair of reads for their overlap in 4 scenarios: x,y; y,x; xc,y; y,xc
	suffList=[]
	suf=CompareSuffixByBit(x,y,n=n,err=err)
	suffList.append(suf)
	suf=CompareSuffixByBit(y,x,n=n,err=err)	
	suffList.append(suf)
	xc=ReverseComp(x,n=n)
	suf=CompareSuffixByBit(xc,y,n=n,err=err)
	suffList.append(suf)
	suf=CompareSuffixByBit(y,xc,n=n,err=err)	
	suffList.append(suf)
	return suffList

def BreakSeqIntoAAcodes(seq,frame=0):
        ## Break a sequence into consecutive 3 letters, with arbitrarily selected frameshift bases
        n=len(seq)
        AAseqs=[]
        i=frame
        while i<n:
            if i+3>=n:
                break
            AAseqs.append(seq[i:(i+3)])
            i+=3
        return AAseqs

def TranslateAA(seq):
        AAList=[]
        for ff in [0,1,2]:
            AAseqs=BreakSeqIntoAAcodes(seq,frame=ff)
            AA=''
            for i in AAseqs:
                if len(i)!=3:
                    continue
                AA+=AAcode[i]
            AAList.append((AA,ff+1))
        for ff in [0,1,2]:
            AAseqsRC=BreakSeqIntoAAcodes(ReverseCompSeq(seq),frame=ff)
            AA=''
            for i in AAseqsRC:
                if len(i)!=3:
                    continue
                AA+=AAcode[i]
            AAList.append((AA,-ff-1))        
        return AAList

def DetectCDR3(seq,VFaDict,JFaDict,Bcell=False):
        AAList=TranslateAA(seq)
        if Bcell:
            CDR3patV0=CDR3patVB
            CDR3patJ0=CDR3patJB
        else:
            CDR3patV0=CDR3patV
            CDR3patJ0=CDR3patJ
        CDR3seq=[]
        for AA in AAList:
            AAseq=AA[0]
            if "*" in AAseq:
                continue
            Vflag=0
            Vpos=[]
            for ppk in CDR3patV0:
                pp=CDR3patV0[ppk]
                tmp=re.findall(pp,AAseq)
                if len(tmp)>0:
                    Vpos.append((ppk,AAseq.find(tmp[0])))
                    if Vflag==0:
                        Vflag=1
                    break
            Jflag=0
            Jpos=[]
            for ppk in CDR3patJ0:
                pp=CDR3patJ0[ppk]
                tmp=re.findall(pp,AAseq)
                if len(tmp)>0:
                    Jpos.append((ppk,AAseq.find(tmp[0])))
                    if Jflag==0:
                        Jflag=1
                    break
            if Vflag==0 and Jflag==0:
                CDR3=[AA,('',''),-1,('',''),-1,0,0]
            if Jflag==1 and Vflag==0:
                ## to save time, first check if J gene is correct
                bestMatchList=[]
                bestScoreList=[]
                for jj in Jpos:
                    ss=AAseq[max(0,jj[1]-4):]
                    Nj=len(ss)
                    NN=Nj
                    reverse=False
                    if Nj>10:
                        reverse=True
                        NN=10
                    bestMatch=''
                    bestScore=0
                    if Nj<6:
                        bestMatchList.append((bestMatch,jj[0]))
                        bestScoreList.append(bestScore)
                        continue
                    for Jgene in JFaDict:
                        Jseq=JFaDict[Jgene]
                        if reverse:
                            aln=pairwise2.align.localms(ss,Jseq,1,-1,-9,0)
                            if len(aln)==0:
                                continue
                            else:
                                aln=aln[0]
                        else:
                            aln=pairwise2.align.localms(Jseq,ss,1,-1,-9,0)
                            if len(aln)==0:
                                continue
                            else:
                                aln=aln[0]
                        if aln[2]>bestScore:
                            bestMatch=Jgene
                            bestScore=aln[2]
                    if bestScore<NN-4:
                        # allow for 2 mismatches or single gaps
                        bestMatch=''
                        bestScore=0
                    bestMatchList.append((bestMatch,jj[0]))
                    bestScoreList.append(bestScore)
                bSj=max(bestScoreList)
                if bSj==0:
                        CDR3 = [AA,('',''),-1,('',''),Jpos[0][1],0,0]                       
                else:
                    ppj=bestScoreList.index(bSj)
                    matched_Jgene=bestMatchList[ppj]
                    matched_Jpos=Jpos[ppj]
                    Nj=len(AAseq[matched_Jpos[1]:])
                    CDR3 = [AA,('',''),-1,matched_Jgene,matched_Jpos[1],bSj,Nj]
            if Vflag==1 and Jflag==1:
                bestMatchListJ=[]
                bestScoreListJ=[]
                for jj in Jpos:
                    ss=AAseq[max(0,jj[1]-4):]
                    Nj=len(ss)
                    NN=Nj
                    reverse=False
                    if Nj>10:
                        reverse=True
                        NN=10
                    bestMatch=''
                    bestScore=0
                    if Nj<6:
                        bestMatchListJ.append((bestMatch,jj[0]))
                        bestScoreListJ.append(bestScore)
                        continue
                    for Jgene in JFaDict:
                        Jseq=JFaDict[Jgene]
                        if reverse:
                            aln=pairwise2.align.localms(ss,Jseq,1,-1,-9,0)
                            if len(aln)==0:
                                continue
                            else:
                                aln=aln[0]
                        else:
                            aln=pairwise2.align.localms(Jseq,ss,1,-1,-9,0)
                            if len(aln)==0:
                                continue
                            else:
                                aln=aln[0]                            
                        if aln[2]>bestScore:
                            bestMatch=Jgene
                            bestScore=aln[2]
                    if bestScore<NN-4:
                        # allow for 2 mismatches or single gaps
                        bestMatch=''
                        bestScore=0
                    bestMatchListJ.append((bestMatch,jj[0]))
                    bestScoreListJ.append(bestScore)
                bSj=max(bestScoreListJ)
                bestMatchListV=[]
                bestScoreListV=[]
                for vv in Vpos:
                            ss=AAseq[0:(vv[1]+3)]
                            Nv=len(ss)
                            bestMatch=''
                            bestScore=0
                            if Nv<4:
                                bestMatchListV.append((bestMatch,vv[0]))
                                bestScoreListV.append(bestScore)
                                continue
                            matched_vgenes=[]
                            vscores=[]
                            for Vgene in VFaDict:
                                Vseq=VFaDict[Vgene][-(Nv+10):]
                                aln=pairwise2.align.localms(Vseq,ss,1,-1,-9,0)
                                if len(aln)==0:
                                    continue
                                else:
                                    aln=aln[0]
                                if aln[2]>bestScore:
                                    bestMatch=Vgene
                                    bestScore=aln[2]
                                matched_vgenes.append(Vgene)
                                vscores.append(aln[2])
                            if bestScore<Nv-10:
                                # allow for 5 mismatches or single gaps
                                bestMatch=''
                                bestScore=0
                            bestIdx=np.where(np.array(vscores)==bestScore)[0]
                            bestMatch='_'.join(list(np.array(matched_vgenes)[bestIdx]))
                            bestMatchListV.append((bestMatch,vv[0]))
                            bestScoreListV.append(bestScore)                        
                bSv=max(bestScoreListV)
                if bSv==0 and bSj==0:
                        CDR3=[AA,('',''),Vpos[0][1],('',''),Jpos[0][1],0,0]
                if bSv>0 and bSj==0:
                        ppv=bestScoreListV.index(bSv)
                        matched_Vgene=bestMatchListV[ppv]
                        matched_Vpos=Vpos[ppv]
                        Nv=len(AAseq[0:matched_Vpos[1]])
                        CDR3=[AA,matched_Vgene,matched_Vpos[1],('',''),Jpos[0][1],bSv,Nv]
                if bSv==0 and bSj>0:
                        ppj=bestScoreListJ.index(bSj)
                        matched_Jgene=bestMatchListJ[ppj]
                        matched_Jpos=Jpos[ppj]
                        Nj=len(AAseq[matched_Jpos[1]:])
                        CDR3 = [AA,('',''),Vpos[0][1],matched_Jgene,matched_Jpos[1],bSj,Nj]
                if bSv>0 and bSj>0:
                        ppj=bestScoreListJ.index(bSj)
                        matched_Jgene=bestMatchListJ[ppj]
                        matched_Jpos=Jpos[ppj]
                        Nj=len(AAseq[matched_Jpos[1]:])
                        ppv=bestScoreListV.index(bSv)
                        matched_Vgene=bestMatchListV[ppv]
                        matched_Vpos=Vpos[ppv]
                        Nv=len(AAseq[0:matched_Vpos[1]])
                        CDR3 = [AA,matched_Vgene,matched_Vpos[1],matched_Jgene,matched_Jpos[1],bSv+bSj,Nv+Nj]
            if Vflag==1 and Jflag==0:
                ## Vflag=1,Jflag=0
                bestMatchList=[]
                bestScoreList=[]
                for vv in Vpos:
                    ss=AAseq[0:(vv[1]+3)]
                    Nv=len(ss)
                    bestMatch=''
                    bestScore=0
                    if Nv<4:
                        bestMatchList.append((bestMatch,vv[0]))
                        bestScoreList.append(bestScore)
                        continue
                    matched_vgenes=[]
                    vscores=[]
                    for Vgene in VFaDict:
                        Vseq=VFaDict[Vgene][-(Nv+10):]
                        aln=pairwise2.align.localms(Vseq,ss,1,-1,-9,0)
                        if len(aln)==0:
                            continue
                        else:
                            aln=aln[0]
                        if aln[2]>bestScore:
                            bestMatch=Vgene
                            bestScore=aln[2]
                        matched_vgenes.append(Vgene)
                        vscores.append(aln[2])
                    if bestScore<Nv-10:
                        # allow for 5 mismatches or single gaps
                        bestMatch=''
                        bestScore=0
                    bestIdx=np.where(np.array(vscores)==bestScore)[0]
                    bestMatch='_'.join(list(np.array(matched_vgenes)[bestIdx]))
                    bestMatchList.append((bestMatch,vv[0]))
                    bestScoreList.append(bestScore)                        
                bSv=max(bestScoreList)
                if bSv==0:
                    CDR3 = [AA,('',''),Vpos[0][1],('',''),-1,0,0]
                else:
                    ppv=bestScoreList.index(bSv)
                    matched_Vgene=bestMatchList[ppv]
                    matched_Vpos=Vpos[ppv]
                    Nv=len(AAseq[0:matched_Vpos[1]])
                    CDR3 = [AA,matched_Vgene,matched_Vpos[1],('',''),-1,bSv,Nv]
            CDR3seq.append(CDR3)
        newCDR3seq=[]
        for CDR3 in CDR3seq:
            if CDR3[2]==-1 and CDR3[4]==-1:
                continue
            newCDR3seq.append(CDR3)
        return newCDR3seq

def AllocateReadsIntoGenes(rr,geneLoci,geneName,REFs):
	readDict={}
	for read in rr:
		if 'N' in read.seq and read.flag & 4>0: ## Remove unmapped reads with "N"
			continue
		if read.qname not in readDict:
                    if '/1' == read.qname[-2:] or '/2' == read.qname[-2:]:
                        ## explicit mate pair information
			if '/1' == read.qname[-2:]:
				sq=re.sub('/1','/2',read.qname)
			else:
				sq=re.sub('/2','/1',read.qname)
                        if sq in readDict:
				readDict[sq].append(read)
                        else:
				readDict[sq]=[read]
		    elif '.1' == read.qname[-2:] or '.2' == read.qname[-2:]:
                        ## explicit mate pair information
			if '.1' == read.qname[-2:]:
                                sq=read.qname[0:-2]+'.2'
				#sq=re.sub('.1','.2',read.qname)
			else:
				#sq=re.sub('.2','.1',read.qname)
                                sq=read.qname[0:-2]+'.1'
                        if sq in readDict:
				readDict[sq].append(read)
                        else:
				readDict[sq]=[read]                        
                    else:
                        readDict[read.qname]=[read]
		else:
			readDict[read.qname].append(read)
	PairedReadDict={}
	for kk in readDict:
		vv=readDict[kk]
		if len(vv)==1:
			continue
		if len(vv)>2:
                    ## multiple hits
                    continue
		if vv[0].flag & 4 > 0 and vv[1].flag & 4 == 0:
			vv=(1,vv)
		elif vv[1].flag & 4 > 0 and vv[0].flag & 4 == 0:
			vv=(2,vv)
		else:
			continue
		PairedReadDict[kk]=vv
	geneDict={}
	CHRlist=[]
	stList=[]
	edList=[]
	for gene in geneLoci:		
		tmp=gene.split(':')
		CHR=tmp[0]
		tmp=tmp[1].split('-')
		st=int(tmp[0])
		ed=int(tmp[1])
		if 'chr7' in REFs:
			CHRlist.append(CHR)
		else:
			CHRlist.append(re.sub('chr','',CHR))
		stList.append(st)
		edList.append(ed)
	ng=len(CHRlist)
	CHRlist=np.array(CHRlist)
	stList=np.array(stList)
	edList=np.array(edList)
	for kk in PairedReadDict:
		vv=PairedReadDict[kk]
		if vv[0]==1:
			CHR0=REFs[vv[1][1].rname]
			pos=vv[1][1].pos
		if vv[0]==2:
			CHR0=REFs[vv[1][0].rname]
			pos=vv[1][0].pos
		ssCHR=set(list(np.where(CHRlist==CHR0)[0]))
		ss_st=set(list(np.where(stList<=pos)[0]))
		ss_ed=set(list(np.where(edList>=pos)[0]))
		ss=list(ssCHR & ss_st & ss_ed)
		if len(ss)==0:
			continue
		for s in ss:
			gg=geneName[s]+'|'+str(CHRlist[s])+':'+str(stList[s])+'-'+str(edList[s])
			if gg not in geneDict:
				geneDict[gg]=[vv]
			else:
				geneDict[gg].append(vv)
	return readDict,PairedReadDict,geneDict	

def OrderUnmappedReads(geneObj,nr):
	## Preliminary ordering of unmapped reads using their mapped pair
	urList=[]
	for vv in geneObj:
		if vv[0]==1:
			urList.append((vv[1][1].pos,vv[1][0]))
		if vv[0]==2:
			urList.append((vv[1][0].pos,vv[1][1]))
	urList_sorted=sorted(urList,key=lambda x:x[0])
	pos_min=urList_sorted[0][0]
	pos_max=urList_sorted[-1][0]
	SLICE=(pos_max-pos_min)/nr
	print "-----dividing unmapped reads into %d slices" %(SLICE)
	UR_sliced=[]
	cur_pos=pos_min
	tmpList=[]
	for kk in urList_sorted:
		if kk[0]<= cur_pos+nr:
			tmpList.append(kk[1])
		else:
			UR_sliced.append(tmpList)
			cur_pos=kk[0]	
			tmpList=[kk[1]]
	return UR_sliced

def GetSeqOverlap(ContigObj,ContigNames,overlap_thr=10,err=1):
        ## compare contigs of assembled sequences, currently do not tolerate errors
        OverlapInfo={}
        nc=len(ContigObj)
        temp=[]
        for i in xrange(0,nc):
            x=ContigObj[i]
            for j in xrange(i,nc):
                if j==i:
                    continue
                ## if both contigs belong to V/D/J/C category but different genes, overlap does not result in merging
                gg1=ContigNames[i].split('|')
                gg2=ContigNames[j].split('|')
                if gg1[0][0:3] != gg2[0][0:3]:
                        ## different TCR gene, very unlikely to happen
                    continue
                if (gg1[0][0:4] == gg2[0][0:4] and gg1[0] != gg2[0]) or (gg1[0]==gg2[0] and gg1[1]!=gg2[1]):
                        ## same TCR gene, same VDJC category, different genes
                    continue
                y=ContigObj[j]
                OP=CompareSuffixSeq(x,y,err=err)
                if OP[0][0]>=overlap_thr:
                    kk=ContigNames[i]+'\t'+ContigNames[j]
                    if ContigNames[i] not in temp:
                        temp.append(ContigNames[i])
                    if ContigNames[j] not in temp:
                        temp.append(ContigNames[j])
                    OverlapInfo[kk]=OP
        for ii in xrange(0,len(ContigNames)):
            if ContigNames[ii] not in temp:
                OverlapInfo[ContigNames[ii]+'\t'+ContigNames[ii]]=((len(ContigObj[ii]),0),0)
        return OverlapInfo

def GetReadsOverlapByGene(geneObj,nr=50,err=1,overlap_thr=10):
	UnmappedReads=[]
	temp=[]
	for vv in geneObj:
		if vv[0]==3:
			continue
		if vv[0]==1:
			UnmappedReads.append(vv[1][0])
		if vv[0]==2:
			UnmappedReads.append(vv[1][1])
	if len(UnmappedReads)>=200:
		ns=len(UnmappedReads)/200
		if ns>1:
			nr0=nr/ns+1
		else:
			nr0=nr
		UR_sliced=OrderUnmappedReads(geneObj,nr0)
	else:
		UR_sliced=[UnmappedReads]	
	OverlapInfo={}
	for UR in UR_sliced:
		print len(UR)
		BitReads=[]
		qNames=[]
		for read in UR:
			if len(read.seq) != nr:
				continue
			qNames.append(read.qname)
			BitReads.append(ConvertDNAtoBinary(read.seq))
		Nu=len(BitReads)
		for i in xrange(0,Nu):
			x=BitReads[i]
			for j in xrange(i,Nu):
				if j==i:
					continue
				y=BitReads[j]
				sL=PairwiseReadComparison(x,y,n=nr,err=err)
				tmpMax=-1
				tmpIdx=-1
				for k in xrange(0,4):
					if sL[k][0]>tmpMax:
						tmpMax=sL[k][0]
						tmpIdx=k
				if tmpMax>=overlap_thr+sL[tmpIdx][1]*3:
					qnn=qNames[i]+'\t'+qNames[j]
					if qNames[i] not in temp:
                                                temp.append(qNames[i])
					if qNames[j] not in temp:
                                                temp.append(qNames[j])
					OverlapInfo[qnn]=(sL[tmpIdx],tmpIdx)
        for rr in UnmappedReads:
            if rr.qname not in temp:
                OverlapInfo[rr.qname+'\t'+rr.qname]=((nr,0),0)
        return OverlapInfo,UnmappedReads

def GetReadsOverlapByGene_SE(geneObj,nr=75,err=1,overlap_thr=10):
        UnmappedReads=[]
        temp=[]
        for rr in geneObj:
            UnmappedReads.append(geneObj[rr][0][0])
	OverlapInfo={}
        BitReads=[]
        qNames=[]
        for read in UnmappedReads:
                if len(read.seq) != nr:
                        continue
                qNames.append(read.qname)
                BitReads.append(ConvertDNAtoBinary(read.seq))
        Nu=len(BitReads)
        for i in xrange(0,Nu):
                x=BitReads[i]
                for j in xrange(i,Nu):
                        if j==i:
                                continue
                        y=BitReads[j]
                        sL=PairwiseReadComparison(x,y,n=nr,err=err)
                        tmpMax=-1
                        tmpIdx=-1
                        for k in xrange(0,4):
                                if sL[k][0]>tmpMax:
                                        tmpMax=sL[k][0]
                                        tmpIdx=k
                        if tmpMax>=overlap_thr+sL[tmpIdx][1]*3:
                                qnn=qNames[i]+'\t'+qNames[j]
                                if qNames[i] not in temp:
                                        temp.append(qNames[i])
                                if qNames[j] not in temp:
                                        temp.append(qNames[j])
                                OverlapInfo[qnn]=(sL[tmpIdx],tmpIdx)
        for rr in UnmappedReads:
            if rr.qname not in temp:
                OverlapInfo[rr.qname+'\t'+rr.qname]=((nr,0),0)
        return OverlapInfo

def FindDisjointCommunities(OverlapInfo):
	kk=OverlapInfo.keys()
	OpDict={}
	for k in kk:
		k=k.split('\t')
		if k[0] in OpDict:
			OpDict[k[0]].append(k[1])
		else:
			OpDict[k[0]]=[k[1]]
		if k[1] in OpDict:
			OpDict[k[1]].append(k[0])
		else:
			OpDict[k[1]]=[k[0]]
	uniqReads=OpDict.keys()
	Nu=len(uniqReads)
	DisComm=[]
	tmpL=list(chain(*DisComm))
	def LoadComm(STACK,cur_read):
		if cur_read in STACK:
			return
		else:
			STACK.append(cur_read)
			vv=OpDict[cur_read]
			for v in vv:
				LoadComm(STACK,v)
		return STACK
	for read in uniqReads:
		if read in tmpL:
			continue
		else:
			STACK=LoadComm([],read)
			DisComm.append(STACK)
			tmpL=list(chain(*DisComm))
	return DisComm,OpDict

def MergeCommReads(assembleReads,assembleSeqs,OverlapInfo,nr=50,thr_overlap=10):
        assembleReadsNew=deepcopy(assembleReads)
        assembleSeqsNew=deepcopy(assembleSeqs)
        while True:
            m=len(assembleReadsNew)
            if m<=1:
                break
            flag=0
            for ii in xrange(0,m):
                flag_innerloop=0
                for jj in xrange(ii,m):
                    if jj==ii:
                        continue
                    ss1=assembleReadsNew[ii]
                    ss2=assembleReadsNew[jj]
                    sq1=assembleSeqsNew[ii]
                    sq2=assembleSeqsNew[jj]
                    op=list(set(ss1)&set(ss2))
                    if len(op)==0:
                        continue
                    ## find overlap, calculate the start positions for each read in each set
                    rr=op[0]
                    idx1=ss1.index(rr)
                    idx2=ss2.index(rr)
                    sign1=sq1[idx1][1]
                    sign2=sq2[idx2][1]
                    ss2o=ss2
                    sq2o=sq2
                    if sign1!=sign2:
                        ss2=ss2[::-1]
                        sq2_new=[]
                        for tmp in sq2:
                            sq2_new=[(ReverseCompSeq(tmp[0]),-tmp[1])]+sq2_new
                        sq2=sq2_new
                    if len(op)>1:
                        ## find if the reads are forming a loop
                        flag_loop=0
                        for rr1 in op[1:]:
                            if (ss1.index(rr1)-idx1)*(ss2.index(rr1)-idx2)<0:
                                ## find a loop, no merging
                                flag_loop=1
                                break
                        if flag_loop==1:
                            continue
                    n1=len(ss1)
                    n2=len(ss2)
                    idx1=ss1.index(rr)
                    idx2=ss2.index(rr)                   
                    ss1_pos=[0]
                    ss2_pos=[0]
                    for kk in xrange(idx1,n1):
                        if kk==idx1:
                            continue
                        tmpKey=ss1[kk-1]+'\t'+ss1[kk]
                        if tmpKey not in OverlapInfo:
                            tmpKey=ss1[kk]+'\t'+ss1[kk-1]
                        tmpV=OverlapInfo[tmpKey]
                        ss1_pos+=[ss1_pos[-1]+nr-tmpV[0][0]]
                    for kk in xrange(idx1,-1,-1):
                        if kk==idx1:
                            continue
                        tmpKey=ss1[kk+1]+'\t'+ss1[kk]
                        if tmpKey not in OverlapInfo:
                            tmpKey=ss1[kk]+'\t'+ss1[kk+1]
                        tmpV=OverlapInfo[tmpKey]
                        ss1_pos=[ss1_pos[0]-nr+tmpV[0][0]]+ss1_pos
                    for kk in xrange(idx2,n2):
                        if kk==idx2:
                            continue
                        tmpKey=ss2[kk-1]+'\t'+ss2[kk]
                        if tmpKey not in OverlapInfo:
                            tmpKey=ss2[kk]+'\t'+ss2[kk-1]
                        tmpV=OverlapInfo[tmpKey]
                        ss2_pos+=[ss2_pos[-1]+nr-tmpV[0][0]]
                    for kk in xrange(idx2,-1,-1):
                        if kk==idx2:
                            continue
                        tmpKey=ss2[kk+1]+'\t'+ss2[kk]
                        if tmpKey not in OverlapInfo:
                            tmpKey=ss2[kk]+'\t'+ss2[kk+1]
                        tmpV=OverlapInfo[tmpKey]
                        ss2_pos=[ss2_pos[0]-nr+tmpV[0][0]]+ss2_pos
                    ## try to merge two sets
                    flag_NoOverlap=0
                    i=idx1
                    j=idx2
                    while i<n1 and j<n2:
                        if ss1[i]==ss2[j]:
                            i+=1
                            j+=1
                            continue
                        if ss1_pos[i]<=ss2_pos[j]:
                            tmpKey=ss1[i]+'\t'+ss2[j]
                            if tmpKey not in OverlapInfo:
                                tmpKey=ss2[j]+'\t'+ss1[i]
                            if tmpKey not in OverlapInfo:
                                flag_NoOverlap=1
                                break
                            i+=1
                        else:
                            tmpKey=ss1[i]+'\t'+ss2[j]
                            if tmpKey not in OverlapInfo:
                                tmpKey=ss2[j]+'\t'+ss1[i]
                            if tmpKey not in OverlapInfo:
                                flag_NoOverlap=1
                                break
                            j+=1
                    i=idx1
                    j=idx2
                    while i>=0 and j>=0:
                        if ss1[i]==ss2[j]:
                            i-=1
                            j-=1
                            continue
                        if ss1_pos[i]<=ss2_pos[j]:
                            tmpKey=ss1[i]+'\t'+ss2[j]
                            if tmpKey not in OverlapInfo:
                                tmpKey=ss2[j]+'\t'+ss1[i]
                            if tmpKey not in OverlapInfo:
                                flag_NoOverlap=1
                                break
                            j-=1
                        else:
                            tmpKey=ss1[i]+'\t'+ss2[j]
                            if tmpKey not in OverlapInfo:
                                tmpKey=ss2[j]+'\t'+ss1[i]
                            if tmpKey not in OverlapInfo:
                                flag_NoOverlap=1
                                break
                            i-=1                                        
                    if flag_NoOverlap==0:   ## Merging happens
                        flag=1
                        print ii,jj,'merging'
                        newComm=[ss1[idx1]]
                        newSeqs=[sq1[idx1]]
                        i=idx1
                        j=idx2
                        while i<n1 and j<n2:
                            if ss1_pos[i]<=ss2_pos[j]:
                                if ss1[i] not in newComm:
                                    newComm+=[ss1[i]]
                                    newSeqs+=[sq1[i]]
                                else:
                                    i+=1
                                    if i>=n1:
                                        for j0 in xrange(j,n2):
                                            if ss2[j0] not in newComm:
                                                newComm+=[ss2[j0]]
                                                newSeqs+=[sq2[j0]]
                                        j=n2
                            else:
                                if ss2[j] not in newComm:
                                    newComm+=[ss2[j]]
                                    newSeqs+=[sq2[j]]
                                else:
                                    j+=1
                                    if j>=n2:
                                        for i0 in xrange(i,n1):
                                            if ss1[i0] not in newComm:
                                                newComm+=[ss1[i0]]
                                                newSeqs+=[sq1[i0]]
                                        i=n1
                        i=idx1
                        j=idx2
                        while i>=0 and j>=0:
                            if ss1_pos[i]<=ss2_pos[j]:
                                if ss2[j] not in newComm:
                                    newComm=[ss2[j]]+newComm
                                    newSeqs=[sq2[j]]+newSeqs
                                else:
                                    j-=1
                                    if j==-1:
                                        for i0 in xrange(i,-1,-1):
                                            if ss1[i0] not in newComm:
                                                newComm=[ss1[i0]]+newComm
                                                newSeqs=[sq1[i0]]+newSeqs
                                        i=-1
                            else:
                                if ss1[i] not in newComm:
                                    newComm=[ss1[i]]+newComm
                                    newSeqs=[sq1[i]]+newSeqs
                                else:
                                    i-=1
                                    if i==-1:
                                        for j0 in xrange(j,-1,-1):
                                            if ss2[j0] not in newComm:
                                                newComm=[ss2[j0]]+newComm
                                                newSeqs=[sq2[j0]]+newSeqs
                                        j=-1
                        assembleReadsNew.remove(ss1)
                        assembleReadsNew.remove(ss2o)
                        assembleReadsNew.append(newComm)
                        assembleSeqsNew.remove(sq1)
                        assembleSeqsNew.remove(sq2o)
                        assembleSeqsNew.append(newSeqs)
                        flag_innerloop=1
                        break
                if flag_innerloop==1:
                    break
            if flag==0:
                break
        return assembleReadsNew,assembleSeqsNew                           

def EMcount(contigReads,countDict,max_iter=1000,thr_s=0.001):
        nc=len(contigReads)
        uniqueReads=[]
        sharedReads=[]
        for rr in countDict:
            if countDict[rr]>1:
                sharedReads.append(rr)
            else:
                uniqueReads.append(rr)
        uniqueCounts=[]
        sharedCounts=[]
        for cc in contigReads:
            temp=list(set(uniqueReads) & set(cc))
            uniqueCounts.append(len(temp))
            temp=list(set(sharedReads) & set(cc))
            sharedCounts.append(temp)
        rsDict={}
        for read in sharedReads:
            for ii in xrange(0,nc):
                if read in contigReads[ii]:
                    if read not in rsDict:
                        rsDict[read]=[ii]
                    else:
                        rsDict[read].append(ii)
        contigReadCount0=[]
        for cc in contigReads:
            contigReadCount0.append(len(cc))
        ii=0
        flag=0
        while ii<max_iter:
            contigReadCount=[]
            for jj in xrange(0,nc):
                nu=uniqueCounts[jj]
                ss=sharedCounts[jj]
                ns=0
                for rs in ss:
                    shared_contig=rsDict[rs]
                    Ns=0
                    for sc in shared_contig:
                        Ns+=contigReadCount0[sc]
                    ns+=float(contigReadCount0[jj])/Ns
                contigReadCount.append(nu+ns)
            max_diff=-1
            for jj in xrange(0,nc):
                diff=abs(contigReadCount0[jj]-contigReadCount[jj])
                if diff > max_diff:
                    max_diff=diff
            if max_diff<=thr_s:
                flag=1
                break
            contigReadCount0=contigReadCount
            ii+=1
        if flag==0:
            print "No convergence"
        return contigReadCount

def AssembleCommReads(readComm,OpDict,OverlapInfo,PairedReadDict,readDict,contigDict={},nr=50,thr_overlap=10,mode="read"):
        if mode not in ['read','contig']:
            print 'wrong mode input!'
            raise
	kks=OverlapInfo.keys()
	vDict={}
	for qn in list(set(readComm)):
		vDict[qn]=1
	assembleReads=[]
	assembleSeqs=[]
	while vDict.values().count(1)>0:
		for qn in vDict.keys():
			if vDict[qn]==1:
				break
		if mode=='read':
#                    if qn not in tmp:
#                        if '\1' in qn or '\2' in qn:
#                            if '\1' in qn:
#                                qn=re.sub('\1','\2',qn)
#                            else:
#                                qn=re.sub('\2','\1',qn)
#                        if '.1' == qn[-2:] or '.2'==qn[-2:]:
#                            if '.1'==qn[-2:]:
#                                qn[-2:]='.2'
#                            else:
#                                qn[-2:]='.1'
                    tmp=PairedReadDict[qn]
                    if tmp[0]==1:
                            tmpSeq=tmp[1][0].seq
                    if tmp[0]==2:
                            tmpSeq=tmp[1][1].seq
                    if tmp[0]==3:
                            vDict[qn]=0
                            continue
                    if len(tmpSeq) is not nr:
                        continue
                    tmpSeq0=ConvertBitToDNA(ReverseComp(ConvertDNAtoBinary(tmpSeq),n=nr),n=nr)
                if mode=='contig':
                    tmpSeq=contigDict[qn]
                    tmpSeq0=ReverseCompSeq(tmpSeq)
		if len(assembleReads)==0:
			assembleReads.append([qn])
			assembleSeqs.append([(tmpSeq,1)])
			vDict[qn]=0
		else:
			flag=0
			vDict[qn]=0
			na=len(assembleReads)
			for aa in xrange(0,na):
				tmpContig=assembleReads[aa]
				tmpSeqs=assembleSeqs[aa]
				ns=len(tmpContig)
				ns0=ns
				for ii in xrange(0,ns):
					pp=tmpContig[ii]
					pp_seq=tmpSeqs[ii]
					if pp in OpDict[qn]:
						#print OpDict[qn]
						flag=1
						tmp=pp+'\t'+qn
						direct=0
						if tmp in OverlapInfo:
							vv=OverlapInfo[tmp]
							direct=1
						else:
							vv=OverlapInfo[qn+'\t'+pp]
						#print direct,vv[1],tmpSeqs
						if direct==0:
							## x=qn,y=pp
							if vv[1]==0:
								if pp_seq[1]==1:
									if ii==0:
										tmpContig=[qn]+tmpContig
										ns0=ns+1
										tmpSeqs=[(tmpSeq,1)]+tmpSeqs
								if pp_seq[1]==-1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq0,-1)]
							if vv[1]==1:
								if pp_seq[1]==1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq,1)]
								if pp_seq[1]==-1:
									if ii==0:
										ns0=ns+1
										tmpContig=[qn]+tmpContig
										tmpSeqs=[(tmpSeq0,-1)]+tmpSeqs
							if vv[1]==2:
								if pp_seq[1]==1:
									if ii==0:
										ns0=ns+1
										tmpContig=[qn]+tmpContig
										tmpSeqs=[(tmpSeq0,-1)]+tmpSeqs
								if pp_seq[1]==-1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq,1)]
							if vv[1]==3:
								if pp_seq[1]==1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq0,-1)]
								if pp_seq[1]==-1:
									if ii==0:
										ns0=ns+1
										tmpContig=[qn]+tmpContig
										tmpSeqs=[(tmpSeq,1)]+tmpSeqs
						if direct==1:
							#  x=pp, y=qn						
							if vv[1]==1:
								if pp_seq[1]==1:
									if ii==0:
										tmpContig=[qn]+tmpContig
										ns0=ns+1
										tmpSeqs=[(tmpSeq,1)]+tmpSeqs
								if pp_seq[1]==-1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq0,-1)]
							if vv[1]==0:
								if pp_seq[1]==1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq,1)]
								if pp_seq[1]==-1:
									if ii==0:
										ns0=ns+1
										tmpContig=[qn]+tmpContig
										tmpSeqs=[(tmpSeq0,-1)]+tmpSeqs
							if vv[1]==3:
								if pp_seq[1]==-1:
									if ii==0:
										ns0=ns+1
										tmpContig=[qn]+tmpContig
										tmpSeqs=[(tmpSeq,1)]+tmpSeqs
								if pp_seq[1]==1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq0,-1)]
							if vv[1]==2:
								if pp_seq[1]==-1:
									if ii==ns0-1:
										tmpContig+=[qn]
										tmpSeqs+=[(tmpSeq,1)]
								if pp_seq[1]==1:
									if ii==0:
										ns0=ns+1
										tmpContig=[qn]+tmpContig
										tmpSeqs=[(tmpSeq0,-1)]+tmpSeqs
				assembleReads[aa]=tmpContig
				assembleSeqs[aa]=tmpSeqs
			if flag==0:
				assembleReads.append([qn])
				assembleSeqs.append([(tmpSeq,1)])
	try:
            assembleReadsNew,assembleSeqsNew=MergeCommReads(assembleReads,assembleSeqs,OverlapInfo,nr=nr,thr_overlap=thr_overlap)
        except KeyError:
            print "Warning: recurring reads in different positions. Increase overlap threshold length (-l) to reduce randomness."
            assembleReadsNew=assembleReads
            assembleSeqsNew=assembleSeqs
        ## Add remaining reads here
        remainReads=[[] for x in xrange(len(assembleReadsNew))]
        for rr in readComm:
            for ii in xrange(0,len(assembleReadsNew)):
                readSet=assembleReadsNew[ii]
                if rr in readSet:
                    continue
                else:
                    for read in readSet:
                        tmpKey=rr+'\t'+read
                        if tmpKey not in OverlapInfo:
                            tmpKey=read+'\t'+rr
                        if tmpKey not in OverlapInfo:
                            continue
                        else:
                            if rr not in remainReads[ii]:
                                remainReads[ii].append(rr)
        ## Count reads in each contig                        
        temp=list(chain(*assembleReadsNew))+list(chain(*remainReads))
        countDict=defaultdict(int)
        for ww in temp:
            countDict[ww]+=1
        contigReads=[]
        for ii in xrange(len(assembleReadsNew)):
            temp=assembleReadsNew[ii]+remainReads[ii]
            temp=list(set(temp))
            contigReads.append(temp)
        contigReadCount=EMcount(contigReads,countDict)
	## stitch up different small contigs
	AssSeqs={}
	for cc in xrange(0,len(assembleReadsNew)):
		contig=assembleReadsNew[cc]
		Seqs=assembleSeqsNew[cc]
		rR=remainReads[cc]
		rC=contigReadCount[cc]
		seq=''
		nn=len(contig)
		if nn+len(rR)<1:
			continue
		i=0
		rr_contig=[]
		rR_info=[]
		for rr in rR:
                        if mode=='read':
                            tmp=PairedReadDict[rr]      ## major change. was: tmp=readDict[rr]
                            if tmp[0]==1:   ## was tmp[0].rname==-1
                                    rR_info.append(tmp[1][0])
                            elif tmp[0]==2:
                                    rR_info.append(tmp[1][1])
                        if mode=='contig':
                            rR_info.append(rr)
		while i<nn:
			p1=contig[i]
			tag=Seqs[i][1]
			if i==0:
                                if mode=='read':
                                    tmp=PairedReadDict[p1]
                                    if tmp[0]==1:
                                        if tag==1:
                                            seq=tmp[1][0].seq
                                        else:
                                            seq=ConvertBitToDNA(ReverseComp(ConvertDNAtoBinary(tmp[1][0].seq),n=nr),n=nr)
                                            rr_contig.append(tmp[1][0])
                                    else:
                                        if tag==1:
                                            seq=tmp[1][1].seq
                                        else:
                                            seq=ConvertBitToDNA(ReverseComp(ConvertDNAtoBinary(tmp[1][1].seq),n=nr),n=nr)
                                            rr_contig.append(tmp[1][1])
                                if mode=='contig':
                                    if tag==1:
                                        seq=contigDict[p1]
                                    else:
                                        seq=ReverseCompSeq(contigDict[p1])
                                    rr_contig.append(p1)
                                i+=1
                                continue
			p0=contig[i-1]
			tmp=p0+'\t'+p1
			if tmp in OverlapInfo:
				vv=OverlapInfo[tmp]
			else:
                                tmpKey=p1+'\t'+p0
                                if tmpKey not in OverlapInfo:
                                    ## same read assembled in different contigs, rare and potentially erroneous
                                    break
				vv=OverlapInfo[p1+'\t'+p0]
			if mode=='read':
                            tmp=PairedReadDict[p1]
                            if tmp[0]==1:
                                    rr_contig.append(tmp[1][0])
                            elif tmp[0]==2:
                                    rr_contig.append(tmp[1][1])
                        if mode=='contig':
                            rr_contig.append(p1)
			seq+=Seqs[i][0][vv[0][0]:]			
			i+=1
		AssSeqs[seq]=(rr_contig,rR_info,rC,contig)
	return AssSeqs

def ParseGeneLoci(geneFile,header=False):
	## input is bed file with gene name as the last column
	gg=open(geneFile)
	if header:
		gg.readline()	## get rid of header line
	geneLoci=[]
	geneName=[]
	for line in gg.readlines():
		ww=line.strip().split()
		LL=ww[0]+':'+ww[1]+'-'+ww[2]
		gene=ww[3]
		geneLoci.append(LL)
		geneName.append(gene)
	return geneLoci,geneName

def AnnotateCDR3(Seq,pRD,ContigReads=[],geneType='',error=1,overlap_thr=10,Bcell=False):
        if Bcell:
            CDR3seq_all=DetectCDR3(Seq,bHVFaDict,bHJFaDict,Bcell=True)[0:]
        else:
            if geneType=='':
                CDR3seqA=DetectCDR3(Seq,AVFaDict,AJFaDict)
                CDR3seqB=DetectCDR3(Seq,BVFaDict,BJFaDict)
                CDR3seqD=DetectCDR3(Seq,DVFaDict,DJFaDict)
                CDR3seqG=DetectCDR3(Seq,GVFaDict,GJFaDict)
                CDR3seq_all=CDR3seqA[0:]+CDR3seqB[0:]+CDR3seqD[0:]+CDR3seqG[0:]
            else:
                if geneType=='A':
                    CDR3seqs=DetectCDR3(Seq,AVFaDict,AJFaDict)
                    CDR3seq_all=CDR3seqs[0:]
                if geneType=='B':
                    CDR3seqs=DetectCDR3(Seq,BVFaDict,BJFaDict)
                    CDR3seq_all=CDR3seqs[0:]
                if geneType=='D':
                    CDR3seqs=DetectCDR3(Seq,DVFaDict,DJFaDict)
                    CDR3seq_all=CDR3seqs[0:]
                if geneType=='G':
                    CDR3seqs=DetectCDR3(Seq,GVFaDict,GJFaDict)
                    CDR3seq_all=CDR3seqs[0:]
        BestScore=-1
        BestII=-1
        for ii in xrange(0,len(CDR3seq_all)):
            tmp=CDR3seq_all[ii]
            if tmp[-2]>BestScore:
                BestScore=tmp[-2]
                BestII=ii
        if BestII==-1:
            #print CDR3seq_all
            return []
        CDR3seq=[CDR3seq_all[BestII]]
        if len(CDR3seq)==0 or len(CDR3seq)>1:
            return []
        else:
            CDR3=CDR3seq[0]
            if (CDR3[2]==-1 or CDR3[4]==-1) and CDR3[5]==0:
                return []
            tag0=CDR3[0][1]
            if tag0>0:
                offset=tag0-1
            else:
                offset=-tag0-1
            AAseq=CDR3[0][0]
            Vinfo=CDR3[1][0]
            Vgene=''
            if len(Vinfo)>0:
                VgeneList=Vinfo.split('_')
                for VgeneC in VgeneList:
                    Vgene+='_'+VgeneC.split('|')[1]
            Vgene=Vgene[1:]
            st=int(CDR3[2])
            Jinfo=CDR3[3][0]
            if len(Jinfo)>0:
                Jgene=Jinfo.split('|')[1]
            else:
                Jgene=''
            ed=int(CDR3[4])
            if st>=ed and len(CDR3[3][0])>0:
                ## erroneous assembly due to repetitive sequence or sequence error
                return []
            cdr3_seq=''
            seq_comp=ReverseCompSeq(Seq)
            if st>=0 and ed>=0:
                cdr3_seq=AAseq[st:(ed+4)]
                if tag0>0:
                    DNA_seq=Seq[(offset+st*3):(offset+ed*3)]
                else:
                    DNA_seq=seq_comp[(offset+st*3):(offset+ed*3)]
            if st>=0 and ed==-1:
                cdr3_seq=AAseq[st:]
                if tag0>0:
                    DNA_seq=Seq[(offset+st*3):]
                else:
                    DNA_seq=seq_comp[(offset+st*3):]
            if ed>=0 and st==-1:
                cdr3_seq=AAseq[:(ed+4)]
                if tag0>0:
                    DNA_seq=Seq[:(offset+ed*3)]
                else:
                    DNA_seq=seq_comp[:(offset+ed*3)]
            if len(DNA_seq)==0:
                ff='NA'
            else:
                DNA_seq_rc=ReverseCompSeq(DNA_seq)
                rr_reads=list(set(ContigReads))
                count0=0
                tmp_seq_list=[]
                if len(rr_reads)<=2:
                    count0=len(rr_reads)
                else:
                    for rr0 in rr_reads:
                        tmp_rr0=pRD[rr0]    
                        if tmp_rr0[0]==1:
                            tmp_seq=tmp_rr0[1][0].seq
                        else:
                            tmp_seq=tmp_rr0[1][1].seq
                        tmp=CompareSuffixSeq(DNA_seq,tmp_seq,err=error)
                        tmp_seq_list.append(tmp_seq)
                        if tmp[0][0]>=min(overlap_thr,len(DNA_seq)):
                            count0+=1
                        elif DNA_seq in tmp_seq or tmp_seq in DNA_seq or DNA_seq_rc in tmp_seq or tmp_seq in DNA_seq_rc:
                            count0+=1
                ff=float(count0)/len(DNA_seq)
            S=CDR3[-2]
            m=CDR3[-1]
            if m>0:
                mLogE=S-np.log(m)
            else:
                mLogE=0
            if cdr3_seq.startswith('FCA'):
                cdr3_seq=cdr3_seq[1:]
                DNA_seq=DNA_seq[3:]
        return [cdr3_seq,DNA_seq,Vgene,Jgene,ff,mLogE]

def MergeMasterContig(MasterContig,pairedReadDict,readDict,overlap_thr=10,error=1):
        ContigDict={}
        for kk in MasterContig:
            vv=MasterContig[kk]
            for cc in xrange(0,len(vv)):
                    nn=kk+'_'+str(cc)
                    ContigDict[nn]=vv[cc][0][0]
        OI_c=GetSeqOverlap(ContigDict.values(),ContigDict.keys(),2*overlap_thr,err=error)
        DC_c,OD_c=FindDisjointCommunities(OI_c)
        ContigFinalList=[]
        for dc_c in DC_c:
            if len(dc_c)>1:
                    assContig=AssembleCommReads(dc_c,OD_c,OI_c,pairedReadDict,readDict,ContigDict,mode='contig',thr_overlap=int(1.5*overlap_thr))
                    for kka in assContig:
                        ggNs=[]
                        vva=assContig[kka]
                        vv=vva[0]+vva[1]
                        Nreads=0
                        ReadsInfo=[]
                        for contig_ID in vv:
                            temp=contig_ID.split('_')
                            ggNs.append(temp[0])
                            Nreads+=MasterContig[temp[0]][int(temp[1])][0][1]
                            ReadsInfo+=MasterContig[temp[0]][int(temp[1])][2]
                        ContigFinalList.append((kka,Nreads,'__'.join(list(set(ggNs))),ReadsInfo))
            else:
                    temp=dc_c[0].split('_')
                    vv=MasterContig[temp[0]][int(temp[1])]
                    ContigFinalList.append((vv[0][0],vv[0][1],temp[0],vv[2]))                                                
        return ContigFinalList

def MergeUnmappedPairs(Contigs):
        readNames=[]
        NewContigs=[]
        for cc in Contigs:
            if len(cc[-1])>1:
                NewContigs.append(cc)
                continue
            qname=cc[-1][0]
            if '/1' in qname or '/2' in qname:
                qname=qname[0:-2]
            if qname not in readNames:
                NewContigs.append(cc)
                readNames.append(qname)
        return NewContigs

def ProcessSingleEndReads(fname,LocusFile,HeavyChain=True,err=1,overlap_thr=10,fasta=True,unmapped=False,Bcell=False):
        time1=time.time()
        time0=time1
        CDR3_rr={}
        count=1
        fHandle=pysam.Samfile(fname)
        nr=0
        N_all=0
        if not unmapped:
            gLocus=open(LocusFile)
            gLocus.readline()
            print "Estimate Read count in Locus"
            time1=time.time()
            Locs=[]
            for line in gLocus.readlines():
                Locs.append(line.strip().split('\t'))
            for Loc in Locs:
                    CHR=Loc[1]
                    st=int(Loc[2])
                    ed=int(Loc[3])
                    fHandle=pysam.Samfile(fname)
                    if 'chr7' in fHandle.references:
                            if 'chr' not in CHR:
                                CHR='chr'+CHR
                    else:
                            CHR=re.sub('chr','',CHR)
                    for read in fHandle.fetch(CHR,st,ed):
                            N_all+=1
        time2=time.time()
        print "... time elapsed %f" %(time2-time1)
        print N_all
        print "Extract reads with joining gene DNA motif"
        for rr in fHandle.fetch(until_eof=True):
            if nr==0:
                nr=rr.rlen
            if rr.flag & 4 ==0:
                continue
            if 'N' in rr.seq:
                continue
            if count % 1000000==0:
                time2=time.time()
                print count, time2-time1
            count+=1
            if Bcell:
                for geneKey in patDictB:
                    patDNA=patDictB[geneKey]
                    mm=patDNA.findall(rr.seq)
                    if len(mm)>0:
                        if geneKey not in CDR3_rr:
                            CDR3_rr[geneKey]=[rr]
                        else:
                            CDR3_rr[geneKey].append(rr)
#                        break
            else:
                if HeavyChain:
                    for geneKey in patDict:
                        patDNA=patDict[geneKey]
                        mm=patDNA.findall(rr.seq)
                        if len(mm)>0:
                            if geneKey not in CDR3_rr:
                                CDR3_rr[geneKey]=[rr]
                            else:
                                CDR3_rr[geneKey].append(rr)
#                            break
                else:
                    for geneKey in patDictAll:
                        patDNA=patDictAll[geneKey]
                        mm=patDNA.findall(rr.seq)
                        if len(mm)>0:
                            if geneKey not in CDR3_rr:
                                CDR3_rr[geneKey]=[rr]
                            else:
                                CDR3_rr[geneKey].append(rr)
#                                break
        time2=time.time()
        print "... time elapsed %f" %(time2-time1)
        print count
        print "Convert DNA sequence into Amino Acid"
        time1=time.time()
        CDR3_AA={}
        for geneKey in CDR3_rr:
            tmp_rr=CDR3_rr[geneKey]
            CDR3_AA[geneKey]=[]
            for rr in tmp_rr:
                CDR3_AA[geneKey].append((rr,TranslateAA(rr.seq)))
        time2=time.time()
        print "... time elapsed %f" %(time2-time1)
        print len(CDR3_AA)
        print "Filter in reads with gene AA motif"
        time1=time.time()
        filteredCDR3_AA={}
        if Bcell:
            MotifDict0=MotifDictB
            DNAFaDict0=DNAFaDictB
        else:
            MotifDict0=MotifDict
            DNAFaDict0=DNAFaDict
        for geneKey in CDR3_AA:
            Motifs=MotifDict0[geneKey[0:2]]
            AApats=Motifs.keys()
            AApats.sort(key=len,reverse=True)
            filteredCDR3_AA[geneKey]=[]
            vv=CDR3_AA[geneKey]
            for tmpAA in vv:
                flag=0
                for AA in tmpAA[1]:
                    if '*' in AA[0]:
                        continue
                    for pp in AApats:
                        pat=re.compile(pp)
                        tmp=pat.findall(AA[0])
                        if len(tmp)>0:
                            filteredCDR3_AA[geneKey].append((tmpAA[0],pp,AA))
                            flag=1
                            break
                    if flag==1:
                        break
        time2=time.time()
        print "... time elapsed %f" %(time2-time1)
        print len(filteredCDR3_AA)
        print "Realign kept CDR3 DNA sequences to IMGT reference genes"
        time1=time.time()
        MatchedCDR3s={}
        for geneKey in filteredCDR3_AA:
            fCDR3=filteredCDR3_AA[geneKey]
            Motifs=MotifDict0[geneKey[0:2]]
            MatchedCDR3s[geneKey]=[]
            DNAFa=DNAFaDict0[geneKey[0:2]]
            for AA in fCDR3:
                if 'reverse' in geneKey:
                    seq=ReverseCompSeq(AA[0].seq)
                else:
                    seq=AA[0].seq
                Pats=Motifs[AA[1]]
                matched_aln=[]
                MaxS=-1
                for pat in Pats:
                    gene=pat.split('|')[1]
                    for kk in DNAFa:
                        if gene in kk:
                            break
                    SEQ=DNAFa[kk].upper()
                    if len(SEQ)>100:    ## For variable genes, only use the 3' end
                        SEQ=SEQ[-45:]
                    aln=pairwise2.align.localms(seq,SEQ,1,-1,-9,0)
                    if len(aln)==0:
                        continue
                    if aln[0][2]>= 0.8 * (aln[0][4]-aln[0][3]) and aln[0][2]>=20:
                        if aln[0][2]>MaxS:
                            MaxS=aln[0][2]
                            matched_aln=[gene,MaxS]
                if len(matched_aln)>0:
                    MatchedCDR3s[geneKey].append((AA[0],AA[2],matched_aln))
        time2=time.time()
        print "... time elapsed %f" %(time2-time1)
        print len(MatchedCDR3s)
        print "Allocate unmapped SE reads into TCR genes"
        geneDict={} ## dictionary of dictionary: gene name first, then qname
        rD={}
        pRD={}
        for kk in MatchedCDR3s:
            vv=MatchedCDR3s[kk]
            for v0 in vv:
                g0=v0[2]
                qname=v0[0].qname
                if qname not in rD:
                    rD[qname]=[v0[0]]
                    pRD[qname]=(1,[v0[0]])
                else:
                    qname=qname+'_1'
                    rD[qname]=v0[0]
                    pRD[qname]=(1,[v0[0]])
                if g0[0] not in geneDict:
                    geneDict[g0[0]]={}
                    geneDict[g0[0]][qname]=[(v0[0],v0[1],v0[2])]
                else:
                    geneDict[g0[0]][qname]=[(v0[0],v0[1],v0[2])]
        print "Assemble reads allocated for each gene"
        time1=time.time()
        MasterContig={}
        for kk in geneDict:
            vv=geneDict[kk]
            OI=GetReadsOverlapByGene_SE(vv,nr=nr,err=err,overlap_thr=overlap_thr)
            DC,OD=FindDisjointCommunities(OI)
            if fasta:
                    Contigs=[]
            else:
                    Contigs=''
            for dc in DC:
                    ASs=AssembleCommReads(dc,OD,OI,pRD,rD,nr=nr,thr_overlap=overlap_thr)
                    for uu in ASs:
                            nk=ASs[uu][2]
                            if fasta:
                                    Contigs.append(([uu,nk],ASs[uu][0]+ASs[uu][1],ASs[uu][3]))
                            else:
                                    Contigs+=uu+':'+str(nk)+'_'
            MasterContig[kk]=Contigs
        time2=time.time()
        print "... time elapsed %f" %(time2-time1)
        print len(MasterContig)
        print "Merging contigs"
        time1=time.time()
        ContigFinalList=MergeMasterContig(MasterContig,pRD,rD,overlap_thr,err)
        time2=time.time()
        print "... time elapsed %f" %(time2-time1)
        if unmapped:
            ## Merge final contigs based on paired-end information
            ContigFinalList=MergeUnmappedPairs(ContigFinalList)
        time1=time.time()
        print "Annotate CDR3 contigs"
        if Bcell:
            PATV=CDR3patVB
            PATJ=CDR3patJB
        else:
            PATV=CDR3patV
            PATJ=CDR3patJ
        annList=[]
        for cc in ContigFinalList:
            rn=cc[3][0]
            for kk in geneDict:
                if rn in geneDict[kk]:
                    v0=geneDict[kk][rn]
                    break
            if False:
#            if len(cc[3])==1:
                AAseq=v0[0][1][0]
                direction=v0[0][1][1]
                if 'V' in kk:
                    for pat in PATV.values():
                        tmp=re.findall(pat,AAseq)
                        if len(tmp)>0:
                            break
                    pos=AAseq.find(tmp[0])
                    CDR3seq=AAseq[pos:]
                    if direction>0:
                        DNAseq=cc[0][pos*3+direction-1:]
                    else:
                        DNAseq=ReverseCompSeq(cc[0])[pos*3-direction-1:]
                    Vgene=kk
                    Jgene=''
                if 'J' in kk:
                    for pat in PATJ.values():
                        tmp=re.findall(pat,AAseq)
                        if len(tmp)>0:
                            break
                    pos=AAseq.find(tmp[0])+4
                    CDR3seq=AAseq[0:pos]
                    if direction>0:
                        DNAseq=cc[0][direction-1:pos*3]
                    else:
                        DNAseq=ReverseCompSeq(cc[0])[-direction-1:pos*3]
                    Vgene=''
                    Jgene=kk
                Smax=v0[0][-1][1]
                Smax=Smax/3.0
                mlogE=Smax-np.log(Smax)
                if CDR3seq.startswith('FCA'):
                    CDR3seq=CDR3seq[1:]
                    DNAseq=DNAseq[3:]
                ann=[CDR3seq,DNAseq,Vgene,Jgene,1.0/len(cc[0]),mlogE]                
            else:
                #ann=AnnotateCDR3(cc[0],pRD,cc[3],kk[2],err,overlap_thr)
                ann=AnnotateCDR3(cc[0],pRD,cc[3],'',err,overlap_thr,Bcell=Bcell)
            annList.append(ann)
        time2=time.time()
        print "... time elapse %f" %(time2-time1)
        return annList,ContigFinalList,N_all

def ScreenGenome(fname,LocusFile,InsThr=10,saveDir='./'):
	gLocus=open(LocusFile)
	gLocus.readline()   ## get rid of header line
	Locs=[]
	for line in gLocus.readlines():
            Locs.append(line.strip().split('\t'))
	seqDir={}
	CHRlist=[]
	stlist=[]
	edlist=[]
	nr=0
	for Loc in Locs:
		CHR=Loc[1]
		st=int(Loc[2])
		ed=int(Loc[3])
		CHRlist.append(CHR)
		stlist.append(st)
		edlist.append(ed)
        	handle=pysam.Samfile(fname)
		if 'chr7' in handle.references:
			if 'chr' not in CHR:
                            CHR='chr'+CHR
		else:
			CHR=re.sub('chr','',CHR)
        	for read in handle.fetch(CHR,st,ed):
                        if nr==0:
                            nr=read.rlen
                	if read.mapq<=30:
                        	continue
                	if read.qname not in seqDir:
                        	seqDir[read.qname]=[read]
                	else:
                        	seqDir[read.qname].append(read)
        print '''Pair reads in the region'''
        count=1
        count_r=0
	InsertSize=[]
        refs=handle.references
        unmapped_reads=[]
        paired_flag=0
        for rr in handle.fetch(until_eof=True):
                if rr.flag & 1 ==1 and paired_flag==0:
                    paired_flag=1
                if rr.flag & 4 >0:
                        unmapped_reads.append(rr)
                if count % 1000000 ==0:
                        print "--processed %d reads" %(count)
                        if paired_flag==0:
                            return 1
                count+=1
		if count %10000==0:
			insize=np.fabs(rr.pnext-rr.pos)+rr.rlen
			if insize<1000:
				InsertSize.append(insize)
			if len(InsertSize)>=1000:
                            medIns=np.median(np.array(InsertSize))
                            if medIns<= 2*nr - InsThr:      ## This is the case when insert size is too small to apply paired-end algorithm. Switch to single end mode automatically.
                                return 1
                if rr.qname in seqDir:
                        vv=seqDir[rr.qname]
                        flag=0
                        for v in vv:
                                if rr.seq==v.seq:
                                        flag=1
                        if flag==0:
                                seqDir[rr.qname].append(rr)
		else:
			if '/1' in rr.qname or '/2' in rr.qname:
                                if '/1' in rr.qname:
                                        qname_paired=re.sub('/1','/2',rr.qname)
                                else:
                                        qname_paired=re.sub('/2','/1',rr.qname)
			elif '.1' in rr.qname or '.2' in rr.qname:
                                if '.1' in rr.qname:
                                        qname_paired=re.sub('.1','.2',rr.qname)
                                else:
                                        qname_paired=re.sub('.2','.1',rr.qname)
                        else:
                                qname_paired=rr.qname
			if qname_paired in seqDir:	
                        	seqDir[qname_paired].append(rr)
                        	count_r+=1
                        	if count_r % 100 ==0:
                                #print refs[rr.rname], rr.pos
                                	print "---retrived %d unmapped reads" %(count_r)
        if paired_flag==0:  ## In case the total library has fewer than 1M reads
            return 1
	nfname=os.path.join(saveDir, fname.split('/')[-1]+'-Locs.bam')
	HH=handle.header
        ghandle=pysam.Samfile(nfname,mode='wb',header=HH,referencenames=handle.references,referencelengths=handle.nreferences)
        for kk in seqDir:
                for read in seqDir[kk]:
                        ghandle.write(read)
        pDict={}
        for rr in unmapped_reads:
             if rr.qname not in pDict:
                     if '/1' in rr.qname:
                             newName=re.sub('/1','/2',rr.qname)
                     else:
                             newName=re.sub('/2','/1',rr.qname)
                     if newName in pDict:
                             pDict[newName].append(rr)
                     else:
                             pDict[rr.qname]=[rr]
             else:
                     pDict[rr.qname].append(rr)
        paired_unmapped=[]
        for kk in pDict:
             if len(pDict[kk])==2:
                     paired_unmapped.append(pDict[kk])
        unmapped_file_name=os.path.join(saveDir, fname.split('/')[-1]+'-unmapped.bam')
        uhandle=pysam.Samfile(unmapped_file_name,mode='wb',header=HH,referencenames=handle.references,referencelengths=handle.nreferences)
        for kk in paired_unmapped:
            for read in kk:
                uhandle.write(read)
	return 0

def CommandLineParser():
	parser=OptionParser()
	print '''
===================================================================
Tcr Repertoire Utilities for Solid Tissue, or TRUST, is a toolbox
for analyzing T cell receptors in solid tumors using unselected
RNA-seq data. TRUST performs de novo assembly of informative unmapped
reads to estimate the TCR transcripts. It also estimates the CDR3 
sequences and the relative fractions of different T cell 
clonotypes. TRUST is developed by Bo Li (bli@jimmy.harvard.edu), 
with all rights reserved. TRUST source code shall not be distributed
without the consent of the author.
===================================================================
'''
	parser.add_option("-d","--directory",dest="Directory",help="Input bam directory",default="")
	parser.add_option("-f","--file",dest="file",default='',help="Input bam file: if given, overwite -d option")
	parser.add_option("-F","--fileList",dest="files",default='',help='Alternative input: a file containing the full path to all the files. If given, overwrite -d and -f option')
	parser.add_option("-m","--mode",dest="RunningMode",default="Full",help="Running mode. Accept Cov and Full. Cov: only report coverage information on each gene; Full: run full analysis, slow. Default: Full")
	parser.add_option("-e","--error",dest="Error",type="int",default=1,help="Maximum number of sequencing error per repeating unit tolerated. Default: 1")
	parser.add_option("-l","--overlaplength",dest="Length",type="int",default=10,help="Minimum length of overlap sequence for calling reads overlap. Default 10")
	parser.add_option('-a',"--fasta",dest="fasta",default=False,action="store_true",help='Whether or not output fasta format, only in Full mode. Default False')
        parser.add_option('-s',"--Single",dest="SingleEnd",default=False,action="store_true",help="If set True, TRUST will always run in single end mode")
        parser.add_option('-H',"--HeavyChainOnly",dest="HeavyChain",default=True,action="store_false",help="To save time, in single-end mode, TRUST only search for beta and delta chains unless this flag is set.")
        parser.add_option('-I','--InsertThreshold',dest="InsertThreshold",default=10,type='int',help="For PE library, when two mates overlap, TRUST cannot properly detect CDR3s based on mapped mates. Set larger value to force TRUST to run in PE mode despite small library insert size. Default 10.")
        parser.add_option('-o','--OutputDirectory',dest="WD",default="./",help="Directory to store intermediate files and final TRUST report. User must have write privilege. If omitted, the current directory will be applied.")
        parser.add_option('-B','--Bcell',dest="Bcell",default=False,action="store_true",help="B cell receptor inference is currently under active development.")
	return parser.parse_args()
	
def main():
	(opt,_)=CommandLineParser()
	bamdir=opt.Directory
	if len(bamdir)>0:
		files=os.listdir(bamdir)
		files0=[]
		for ff in files:
			if ff[-4:]!='.bam':
				continue
			ff=bamdir+'/'+ff
			files0.append(ff)
		files=files0
	else:
		files=[]
	File=opt.file
	if len(File)>0:
		files=[File]
	FileList=opt.files
	if len(FileList)>0:
		files=[]
		fL=open(FileList)
		for ff in fL.readlines():
                    	if ff.strip()[-4:]!='.bam':
			    continue
			files.append(ff.strip())
	mode=opt.RunningMode
	#geneFile=opt.GeneLoci
	#if len(geneFile)==0:
	#	raise "No gene annotation file specified!"
	Err=opt.Error
	thr_L=opt.Length
	fasta=opt.fasta
	IT=opt.InsertThreshold
	Bcell=opt.Bcell
	SE=opt.SingleEnd
	HC=opt.HeavyChain
	WD=opt.WD
	#if len(WD)==0:
        #    WD='./'
        if Bcell:
            print "----- Run BCR analysis -----"
            LocusFile=os.path.join(_static_path_dir,'Bcell.bed')
            geneFile=os.path.join(_static_path_dir,'BCRall.bed')
        else:
            LocusFile=os.path.join(_static_path_dir,'Tcell.bed')
            geneFile=os.path.join(_static_path_dir,'TCRall.bed')
	gL,gN=ParseGeneLoci(geneFile)
	for ff in files:
		print ff
		fname=ff.split('/')[-1]
		if not fasta:
                        gC=open(WD+'/'+fname+'_coverage','w')
                else:
                        gC=open(WD+'/'+fname+'.fa','w')
                print "Paired-end mode"
                print "Screen for informative reads"
                if SE:
                    sr=1
                else:
                    try:
                        sr=ScreenGenome(ff,LocusFile,IT,WD)
                    except ValueError:
                        continue
                if sr==0:   ## Paired end mode checks out
                    ffm=WD+'/'+fname+'-Locs.bam'
                    ffu=WD+'/'+fname+'-unmapped.bam'
                    hh=pysam.Samfile(ffm)
                    REFs=hh.references
                    # Extract all reads from the bam file without index
                    all_reads=[]
                    count=0
                    oversize=0
                    for rr in hh.fetch(until_eof=True):
                            all_reads.append(rr)
                            count+=1
                            if count>=500000:
                                    oversize=1
                                    break
                    if oversize==1:
                            continue
                    N_all=len(all_reads)
                    if N_all==0:
                        continue
                    nr=all_reads[0].rlen
                    if mode not in ['Cov','Full']:
                            raise ''' the value of mode can only be Cov or Full, Cov: return coverage on each gene only; Full: assemble sequences ''' 
                    rD,pRD,gD=AllocateReadsIntoGenes(all_reads,gL,gN,REFs)
                    LrD=0
                    for kk in rD:
                        LrD+=len(rD[kk])
                    print "Total TCR reads extracted %d" %(LrD)
                    print "Informative read pairs kept %d" %(len(pRD))
                    if mode=='Cov':
                            gC.write('#Gene\tCoverage\n')
                            for kk in gD:
                                    gC.write(kk+'\t'+str(len(gD[kk]))+'\n')
                    if mode=='Full':
                            MasterContig={}
                            UR_all=[]
                            if not fasta:
                                    gC.write('#File name\tRelative Expression\tContig length\tRead count in Locus\tVgene\tJgene\tMappedGene\t+CDR3 amino acid\tQuality Score\tCDR3 DNA\tContig DNA sequence\n')
                            else:
                                    gC.write('''## Command: python TRUST.py -f %s -d %s -F %s -m %s -e %d -l %d -a %s -s %s -H %s -I %d -o %s -B %s
## Information line contains the following fields:
# File name
# Normalized read count, or relative expression
# Contig sequence length
# Total TCR reads count
# TRUST annotated variable gene
# TRUST annotated joining gene
# Aligner reported gene (PE mode only)
# CDR3 amino acid sequence
# -log(E value), QC measure for mapping CDR3 contig to IMGT reference
# CDR3 DNA sequence
''' %(File,bamdir,FileList,mode,Err,thr_L,fasta,SE,HC,IT,WD,Bcell))
                            for kk in gD:
                                    print "---Processing %s." %(kk)
                                    vv=gD[kk]
                                    OI,UR=GetReadsOverlapByGene(vv,nr,Err,thr_L)
                                    UR_all+=UR
                                    DC,OD=FindDisjointCommunities(OI)
                                    if fasta:
                                            Contigs=[]
                                    else:
                                            Contigs=''
                                    for dc in DC:
                                            ASs=AssembleCommReads(dc,OD,OI,pRD,rD,nr=nr,thr_overlap=thr_L)
                                            for uu in ASs:
                                                    nk=ASs[uu][2]
                                                    if fasta:
                                                            Contigs.append(([uu,nk],ASs[uu][0]+ASs[uu][1],ASs[uu][3]))
                                                    else:
                                                            Contigs+=uu+':'+str(nk)+'_'
                                    MasterContig[kk]=Contigs
                            time1=time.time()
                            ContigFinalList=MergeMasterContig(MasterContig,pRD,rD,thr_L,Err)
                            time2=time.time()
                            print '...time elapsed %f' %(time2-time1)
                            print "Final contigs assembled %d" %(len(ContigFinalList))
                            annList=[]
                            for cc in ContigFinalList:
                                tmpSeq=cc[0]
                                ann=AnnotateCDR3(tmpSeq,pRD,cc[3],error=Err,overlap_thr=thr_L,Bcell=Bcell)
                                annList.append(ann)
                            print "Process unmapped read pairs"
                            annListU,ContigFinalListU,N_all_u=ProcessSingleEndReads(ffu,LocusFile,HeavyChain=HC,err=Err,overlap_thr=thr_L,fasta=fasta,unmapped=True,Bcell=Bcell)
                            annList+=annListU
                            ContigFinalList+=ContigFinalListU
                            N_all+=N_all_u
                if sr==1:
                    print "Switching to single end mode"
                    annList,ContigFinalList,N_all=ProcessSingleEndReads(ff,LocusFile,HeavyChain=HC,err=Err,overlap_thr=thr_L,fasta=fasta,Bcell=Bcell)
                    print "Number of annList %d" %(len(annList))
                for ii in xrange(0,len(annList)):
                    ann=annList[ii]
                    if len(ann)>0:
                        cc=ContigFinalList[ii]
                        if sr==1 or SE:
                            motifGene=''
                        else:
                            motifGene=cc[2]
                        if fasta:
                            tmpInfo='>'+fname+'+est_clonal_freq='+str(ann[4])+'+seq_length='+str(len(cc[0]))+'+est_lib_size='+str(N_all)+'+'+ann[2]+'+'+ann[3]+'+'+motifGene+'+'+ann[0]+'+minus_log_Eval='+str(ann[5])+'+'+ann[1]+'\n'
                            gC.write(tmpInfo)
                            gC.write(cc[0]+'\n')
                        else:
                            tmpInfo=fname+'\t'+'\test_clonal_freq='+str(ann[4])+'\tseq_length='+str(len(cc[0]))+'\test_lib_size='+str(N_all)+'\t'+ann[2]+'\t'+ann[3]+'\t'+motifGene+'\t'+ann[0]+'\tminus_log_Eval='+str(ann[5])+'\t'+ann[1]+'\t'+cc[0]+'\n'
                            gC.write(tmpInfo)
		gC.close()

if __name__=='__main__':
	main()
	print "Total time elapsed: %f" %(time.time()-t0)
	print "Maximum memory usage: %f MB" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000)
					
							
