#!/usr/bin/python3.6
#######################################
########qijiwei love wangrui###########
#######################################
Usage="""
culculateLOF.py - version 1
read two vcf file(LOF.vcf from grep "LOF";intergenic.vcf from grep "intergenic") and culculate RA/B by gorilla paper and give jackknife file.

Usage:
             python culculateLOF.py -lof LOF.vcf -intergenic intergenic.vcf -o output.txt -o2 jackknife.txt"""

import argparse
import os
import operator
import math
parser=argparse.ArgumentParser(description="choose lof file and intergenic file")
parser.add_argument('-lof',help='lof vcf from grep SNPEFF result')
parser.add_argument('-intergenic',help='intergenic vcf from grep SNPEFF reuslt')
parser.add_argument('-o',help='output file')
parser.add_argument('-o2',help='output jackknife')
args=parser.parse_args()
LLOnumerator=[]
LLOdenominator=[]
LLTnumerator=[]
LLTdenominator=[]
LTOnumerator=[]
LTOdenominator=[]
LOLnumerator=[]
LOLdenominator=[]
LTLnumerator=[]
LTLdenominator=[]
LOTnumerator=[]
LOTdenominator=[]
with open(args.lof,'r') as f1:
    f1 = f1.read().split('\n')
    for line in f1[1:-1]:
        other=[]
        TH=[]
        line = line.split('\t')
        for i in line[9:16]:
            if i[0:3] == '0/1':
                other.append(int(1))
            elif i[0:3] == '1/1':
                other.append(int(2))
        for g in line[16:24]:
            if g[0:3] == '0/1':
                TH.append(int(1))
            elif g[0:3] == '1/1':
                TH.append(int(2))
        for h in line[24:27]:
            if h[0:3] == '0/1':
                other.append(int(1))
            elif h[0:3] == '1/1':
                other.append(int(2))
        for i in line[27:29]:
            if h[0:3] == '0/1':
                TH.append(int(1))
            elif h[0:3] == '1/1':
                TH.append(int(2))
        for h in line[29:35]:
            if h[0:3] == '0/1':
                other.append(int(1))
            elif h[0:3] == '1/1':
                other.append(int(2))
        Fother=sum(other)/32
        FTH=sum(TH)/20
        fsiteOT=Fother*(1-FTH)
        fsiteTO=FTH*(1-Fother)
        LTOnumerator.append(fsiteTO)
        LOTnumerator.append(fsiteOT)
with open(args.intergenic,'r') as f2:
    f2 = f2.read().split('\n')
    for line in f2[:-1]:
        other=[]
        TH=[]
        line = line.split('\t')
        for i in line[9:16]:
            if i[0:3] == '0/1':
                other.append(int(1))
            elif i[0:3] == '1/1':
                other.append(int(2))
        for g in line[16:24]:
            if g[0:3] == '0/1':
                TH.append(int(1))
            elif g[0:3] == '1/1':
                TH.append(int(2))
        for h in line[24:27]:
            if h[0:3] == '0/1':
                other.append(int(1))
            elif h[0:3] == '1/1':
                other.append(int(2))
        for i in line[27:29]:
            if h[0:3] == '0/1':
                TH.append(int(1))
            elif h[0:3] == '1/1':
                TH.append(int(2))
        for h in line[29:35]:
            if h[0:3] == '0/1':
                other.append(int(1))
            elif h[0:3] == '1/1':
                other.append(int(2))
        Fother=sum(other)/32
        FTH=sum(TH)/20
        fsiteOT=Fother*(1-FTH)
        fsiteTO=FTH*(1-Fother)
        LTOdenominator.append(fsiteTO)
        LOTdenominator.append(fsiteOT)
LTO=sum(LTOnumerator)/sum(LTOdenominator)
LOT=sum(LOTnumerator)/sum(LOTdenominator)
RTO=LTO/LOT
with open(args.o,'w') as f3:
    f3.write('RTO=\t'+str(RTO)+'\n')
LTOnumerator=[]
LOTnumerator=[]
with open(args.lof,'r') as f4:
    f4 = f4.read().split('\n')
    linenum = 0
    all=os.popen('wc -l '+args.lof).read().split(' ')[0]
    part=int(all)/100
    dictTO={}
    dictOT={}
    for i in range(0,100):
        dictTO[i]=[]
        dictOT[i]=[]
    for line in f4[1:-1]:
        linenum += 1
        block=math.floor(linenum/part)
        other=[]
        TH=[]
        line = line.split('\t')
        for i in line[9:16]:
            if i[0:3] == '0/1':
                other.append(int(1))
            elif i[0:3] == '1/1':
                other.append(int(2))
        for g in line[16:24]:
            if g[0:3] == '0/1':
                TH.append(int(1))
            elif g[0:3] == '1/1':
                TH.append(int(2))
        for h in line[24:27]:
            if h[0:3] == '0/1':
                other.append(int(1))
            elif h[0:3] == '1/1':
                other.append(int(2))
        for i in line[27:29]:
            if h[0:3] == '0/1':
                TH.append(int(1))
            elif h[0:3] == '1/1':
                TH.append(int(2))
        for h in line[29:35]:
            if h[0:3] == '0/1':
                other.append(int(1))
            elif h[0:3] == '1/1':
                other.append(int(2))
        Fother=sum(other)/32
        FTH=sum(TH)/20
        fsiteOT=Fother*(1-FTH)
        fsiteTO=FTH*(1-Fother)
        LTOnumerator.append(fsiteTO)
        LOTnumerator.append(fsiteOT)
        dictTO[block].append(fsiteTO)
        dictOT[block].append(fsiteOT)
    with open(args.o2,'w') as f4:
        f4.write('RTO\n')
        for f in range(0,100):
            LTOup=0
            LOTup=0
            for m in range(0,f):
                LTOup=LTOup+sum(dictTO[m])
                LOTup=LOTup+sum(dictOT[m])
            for o in range(f+1,100):
                LTOup=LTOup+sum(dictTO[o])
                LOTup=LOTup+sum(dictOT[o])
            LTO=LTOup/sum(LTOdenominator)
            LOT=LOTup/sum(LOTdenominator)
            RTO=LTO/LOT
            f4.write(str(RTO)+'\n')
        f4.close()
