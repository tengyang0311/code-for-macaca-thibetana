#!/usr/bin/python

import sys
#import gzip

file1 = sys.argv[1]
file2 = sys.argv[2]

f0 = open('/media/primates/tythi/ref/sortlenthofgenome.fa.fai','r')
dict1 = {}
dict2 = {}
dict3 = {}
for x in f0:
    a = x.strip().split('\t')
    ID = int(a[1])/100000
    for i in range(0,ID+1):
         dict2[a[0]+'_'+str(i)] = 0
         dict3[a[0]+'_'+str(i)] = 0
         if i != ID:
             dict1[a[0]+'_'+str(i)] = a[0]+'\t'+str(i*100000+1)+'\t'+str((i+1)*100000)+'\t100000'
         else:
              dict1[a[0]+'_'+str(i)] = a[0]+'\t'+str(i*100000+1)+'\t'+a[1]+'\t'+str(int(a[1])-i*100000)

f1 = open(file1,'rb')
f2 = open(file2,'w')

for eachline in f1:
    if '#' in eachline:
        continue
    elif 'KQ' in eachline:
        continue
    elif 'JSUE' in eachline:
        continue
    else:
        a = eachline.strip().split('\t')
        if a[6] == 'PASS':
            ID = int(a[1])/100000
            key = a[0]+'_'+str(ID)
            b = a[9].split(':')[0]
            if b == '0/1':
                dict2[key]+= 1
            elif b == '1/1':
                dict3[key]+= 1

for a in dict1:
    c = dict1[a].split('\t')[3]
    print >> f2, a+'\t'+dict1[a]+'\t'+str(dict2[a])+'\t'+str(dict3[a])+'\t'+str(float(dict2[a])/int(c))
