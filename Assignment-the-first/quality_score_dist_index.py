#!/usr/bin/env python

import Bioinfo
import numpy as np
import matplotlib.pyplot as plt 
import gzip 
import argparse


def get_args():             
    """This function passes arguements to a python script to"""
    parser = argparse.ArgumentParser(description="A program to normalize kmerspec")
    parser.add_argument("-f", "--filename", help="Your filename", required=True)               
    parser.add_argument("-o", "--output_filename", help="what is your output file?", required=True)
    return parser.parse_args()
	
args = get_args()    
filename = (args.filename)
output = (args.output_filename)


all_qscores=np.zeros((8),dtype=float)

i=0
with gzip.open(filename, "rt") as fh: 
    rownum=0
    for line in fh:
        line=line.strip('\n')
        i+=1
        if i%4 == 0:
            for y, x in enumerate(line):
                #print(y, Bioinfo.convert_phred(x))
                all_qscores[y]+=Bioinfo.convert_phred(x)
                #print(y)
            rownum+=1

print("#Base Pair"'\t'"Mean Quality Score")
for t, phredsum in enumerate(all_qscores):
    all_qscores[t] = (phredsum / (i/4))
    print(t,"\t",all_qscores[t])            
                
    values = all_qscores

plt.figure(figsize=(15,4))
plt.ylabel('mean quality score')
plt.xlabel('nt position')
plt.bar(range(0, 8), values)
plt.savefig("{}.png".format(output))

