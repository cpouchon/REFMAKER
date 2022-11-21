#!/usr/bin/env python
from Bio import SeqIO

import sys
import os, errno
import argparse
import time


start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering consensus sequences from catalog alignments - step 1. Script was writen by C. Pouchon (2021).')
parser.add_argument("-i","--inpath", help="path to fasta consensus",
                    type=str)
parser.add_argument("-o","--outpath", help="output directory",
                    type=str)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

class ProgressBar:
	'''
	Progress bar
	'''
	def __init__ (self, valmax, maxbar, title):
		if valmax == 0:  valmax = 1
		if maxbar > 200: maxbar = 200
		self.valmax = valmax
		self.maxbar = maxbar
		self.title  = title
	def update(self, val):
		import sys
		perc  = round((float(val) / float(self.valmax)) * 100)
		scale = 100.0 / float(self.maxbar)
		bar   = int(perc / scale)
		out = '\r %20s |%s%s| %3d %% ' % (self.title, '.' * bar, ' ' * (self.maxbar - bar), perc)
		sys.stdout.write(out)
		sys.stdout.flush()

def mkdir(path, overwrite=False):
    '''
    function to create a directory for output fasta
    '''
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            if not overwrite:
                pass#print ("path '%s' already exists" % path)   # overwrite == False and we've hit a directory that exists
        else: raise

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


inpath=args.inpath
outpath=args.outpath

ref_list=list()
ref_len={}

paritions_infos=list()
combined=list()
with open(os.path.join(inpath, "concatenated.info")) as tab:
    for l in tab:
        line=l.rstrip().split("\t")
        ref=line[3]
        start=line[0]
        end=line[1]
        tmp=dict()
        tmp["ref"]=ref
        tmp['start']=start
        tmp['end']=end
        paritions_infos.append(tmp)
        combined.append(str(ref)+":"+str(start)+"-"+str(end))



ref_catalog=to_dict_remove_dups(SeqIO.parse(os.path.join(inpath, "concatenated.fa"), "fasta"))
tot=len(ref_catalog)
size=len(ref_catalog[list(ref_catalog.keys())[0]].seq)

ht_basis=["R","Y","S","W","K","M","B","D","H","V"]


fname="summary.txt"
open(os.path.join(outpath, fname), 'w').close()

with open(os.path.join(outpath,fname), 'a+') as file:
    file.write("Name"+"\t"+"Seqs"+"\t"+"Sites"+"\t"+"Infor"+"\t"+"Invar"+"\t"+"Amb"+"\t"+"Miss"+"\t"+"MissAll"+"\n")

Bar = ProgressBar(len(paritions_infos), 60, '\t ')
barp=0

print("Get summary statistics over Loci alignments")
for g in paritions_infos:
    barp=barp+1
    taxa=0
    totmiss=0
    totlen=0
    totinvar=0
    totinfo=0
    totamb=0
    bmin=int(g["start"])-1
    bmax=int(g["end"])
    window=len(range(bmin,bmax))
    miss_samp=0
    for s in ref_catalog:
        seq=str(ref_catalog[s].seq)[bmin:bmax]
        miss=seq.count("-")
        if miss==len(seq):
            miss_samp+=1
    for b in range(bmin,bmax):
        tmpsite=list()
        tmpsite_without_miss=list()
        for s in ref_catalog:
            totlen+=1
            tmpsite.append(str(ref_catalog[s].seq)[b])
            if str(ref_catalog[s].seq)[b]=="-":
                pass
            else:
                if str(ref_catalog[s].seq)[b].upper() in ht_basis:
                    pass
                else:
                    tmpsite_without_miss.append(str(ref_catalog[s].seq)[b])
        amb=sum([tmpsite.count(k) for k in ht_basis])
        totamb+=amb
        miss=tmpsite.count("-")
        totmiss+=miss
        sumsite=dict((x,tmpsite_without_miss.count(x)) for x in set(tmpsite_without_miss))
        if len(sumsite)>1:
            #totvar+=1
            cond_shared=0
            for p in sumsite:
                if sumsite[p]>1:
                    cond_shared+=1
                else:
                    continue
            if cond_shared>=2:
                totinfo+=1
        else:
            totinvar+=1

    famb=round((totamb/totlen)*100,2)
    fmiss_tot=round((totmiss/totlen*100),2)

    fmiss_sel=round(((totmiss-(miss_samp*window))/totlen)*100,2)
    with open(os.path.join(outpath,fname), 'a+') as file:
        file.write(g["ref"]+"\t"+str(len(ref_catalog)-miss_samp)+"\t"+str(window)+"\t"+str(totinfo)+"\t"+str(totinvar)+"\t"+str(famb)+"\t"+str(fmiss_sel)+"\t"+str(fmiss_tot)+"\n")
    Bar.update(barp)
print("")
print("done")
