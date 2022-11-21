#!/usr/bin/env python

import sys
import os, errno
import argparse
import time

# idée: on extrait les séquences avec la representative et les contigs
start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering catalog from contigs alignments. Script was writen by C. Pouchon (2021).')
parser.add_argument("-p","--pop", help="population file",
                    type=str)
parser.add_argument("-c","--clust", help="cd-hit clustering file",
                    type=str)
parser.add_argument("-o","--outpath", help="outpath directory",
                    type=str)
parser.add_argument("-f","--fasta", help="fasta input file of cd-hit clusters",
                    type=str)
parser.add_argument("-l","--length", help="minimal length of representative cluster",
                    type=int)
parser.add_argument("-r", help="minimum percentage of individuals in a population required to process a locus",
                    type=float)
parser.add_argument("-R", help="minimum percentage of populations required to process a locus",
                    type=float)
if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

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


clust_file=args.clust
outpath=args.outpath
infile=args.fasta
within_clust_pop=args.r
between_clust_pop=args.R
pop_file=args.pop
lmin=args.length

print("min. frequency of overlapping contig within population per locus: %s" % (str(within_clust_pop)))
print("min. frequency of populations sharing a locus: %s" % (str(between_clust_pop)))

mkdir(str(outpath))

populations_inds=dict()
with open(pop_file) as p:
    poptab = p.readlines()
for line in poptab:
    l=line.rstrip()
    ind=l.split("\t")[0]
    pop=l.split("\t")[1]
    if pop in list(populations_inds.keys()):
        populations_inds[pop].append(ind)
    else:
        populations_inds[pop]=[]
        populations_inds[pop].append(ind)


seqs={}
cond_init=0
with open(infile) as f:
    reftab = f.readlines()
for line in reftab:
    if line.startswith(">"):
        l=line.rstrip()
        if cond_init==0:
            cond_init=cond_init+1
        else:
            if seqID not in list(seqs.keys()):
                cons="".join(s)
                seqs[seqID]=cons.replace("N","")
            else:
                pass
            cond_init=cond_init+1
        seqID=l.replace(">","")
        s=[]
    else:
        l=line.rstrip()
        s.append(l)
        if line==reftab[len(reftab)-1]:
            if seqID not in list(seqs.keys()):
                cons="".join(s)
                seqs[seqID]=cons
            else:
                pass
        else:
            pass

totseq=len(list(seqs.keys()))

with open(clust_file) as f:
    clustab = f.readlines()
refnames=list()
clust_dict=dict()
cond_init=0
for line in clustab:
    if line.startswith(">"):
        l=line.rstrip()
        if cond_init==0:
            cond_init=cond_init+1
        else:
            cond_pop=0
            for p in list(populations_inds.keys()):
                ratio=len(list(set(populations_inds[p]).intersection(inds)))/len(populations_inds[p])
                if ratio>=within_clust_pop:
                    cond_pop=cond_pop+1
                else:
                    pass
            ratio_btpop=cond_pop/(len(populations_inds.keys()))
            if ratio_btpop>=between_clust_pop:
                clust_dict[clust_ref]=subnames
            else:
                pass
            cond_init=cond_init+1
        inds=[]
        subnames=[]
    else:
        l=line.rstrip()
        subl=l.split(" ")
        if subl[len(subl)-1]=="*":
            clust_ref=subl[1].split(">")[1].split("...")[0]
        else:
            name=subl[1].split(">")[1].split("...")[0]
            ind=name.split("_NODE")[0]
            if ind in inds:
                pass
            else:
                inds.append(ind)
            subnames.append(name)


sel_tax=len(list(clust_dict.keys()))

with open(os.path.join(outpath, "unclean_catalog.fa"), 'w') as file:
    file.close()

pass_length_filter=0
for ref in list(clust_dict.keys()):
    reflen=int(ref.split("_length_")[1].split("_cov")[0])
    header=">"+ref
    refseq=seqs[ref]
    if reflen>=lmin:
        pass_length_filter=pass_length_filter+1
        with open(os.path.join(outpath, "unclean_catalog.fa"), 'a+') as file:
            file.write(header+'\n')
            file.write(str(refseq)+'\n')
            file.close()
    else:
        pass

print("%s/%s clusters passed the populations thresholds" % (str(sel_tax),str(totseq)))
print("%s/%s clusters passed the minimum length threshold" % (str(pass_length_filter),str(sel_tax)))




# outputs:
# 1) list des reprs
# 2) directory rep.fa and seq_rep.fa pour winnowmap mapping?

# unclean --> with popcondition (tester avec celui-ci le)
# clean catalog --> with winowmapp and pop condition (need overlapp) sur ref overlapp 25%
