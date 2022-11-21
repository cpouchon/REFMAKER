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
parser.add_argument("--outgroups", help="list of outgroup taxa",
                    type=str)
parser.add_argument("-o","--outpath", help="output directory",
                    type=str)

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

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


inpath=args.inpath
outpath=args.outpath
list_outgroups = [str(item) for item in args.outgroups.split(',')]


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


fname="concatenated.nex"
open(os.path.join(outpath, fname), 'w').close()

with open(os.path.join(outpath,fname), 'a+') as file:
    file.write("#NEXUS"+"\n")
    file.write("begin data;"+"\n")
    file.write(str("\t"+"dimensions ntax="+str(tot)+" nchar="+str(size)+";"+"\n"))
    file.write("\t"+"format datatype=nucleotide missing=? gap=-;"+"\n")
    file.write("matrix"+"\n")

for s in ref_catalog:
    seqID=s.replace("-","_")
    cons=str(ref_catalog[s].seq)
    with open(os.path.join(outpath,fname), 'a+') as file:
        file.write("\t"+str(seqID)+"\n")
        file.write(str(cons.replace("N","?"))+"\n")

with open(os.path.join(outpath,fname), 'a+') as file:
    file.write("\t"+";"+"\n")
    file.write("end;"+"\n")
    file.write("\n")
    file.write("begin sets;"+"\n")

for p in list(paritions_infos):
    id=p["ref"]
    s=p["start"]
    e=p["end"]
    with open(os.path.join(outpath,fname), 'a+') as file:
        file.write("\t"+"charset "+str(id)+" = "+str(s)+"-"+str(e)+";"+"\n")

with open(os.path.join(outpath,fname), 'a+') as file:
    file.write("\n")
    file.write("\t"+"charpartition combined = "+str(", ".join(combined) )+";"+"\n")
    file.write("end;"+"\n")
    file.write("\n")
    file.write("begin paup;"+"\n")
    file.write("\t"+"outgroup "+str(" ".join(list_outgroups))+";"+"\n")
    file.write("\t"+"set outroot=mono;"+"\n")
    file.write("end;"+"\n")
