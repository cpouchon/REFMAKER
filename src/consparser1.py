#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord import *
from Bio.Seq import *
from Bio.SeqUtils import *
import sys
import os, errno
import argparse
import time
from scipy import stats
import numpy as np

start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering consensus sequences from catalog alignments - step 1. Script was writen by C. Pouchon (2021).')
parser.add_argument("-f","--fasta", help="path to fasta consensus",
                    type=str)
parser.add_argument("-o","--output_dir", help="output directory",
                    type=str)
parser.add_argument("-c","--catalog", help="reference catalog fasta-file",
                    type=str)
parser.add_argument("-m","--min_seq", help="minimal sequences required",
                    type=int)

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


inpath=args.fasta
ref=args.catalog
outpath=args.output_dir
cond_m=args.min_seq

mkdir(outpath)

# parsing of catalog / add ref seq in the outfile (for alignment)
ref_list=list()
ref_len={}
cond_init=0
count=0
out_catalog={}
seqs={}

print("parsing catalog")

ref_catalog=to_dict_remove_dups(SeqIO.parse(ref, "fasta"))
for s in ref_catalog:
    seqID=s
    if seqID not in list(out_catalog.keys()):
        out_catalog[seqID]={}
        if seqID not in ref_list:
            cons=str(ref_catalog[s].seq)
            ref_len[seqID]=len(cons)
            out_catalog[seqID][seqID]=cons.upper()
            ref_list.append(seqID)


print("done")

print("parsing sample files")
# we parse all fasta files
in_files = []
path_to_seek=os.walk(os.path.join(inpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(".fa"):
            in_files.append(os.path.join(r, file))

for file in in_files:
    sampleID=os.path.basename(file).replace(str("."+"fa"),"")
    sample_catalog=to_dict_remove_dups(SeqIO.parse(file, "fasta"))
    seq_processed=list()
    for s in sample_catalog:
        seqID=s
        if seqID not in seq_processed:
            cons=str(sample_catalog[s].seq)
            seq_mapp=len(cons)-cons.count("N")
            if seq_mapp>0:
                out_catalog[seqID][sampleID]=cons.replace("N","-").upper()
                seq_processed.append(seqID)
            else:
                seq_processed.append(seqID)
        else:
            pass
print("done")

# now we write only sequences with at least 4 samples (+ref) so 5 IDs
print("writing outfiles")
for ref in list(out_catalog.keys()):
    count=len(list(out_catalog[ref].keys()))
    if count>=cond_m:
        fname=ref+".fa"
        for seqs in out_catalog[ref]:
            header=">"+seqs
            seq=out_catalog[ref][seqs]
            if os.path.isfile(os.path.join(outpath,fname)):
                with open(os.path.join(outpath,fname), 'a+') as file:
                    old_headers = []
                    end_file=file.tell()
                    file.seek(0)
                    for line in file:
                        if line.startswith(">"):
                            old_headers.append(line.replace(">",""))
                    if not seqs in old_headers:
                        file.seek(end_file)
                        file.write(header+'\n')
                        file.write(str(seq)+'\n')
                    else:
                        pass
            else :
                with open(os.path.join(outpath,fname), 'w') as out:
                    out.write(header+'\n')
                    out.write(str(seq)+'\n')
print("done")
