#!/usr/bin/env python

import sys
import os, errno
import argparse
import time

# idée: on extrait les séquences avec la representative et les contigs
start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering catalog from contigs alignments. Script was writen by C. Pouchon (2021).')
parser.add_argument("-i","--infile", help="fastq consensus file",
                    type=str)
parser.add_argument("-o","--outpath", help="outpath directory",
                    type=str)
parser.add_argument("--refcov", help="minimum percentage of the reference covered to process a locus",
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


input_file=args.infile
outpath=args.outpath
ref_lenseq=args.refcov


mkdir(str(outpath+"files"))

nameofsample=input_file.split(".fq")[0]

seq_dict=dict()

with open(input_file) as f:
    reftab = f.readlines()
for line in reftab:
    if line.startswith("@") and "NODE_" in line:
        l=line.rstrip()
        refname=l.replace("@","")
        reflength=int(l.split("_length_")[1].split("_cov_")[0])
        seq=list()
    elif line.startswith("+"):
        l=line.rstrip()
        if l=="+":
            cons="".join(seq)
            #sequence=cons.upper()
            diff=reflength-len(cons)
            sequence=str(cons+("N"*diff)).upper()
            ratio=(len(sequence)-sequence.count("N"))/reflength
            if ratio>=ref_lenseq:
                seq_dict[refname]="".join(sequence)
            else:
                pass
        else:
            pass
    else:
        l=line.rstrip()
        seq.append(l)


for ref in list(seq_dict.keys()):
    fname=ref+".fa"
    header=">"+nameofsample
    seq=seq_dict[ref]
    if os.path.isfile(os.path.join(outpath+"/files", fname)):
        with open(os.path.join(outpath+"/files", fname), 'a+') as file:
            old_headers = []
            end_file=file.tell()
            file.seek(0)
            for line in file:
                if line.startswith(">"):
                    old_headers.append(line.rstrip().replace(">","").split(";")[0])
            if not nameofsample in old_headers:
                file.seek(end_file)
                file.write(header+'\n')
                file.write(str(seq)+'\n')
            else:
                pass
    else:
        with open(os.path.join(outpath+"/files", fname), 'w') as out:
            out.write(header+'\n')
            out.write(str(seq)+'\n')
