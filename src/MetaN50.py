#!/usr/bin/env python

import glob
import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Selection of best metassembly according to the N50. Script was writen by C. Pouchon (2020).')
parser.add_argument('--input', help='path to metassembly directories', type=str)
parser.add_argument("--output", help="path to outpout best metacontigs",
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


inpath=args.input
outpath=args.output

in_files = []
path_to_seek=os.walk(os.path.join(inpath))
for r, d, f in path_to_seek:
    for file in f:
        if "unclean_catalog_" in file:
                in_files.append(os.path.join(r, file))

mkdir(outpath)

meta_infos={}
# parsing of metassemblies
for file in in_files:
    meta_infos[file]=[]
    with open(file) as f:
        reftab = f.readlines()
    for line in reftab:
        if line.startswith(">"):
            l=line.rstrip()
            length=int(l.split("_length_")[1].split("_cov")[0])
            meta_infos[file].append(length)

N50={}
L50={}
print("file","N50", "L50", "tot. size","tot. contigs")
# sort assemblies length and N50 computing
for file in meta_infos:
    meta_infos[file].sort(reverse=True)
    totlen=sum(meta_infos[file])
    half=totlen/2
    totcont=len(meta_infos[file])
    #init loop to check the N50 values
    size=0
    i=0
    for s in meta_infos[file]:
        cont_size=s+size
        if cont_size<half:
            size=cont_size
            i=i+1
            continue
        else:
            mN50=s
            mL50=i
            break
    print(file, mN50, mL50,totlen,totcont)
    N50[file]=mN50
    L50[file]=mL50
