#!/usr/bin/env python

import glob
import argparse
import sys
import os, errno
from ete3 import NCBITaxa



parser = argparse.ArgumentParser(description='Parsing blast output of contig mapping into the NCBI Nt database and taxonomic filtering. Script was writen by C. Pouchon (2020).')
parser.add_argument("--blast", help="output from blast of contigs into database",
                    type=str)
parser.add_argument('--rank', help='expected taxonomic rank', type=str)
parser.add_argument('--paln', help='minimal percent of contigs aligned into the ref', type=float)

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

blastout=args.blast
expected=args.rank
ratio=args.paln

contigs={}
contigs_align={}
contigs_pos={}
toremove=list()
cont_to_remove=list()
accessions={}
taxids={}
rawscore_dict={}
best_hit={}

ncbi = NCBITaxa()

# put to 10% of contigs blast

with open(blastout) as blastf:
    for l in blastf:
        contigid=l.rstrip().split("\t")[0]
        contigs.setdefault(contigid, {})
        contigs_align.setdefault(contigid, {})
        contigs_pos.setdefault(contigid, {})
        refid=l.rstrip().split("\t")[8]
        taxids.setdefault(refid, "NA")
        accessions.setdefault(refid, "NA")
        taxid_tmp=l.rstrip().split("\t")[9]
        if ";" in taxid_tmp:
            taxid=taxid_tmp.split(";")[0]
        else:
            taxid=taxid_tmp
        rawscore=float(l.rstrip().split("\t")[10])
        aln=int(l.rstrip().split("\t")[3])
        start=int(l.rstrip().split("\t")[4])
        end=int(l.rstrip().split("\t")[5])

        if contigid in contigs.keys():
            if refid in contigs[contigid].keys():
                for p in contigs_pos[contigid][refid]:
                    min=int(p.split("-")[0])
                    max=int(p.split("-")[1])
                    rhit=set(range(min,max))
                    rc=set(range(start,end))
                    intersection=rhit.intersection(rc)
                    if len(intersection)>0:
                        pass
                    else:
                        old_score=float(contigs[contigid][refid])
                        new_score=float(old_score)+float(rawscore)
                        old_aln=int(contigs_align[contigid][refid])
                        new_aln=int(old_aln)+int(aln)
                        contigs[contigid][refid]=new_score
                        contigs_align[contigid][refid]=new_aln
                        contigs_pos[contigid][refid].append(str(start)+"-"+str(end))
            else:
                contigs[contigid][refid]=rawscore
                contigs_align[contigid][refid]=aln
                contigs_pos[contigid][refid]=[]
                contigs_pos[contigid][refid].append(str(start)+"-"+str(end))
                taxids[refid]=taxid
        else:
            pass

for allcont in contigs.keys():
    for elem in contigs[allcont].keys():
        if allcont in rawscore_dict.keys():
            score=contigs[allcont][elem]
            if float(score) >= float(rawscore_dict[allcont]):
                best_hit[allcont]=elem
                rawscore_dict[allcont]=float(score)
            else:
                pass
        else:
            score=contigs[allcont][elem]
            best_hit[allcont]=elem
            rawscore_dict[allcont]=float(score)

for k in accessions.keys():
    tax=int(taxids[k])
    if len(ncbi.get_lineage(tax))==0:
        pass
    else:
        sub_lineages=ncbi.get_lineage(tax)
        sub_names=ncbi.get_taxid_translator(sub_lineages)
        sub_taxo=[value for value in sub_names.values()]
        taxo=",".join(sub_taxo)
        if expected.lower() in taxo.lower():
            accessions[k]="TRUE"
        else:
            accessions[k]="FALSE"

for cont in best_hit.keys():
    if cont in cont_to_remove:
        continue
    else:
        refcont=best_hit[cont]
        if refcont in accessions.keys():
            toprint=list()
            if accessions[refcont]=="FALSE":
                align_cont=contigs_align[cont][refcont]
                cont_size=int(cont.split("length_")[1].split("_")[0])
                if align_cont/cont_size >=ratio:
                    if cont in cont_to_remove:
                        pass
                    else:
                        toprint.append(cont)
                        toprint.append(refcont)
                        toremove.append(toprint)
                        cont_to_remove.append(cont)
                else:
                    pass
            else:
                pass
        else:
            pass

if len(toremove)>0:
    for l in toremove:
        print("\t".join(l))
else:
    print("none")
