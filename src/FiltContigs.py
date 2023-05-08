#!/usr/bin/env python

import glob
import argparse
import sys
import os, errno


parser = argparse.ArgumentParser(description='Identification of contaminant according to taxonomic assignment of contigs into rRNA database. Script was writen by C. Pouchon (2020).')
parser.add_argument("--blast", help="output from blast of contigs into database",
                    type=str)
parser.add_argument('--database', help='path to databases', type=str)
parser.add_argument('--outpath', help='path to outpout contig infos', type=str)


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

blastout=args.blast
dbpath=args.database
outpath=args.outpath

contigs={}
contigs_align={}
contigs_pos={}
toremove=list()
cont_to_remove=list()
bitscore_dict={}
best_hit={}

in_files = []
path_to_seek=os.walk(os.path.join(dbpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(".fasta"):
            in_files.append(os.path.join(r, file))

rdna=list()
cpdna=list()
mtdna=list()
others=list()
for f in in_files:
    if any(x in f for x in ("rfam-5s","_nucrdna","euk-28s","rfam-5.8s","euk-18s")):
        with open(f) as ffile:
            for l in ffile:
                if l.startswith(">"):
                    id=l.rstrip().replace(">","").split(" ")[0]
                    rdna.append(id)
                else:
                    pass
    elif "mitochondrion" in f:
        with open(f) as ffile:
            for l in ffile:
                if l.startswith(">"):
                    id=l.rstrip().replace(">","").split(" ")[0]
                    mtdna.append(id)
                else:
                    pass
    elif "chloroplast" in f:
        with open(f) as ffile:
            for l in ffile:
                if l.startswith(">"):
                    id=l.rstrip().replace(">","").split(" ")[0]
                    cpdna.append(id)
                else:
                    pass
    else:
        with open(f) as ffile:
            for l in ffile:
                if l.startswith(">"):
                    id=l.rstrip().replace(">","").split(" ")[0]
                    others.append(id)
                else:
                    pass

with open(blastout) as blastf:
    for l in blastf:
        contigid=l.rstrip().split("\t")[0]
        bitscore=float(l.rstrip().split("\t")[11])
        aln=int(l.rstrip().split("\t")[3])
        start=int(l.rstrip().split("\t")[6])
        end=int(l.rstrip().split("\t")[7])
        refid=l.rstrip().split("\t")[1]
        if contigid in list(contigs.keys()):
            if refid in list(contigs[contigid].keys()):
                tmp=contigs_pos[contigid][refid]
                cond_intersect=0
                for p in contigs_pos[contigid][refid]:
                    min=int(p.split("-")[0])
                    max=int(p.split("-")[1])
                    rhit=set(range(min,max))
                    rc=set(range(start,end))
                    intersection=rhit.intersection(rc)
                    if len(intersection)==0:
                        pass
                    else:
                        cond_intersect=cond_intersect+1
                        continue
                if cond_intersect==0:
                    old_score=float(contigs[contigid][refid])
                    new_score=float(old_score)+float(bitscore)
                    old_aln=int(contigs_align[contigid][refid])
                    new_aln=int(old_aln)+int(aln)
                    contigs[contigid][refid]=new_score
                    contigs_align[contigid][refid]=new_aln
                    contigs_pos[contigid][refid].append(str(start)+"-"+str(end))
            else:
                contigs[contigid][refid]=bitscore
                contigs_align[contigid][refid]=aln
                contigs_pos[contigid][refid]=[]
                contigs_pos[contigid][refid].append(str(start)+"-"+str(end))
        else:
            contigs.setdefault(contigid, {})
            contigs_align.setdefault(contigid, {})
            contigs_pos.setdefault(contigid, {})
            contigs[contigid][refid]=bitscore
            contigs_align[contigid][refid]=aln
            contigs_pos[contigid][refid]=[]
            contigs_pos[contigid][refid].append(str(start)+"-"+str(end))

for allcont in contigs.keys():
    for elem in contigs[allcont].keys():
        if allcont in bitscore_dict.keys():
            score=contigs[allcont][elem]
            if float(score) >= float(bitscore_dict[allcont]):
                best_hit[allcont]=elem
                bitscore_dict[allcont]=float(score)
            else:
                pass
        else:
            score=contigs[allcont][elem]
            best_hit[allcont]=elem
            bitscore_dict[allcont]=float(score)

for cont in best_hit.keys():
    if cont in cont_to_remove:
        continue
    else:
        refcont=best_hit[cont]
        toprint=list()
        if cont in cont_to_remove:
            pass
        else:
            toprint.append(cont)
            toprint.append(refcont)
            toremove.append(toprint)
            cont_to_remove.append(cont)
            align_cont=contigs_align[cont][refcont]
            cont_size=int(cont.split("length_")[1].split("_")[0])
            if align_cont/cont_size >=0.0:
                if refcont in rdna:
                    with open(os.path.join(str(outpath),'rdna_contigs.infos'), 'a+') as out:
                        out.write(cont+"\n")
                        out.close()
                elif refcont in cpdna:
                    with open(os.path.join(str(outpath),'cpdna_contigs.infos'), 'a+') as out:
                        out.write(cont+"\n")
                        out.close()
                elif refcont in mtdna:
                    with open(os.path.join(str(outpath),'mtdna_contigs.infos'), 'a+') as out:
                        out.write(cont+"\n")
                        out.close()
                elif refcont in others:
                    with open(os.path.join(str(outpath),'others_contigs.infos'), 'a+') as out:
                        out.write(cont+"\n")
                        out.close()
            else:
                pass

if len(toremove)>0:
    for l in toremove:
        print("\t".join(l))
else:
    print("none")
