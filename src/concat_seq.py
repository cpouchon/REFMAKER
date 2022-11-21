#! /usr/bin/env python

import glob
import sys
import os, errno
import argparse
import fnmatch



parser = argparse.ArgumentParser(description='Concatenation of sequences alignments')
parser.add_argument("-p","--path", help="searching path of sequences files",
                    type=str)
parser.add_argument("-i","--infile", help="input file list of sequence files, sequences have to be aligned before",
                    type=str)
parser.add_argument("-pfind","--pathfind", help="[mode] search all sequences files in given path (-p) otherwise parse a given list (-i)",
                    action="store_true")
parser.add_argument("-o","--outfile", help="output file list of concatenated sequences",
                    type=str)
parser.add_argument("-t","--taxa", help="list of taxa to include in the final concatenated file",
                    type=str)
parser.add_argument("--missingfract", help="maximal missing data threshold allowed to consider final sequence (e.g. 0.5 meaning that final sequence has fewer than 0.5 of missing data)",
                     type=float,default=None)

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
		out = '\r %20s [%s%s] %3d %% ' % (self.title, '.' * bar, ' ' * (self.maxbar - bar), perc)
		sys.stdout.write(out)
		sys.stdout.flush()



if args.pathfind:
    'Search all fasta files files'
    path = args.path
    in_files = []
    path_to_seek=os.walk(os.path.join(path,"files"))
    for r, d, f in path_to_seek:
        for file in f:
            if str("."+"fa") in file:
                in_files.append(os.path.join(r, file))
else:
    'We parse a list of files'
    input_list = args.infile
    with open(input_list) as f:
        in_files = f.readlines()

#cat *.fa | grep '>' | sort | uniq | perl -pe 's/>//' > list_taxa

fname=args.outfile+".fa"
Nthreshold=args.missingfract

giventaxa = args.taxa
with open(giventaxa) as f:
    taxal = f.readlines()

taxa_dict={}
taxa_list = list()
for line in taxal:
    l = line.rstrip()
    if len(l)>0:
        taxa_dict.setdefault(l, {})


end=0
partnumber=0

for file in in_files:
    seqs={}
    # if args.trim:
    #     refname=os.path.basename(file).replace(str("."+"trim"),"")
    # else:
    #     refname=os.path.basename(file).replace(str("."+"fa"),"")
    refname=os.path.basename(file).replace(str("."+"fa"),"")
    cond_init=0
    with open(file) as f:
        reftab = f.readlines()
    for line in reftab:
        if line.startswith(">"):
            l=line.rstrip()
            if cond_init==0:
                cond_init=cond_init+1
            else:
                if seqID not in list(seqs.keys()):
                    seqs[seqID]=dict()
                    if refname not in list(seqs[seqID].keys()):
                        seqs[seqID][refname]="".join(s)
                    else:
                        pass
                else:
                    if refname not in list(seqs[seqID].keys()):
                        seqs[seqID][refname]="".join(s)
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
                    seqs[seqID]=dict()
                    if refname not in list(seqs[seqID].keys()):
                        seqs[seqID][refname]="".join(s)
                    else:
                        pass
                else:
                    if refname not in list(seqs[seqID].keys()):
                        seqs[seqID][refname]="".join(s)
                    else:
                        pass
            else:
                pass


    start=end+1
    end=(start-1)+len(str(seqs[list(seqs.keys())[0]][refname]))
    partnumber=partnumber+1

    splitfile=os.path.split(os.path.abspath(fname))
    partfile=str(splitfile[0]+'/'+splitfile[1].split(".")[0]+".partitions")
    partinfo=str(splitfile[0]+'/'+splitfile[1].split(".")[0]+".info")

    if os.path.isfile(partfile):
        with open(partfile, 'a+') as file:
            file.write(str("DNA, part"+str(partnumber)+" = "+str(start)+"-"+str(end)+"\n"))
    else :
        with open(partfile, 'w') as out:
            out.write(str("DNA, part"+str(partnumber)+" = "+str(start)+"-"+str(end)+"\n"))

    if os.path.isfile(partinfo):
        with open(partinfo, 'a+') as file:
            file.write(str(str(start)+"\t"+str(end)+"\t"+str(refname)+"\t"+"part"+str(partnumber)+"\n"))
    else :
        with open(partinfo, 'w') as out:
            out.write(str(str(start)+"\t"+str(end)+"\t"+str(refname)+"\t"+"part"+str(partnumber)+"\n"))

    for taxa in list(taxa_dict.keys()):
        if taxa in list(seqs.keys()):
            taxa_dict[taxa][refname]=seqs[taxa][refname]
        else:
            'we create a null sequence for missing taxa per gene'
            lenseq = len(str(seqs[list(seqs.keys())[0]][refname]))
            nullseq = "N"*lenseq
            taxa_dict[taxa][refname]=nullseq


'concatenation of sequences in dictionary and output sequences'
for taxa in taxa_dict:
    header= ">"+str(taxa)
    concatenate_seq=""
    for file in in_files:
        # if args.trim:
        #     refname=os.path.basename(file).replace(str("."+"trim"),"")
        # else:
        #     refname=os.path.basename(file).replace(str("."+"fa"),"")
        #
        refname=os.path.basename(file).replace(str("."+"fa"),"")
        concatenate_seq=concatenate_seq+str(taxa_dict[taxa][refname])

    Ncount=int(concatenate_seq.count("N"))
    gappcount=int(concatenate_seq.count("-"))
    SeqLength=len(concatenate_seq)
    Nratio=float(float(Ncount+gappcount)/float(SeqLength))
    AmbCount=Ncount+gappcount
    Ambratio=float(float(AmbCount)/float(SeqLength))
    'we add condition to check if a taxa is completely missing'
    seqtestN='N'*len(concatenate_seq)
    splitfile=os.path.split(os.path.abspath(fname))
    missfile=str(splitfile[0]+'/'+splitfile[1].split(".")[0]+".missingdata")
    with open(missfile, 'a+') as file:
        file.write(str(str(taxa)+"\t"+str(Ambratio)+"\n"))
    if seqtestN == concatenate_seq:
        print ("WARN: %s is missing" % (taxa))
    elif Nratio > Nthreshold:
        print ("WARN: %s has too missing data, not passed filter fraction of %s" % (taxa,Nthreshold))
        if os.path.isfile("missing_seqs.fasta"):
            'we stored only the missing sequences in a file'
            with open("missing_seqs.fasta", 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not taxa in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(concatenate_seq)+'\n')
                else:
                    pass
        else :
            with open("missing_seqs.fasta", 'w') as out:
                out.write(header+'\n')
                out.write(str(concatenate_seq)+'\n')
    else:
        if os.path.isfile(fname):
            with open(fname, 'a+') as file:
                old_headers = []
                end_file=file.tell()
                file.seek(0)
                for line in file:
                    if line.startswith(">"):
                        old_headers.append(line.replace(">","").split(";")[0])
                if not taxa in old_headers:
                    file.seek(end_file)
                    file.write(header+'\n')
                    file.write(str(concatenate_seq)+'\n')
                else:
                    pass
        else :
            with open(fname, 'w') as out:
                out.write(header+'\n')
                out.write(str(concatenate_seq)+'\n')
