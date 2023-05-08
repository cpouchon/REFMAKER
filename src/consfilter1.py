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

parser = argparse.ArgumentParser(description='Filtering consensus sequences from catalog alignments. Script was writen by C. Pouchon (2021).')
parser.add_argument("-d","--cov", help="path to coverage and depth files",
                    type=str)
parser.add_argument("--df", help="filtering sample depth frequency (eg. 0.5)",
                    type=float)
parser.add_argument("-f","--fasta", help="path to fasta consensus",
                    type=str)
parser.add_argument("--ht",help="maximal heterozygous sites frequency (per sample)",
                    type=float)
parser.add_argument("--Ht",help="maximal shared heterozygous site (accross samples)",
                    type=float)
parser.add_argument("-l", help="minimal sequence length",
                    type=int)
parser.add_argument("-M", help="final missing data allowed per sample (across loci)",
                    type=float)
parser.add_argument("-m", help="maximal missing data allowed (locus)",
                    type=float)
parser.add_argument("-m_out", help="maximal missing data allowed for outgroup samples (locus)",
                    type=float)
parser.add_argument("-M_out", help="final missing data allowed allowed for outgroup samples (across loci)",
                    type=float)
parser.add_argument("-o","--output_dir", help="output directory",
                    type=str)
parser.add_argument("-p","--pop", help="population file",
                    type=str)
parser.add_argument("--remove", help="removing of paralogous sites/loci",
                    type=str, choices=["sites","loci"])
parser.add_argument("-r", help="minimum percentage of individuals in a population required to process a locus",
                    type=float)
parser.add_argument("-R", help="minimum percentage of populations required to process a locus",
                    type=float)
parser.add_argument("--window_size", help="sliding window size to identify polymorphic sites (eg. 20bp)",
                    type=int)
parser.add_argument("--window_psites", help="maxmimal polymorphic sites allowed within sliding window",
                    type=int)

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



inpath=args.fasta
covpath=args.cov
outpath=args.output_dir
within_clust_pop=args.r #0.5
between_clust_pop=args.R #0.5
popfile=args.pop
seq_h=args.ht #0.05
samples_h=args.Ht #0.5
mmissing=args.m #0.8
mmissing_out=args.m_out
lmissing_out=args.M_out
lmissing=args.M #0.0
w_site=args.window_psites #5
w_size=args.window_size #20
cond_remove_paralog=args.remove #loci
f_freq=args.df #0.5 (covs chez 50% des indivs) ou 0.25
mlen=args.l #100


mkdir(outpath)
logname="concatenated.log"
open(os.path.join(outpath, logname), 'w').close()

print("[INFOS]: Filtering of consensus sequences")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: Filtering of consensus sequences")+"\n")
print("[INFOS]: params set:")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: params set:")+"\n")
print("min. frequency of overlapping contig within population per locus: %s" % (str(within_clust_pop)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("min. frequency of overlapping contig within population per locus: %s" % (str(within_clust_pop)))+"\n")
print("min. frequency of populations sharing a locus: %s" % (str(between_clust_pop)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("min. frequency of populations sharing a locus: %s" % (str(between_clust_pop)))+"\n")
print("max. frequency of missing data per sample per locus: %s" % (str(mmissing)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("max. frequency of missing data per sample per locus: %s" % (str(mmissing)))+"\n")
print("final frequency of missing data per sample: %s" % (str(lmissing)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("final frequency of missing data per sample: %s" % (str(lmissing)))+"\n")
#print("max. frequency of missing loci: %s" % (str(lmissing)))
print("max. frequency of missing data per sample per locus for outgroups: %s" % (str(mmissing_out)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("max. frequency of missing data per sample per locus for outgroups: %s" % (str(mmissing_out)))+"\n")
print("final frequency of missing data per sample for outgroups: %s" % (str(lmissing_out)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("final frequency of missing data per sample for outgroups: %s" % (str(lmissing_out)))+"\n")
print("max. heterozygosity frequency per sequence: %s" % (str(seq_h)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("max. heterozygosity frequency per sequence: %s" % (str(seq_h)))+"\n")
print("max. frequency of shared heterozygote site per locus: %s" % (str(samples_h)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("max. heterozygosity frequency par locus: %s" % (str(samples_h)))+"\n")
print("sliding window size to detect polymorphic sites in sequence: %s" % (str(w_size)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("sliding window size to detect polymorphic sites in sequence: %s" % (str(w_size)))+"\n")
print("max. of polymorphic sites within the sliding window: %s" % (str(w_site)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("max. of polymorphic sites within the sliding window: %s" % (str(w_site)))+"\n")
print("removing of paralogs by: %s" % (str(cond_remove_paralog)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("removing of paralogs by: %s" % (str(cond_remove_paralog)))+"\n")
print("max. frequency of samples sharing loci with excess coverage depth: %s" % (str(f_freq)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("max. frequency of samples sharing loci with excess coverage depth: %s" % (str(f_freq)))+"\n")
print("min. sequence length of loci: %s" % (str(mlen)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("min. sequence length of loci: %s" % (str(mlen)))+"\n")


fname="concatenated.fa"
open(os.path.join(outpath, fname), 'w').close()

old_settings = np.seterr(all='ignore')

populations_inds={}
sampling=list()
stats_ind={}
with open(popfile) as p:
    poptab = p.readlines()
for line in poptab:
    l=line.rstrip()
    ind=l.split("\t")[0]
    sampling.append(ind)
    if ind in list(stats_ind.keys()):
        pass
    else:
        stats_ind[ind]={}
        stats_ind[ind]["c.missing"]=0
        stats_ind[ind]["c.ht"]=0
        stats_ind[ind]["c.paralogs"]=0
    pop=l.split("\t")[1]
    if pop in list(populations_inds.keys()):
        populations_inds[pop].append(ind)
    else:
        populations_inds[pop]=[]
        populations_inds[pop].append(ind)

cond_outgrp=0
if "OUT" in list(populations_inds.keys()):
    cond_outgrp=cond_outgrp+1
    good_outgp_shared={}
else:
    pass

in_files = []
path_to_seek=os.walk(os.path.join(inpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(".fa"):
            in_files.append(os.path.join(r, file))


#1) parse all files
clean_catalog_1={}
ref_list=list()
ref_len={}
stat_catalog={}
bad_loci_len=list()
bad_loci_pop=list()
bad_loci_ht=list()
ref_pop={}
#ref_pop={}
ht_basis=["R","Y","S","W","K","M","B","D","H","V","N"]
if cond_remove_paralog=="sites":
    removed_sites=0

toprint=str("")
print("")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(toprint+"\n")
toprint=str("[INFOS]: parsing consensus fasta files")
print("[INFOS]: parsing consensus fasta files")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(toprint+"\n")

Bar = ProgressBar(len(in_files), 60, '\t ')
barp=0

print("total number of loci: %s" % (str(len(in_files))))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("total number of loci: %s" % (str(len(in_files))))+"\n")

for file in in_files:
    barp=barp+1
    geneID=os.path.basename(file).replace(str("."+"fa"),"")
    seqs=to_dict_remove_dups(SeqIO.parse(file, "fasta"))
    ref_list.append(geneID)
    tmpdict={}
    # cond1: ind. missing data
    # cond2: ind. ht
    # cond3: polymorphic sites
    # cond4: global ht
    # bad sample ID stored in list
    bad=list()
    cond12_list=list()
    #cond3_dict=dict()
    cond3_list=list()
    cond4_list=list()
    refseq=str(seqs[geneID].seq)
    lenseq=len(refseq)
    ref_len[geneID]=lenseq
    bad_sites=list()
    bad_samp=list()
    if lenseq<mlen:
        bad_loci_len.append(geneID)
        for samp in sampling:
            bad.append(samp)
            stats_ind[samp]["c.missing"]+=1
        #continue
    else:
        for samp in sampling:
            if samp in list(seqs.keys()):
                sequence=str(seqs[samp].seq).upper()
                cond1=sequence.count("-")/lenseq #missing/gappy data
                cond2=sum([sequence.count(k) for k in ht_basis])/lenseq #ht bases count condition
                if cond_outgrp>0:
                    if samp in populations_inds["OUT"]:
                        if cond1>=mmissing_out:
                            bad.append(samp)
                            stats_ind[samp]["c.missing"]+=1 #voir si condition speciale outgroups
                    else:
                        if cond1>=mmissing:
                            bad.append(samp)
                            stats_ind[samp]["c.missing"]+=1 #voir si condition speciale outgroups
                else:
                    if cond1>=mmissing:
                        bad.append(samp)
                        stats_ind[samp]["c.missing"]+=1 #voir si condition speciale outgroups
                if cond2>=seq_h:
                    bad.append(samp)
                    stats_ind[samp]["c.ht"]+=1
                if samp in bad:
                    pass
                else:
                    tmpdict[samp]=sequence
                    #cond3_dict[samp]=0
            else:
                bad.append(samp)
                stats_ind[samp]["c.missing"]+=1

    if len(tmpdict)>0:
        bad_sites=list()
        bad_samp=list()
        sites=list()
        for i in range(lenseq - w_size + 1):
            bmin=i
            bmax=i+w_size
            poly_dict=dict.fromkeys(list(tmpdict.keys()), 0)
            for j in range(bmin,bmax):
                tmpsite=list()
                for s in tmpdict:
                    tmpsite.append(tmpdict[s][j])
                    # polymorphic sites - paralogs
                    if tmpdict[s][j]=="-":
                        continue
                    else:
                        if tmpdict[s][j]!=refseq[j]:
                            poly_dict[s]+=1
                        else:
                            continue
                if len(set(tmpsite))==1:
                    if '-' in tmpsite:
                        if j in bad_sites:
                            pass
                        else:
                            bad_sites.append(j) # only gappy sites
                #cond3 - shared ht site --> TO REMOVE?
                ht_bases=sum([tmpsite.count(k) for k in ht_basis])
                miss_bases=tmpsite.count("-")
                tot_bases=len(tmpsite)
                if tot_bases==miss_bases:
                    pass
                else:
                    cond3_ratio=ht_bases/(tot_bases-miss_bases)
                    if cond3_ratio>=samples_h:
                        if cond_remove_paralog=="loci":
                            if geneID in bad_loci_ht:
                                pass
                            else:
                                bad_loci_ht.append(geneID)
                        else:
                            bad_sites.append(j)
                            removed_sites=+1
            tmp_samp=[k for (k,v) in poly_dict.items() if v >= w_site]
            if len(tmp_samp)>0:
                for s in tmp_samp:
                    if s in bad_samp:
                        pass
                    else:
                        bad_samp.append(s)
                        stats_ind[s]["c.paralogs"]+=1
        Bar.update(barp)
    else:
        pass
        Bar.update(barp)

    if geneID in bad_loci_len or geneID in bad_loci_ht:
        pass
    else:
        if len(refseq)-len(bad_sites)<mlen:
            bad_loci_len.append(geneID)
        else:
            ref_len[geneID]=len(refseq)-len(bad_sites)
            clean_catalog_1[geneID]={}
            ref_pop[geneID]={}
            # condition pour la population (r and R)
            for pop in populations_inds:
                ref_pop[geneID][pop]=0
            for s in tmpdict:
                if s in list(clean_catalog_1[geneID].keys()):
                    pass
                else:
                    #save clean
                    sseq=tmpdict[s]
                    newseq="".join([char for idx, char in enumerate(sseq) if idx not in bad_sites]) #we remove bad sites (gappy + paralog if specified)
                    clean_catalog_1[geneID][s]=newseq
                for pop in populations_inds:
                    if s in populations_inds[pop]:
                        ref_pop[geneID][pop]+=1
                    else:
                        pass

            # count pop
            cond_r=0
            for pop in ref_pop[geneID]:
                if ref_pop[geneID][pop]/len(populations_inds[pop])>=within_clust_pop:
                    cond_r=cond_r+1
                else:
                    pass
            if cond_outgrp==0:
                if cond_r/len(populations_inds.keys())>=between_clust_pop:
                    pass
                    # infos si outgrp
                else:
                    bad_loci_pop.append(geneID)

            else:
                num_pops=len(populations_inds.keys())-1
                if cond_r/num_pops>=between_clust_pop:
                    if ref_pop[geneID]["OUT"]>=1:
                        good_outgp_shared[geneID]=ref_pop[geneID]["OUT"]
                    else:
                        pass
                else:
                    bad_loci_pop.append(geneID)

print("")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("")+"\n")

# parse coverage files to filter out outliers contigs in depth

in_covs = []
path_to_seek=os.walk(os.path.join(covpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(".infos"):
            in_covs.append(os.path.join(r, file))

print("[INFOS]: checking for outliers in locus depths")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: checking for outliers in locus depths")+"\n")

stat_depth={}
stat_cov={}
for file in in_covs:
    sampleID=os.path.basename(file).replace(str("_cov.infos"),"")
    if sampleID in sampling:
        cov_dict={}
        depths=list()
        covs=list()
        ids=list()
        with open(file) as f:
            fcov = f.readlines()
        for line in fcov:
            l=line.rstrip().split("\t")
            if l[0] in ref_list:
                cont=l[0]
                cov=float(l[5])
                depth=float(l[6])
                cov_dict[cont]=depth
                depths.append(depth)
                covs.append(cov)
                ids.append(cont)
                if cont in stat_depth.keys():
                    pass
                else:
                    stat_depth[cont]=list()
                    stat_cov[cont]=list()
            else:
                continue
        z = np.abs(stats.zscore(depths))
        threshold = 3
        # Position of the outlier
        tmp=np.where(z > 3)
        bads_cov=tmp[0].tolist()

        for i in range(len(ids)):
            cid=ids[i]
            stat_cov[cid].append(covs[i])
            if i in bads_cov:
                stat_depth[cid].append(str("B"))
            else:
                stat_depth[cid].append(str("G"))
    else:
        continue

bad_depth_loci=list()
for cont in stat_depth:
    bad_ratio=float(stat_depth[cont].count("B")/(stat_depth[cont].count("B")+stat_depth[cont].count("G")))
    if bad_ratio<f_freq:
        pass
    else:
        bad_depth_loci.append(cont)


tot_loci=len(in_files)
for s in stats_ind:
    stats_ind[s]["f.missing"]=float(stats_ind[s]["c.missing"])/float(tot_loci)

## print stats inds
print("[INFOS]: statistics over filtering")
print("number of removed loci per sample:")
print("sample"+"\t"+"c.missing"+"\t"+"f.missing"+"\t"+"c.heterozygous"+"\t"+"c.paralogous")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: statistics over filtering"))
    file.write(str("number of removed loci per sample:")+"\n")
    file.write("sample"+"\t"+"c.missing"+"\t"+"f.missing"+"\t"+"c.heterozygous"+"\t"+"c.paralogous"+"\n")
for samp in stats_ind:
    print(samp+"\t"+str(stats_ind[samp]["c.missing"])+"\t"+str(round(stats_ind[samp]["f.missing"],2))+"\t"+str(stats_ind[samp]["c.ht"])+"\t"+str(stats_ind[samp]["c.paralogs"]))
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(samp+"\t"+str(stats_ind[samp]["c.missing"])+"\t"+str(round(stats_ind[samp]["f.missing"],2))+"\t"+str(stats_ind[samp]["c.ht"])+"\t"+str(stats_ind[samp]["c.paralogs"])+"\n")


print("%s/%s loci removed according to the depth cutoff" % (len(bad_depth_loci),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("%s/%s loci removed according to the depth cutoff" % (len(bad_depth_loci),len(ref_list)))+"\n")

bads=set(bad_loci_len+bad_loci_ht+bad_depth_loci+bad_loci_pop)

# condition population level
print("%s/%s loci removed according to the minimal length" % (len(bad_loci_len),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("%s/%s loci removed according to the minimal length" % (len(bad_loci_len),len(ref_list)))+"\n")

print("%s/%s loci removed according to the heterozygosity" % (len(bad_loci_ht),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("%s/%s loci removed according to the heterozygosity" % (len(bad_loci_ht),len(ref_list)))+"\n")

print("%s/%s loci removed according to the population level thresholds" % (len(bad_loci_pop),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("%s/%s loci removed according to the population level thresholds" % (len(bad_loci_pop),len(ref_list)))+"\n")

print("remaining loci: %s/%s" % (len(ref_list)-len(bads),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("remaining loci: %s/%s" % (len(ref_list)-len(bads),len(ref_list)))+"\n")

if cond_outgrp > 0:
    shared_notbads={k:v for (k,v) in good_outgp_shared.items() if k not in bads}
    vals=set(list(shared_notbads.values()))
    print("%s/%s loci shared with a least one outgroup" % (len(shared_notbads),len(ref_list)-len(bads)))
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("%s/%s loci shared with a least one outgroup" % (len(shared_notbads),len(ref_list)-len(bads)))+"\n")
    for v in vals:
        lc=list(shared_notbads.values())
        count=lc.count(v)
        print("%s/%s loci shared with %s outgroup taxa" % (count,len(ref_list)-len(bads),v))
        with open(os.path.join(outpath,logname), 'a+') as file:
            file.write(str("%s/%s loci shared with %s outgroup taxa" % (count,len(ref_list)-len(bads),v))+"\n")

# 4) concat proceed outputs

## partition file
end=0
partnumber=0
partfile="concatenated.partitions"
partinfo="concatenated.info"
open(os.path.join(outpath, partfile), 'w').close()
open(os.path.join(outpath, partinfo), 'w').close()

for contig in ref_list:
    if contig in bads:
        pass
    else:
        if contig in list(clean_catalog_1.keys()):
            start=end+1
            end=(start-1)+ref_len[contig]
            partnumber=partnumber+1
            if os.path.isfile(os.path.join(outpath,partfile)):
                with open(os.path.join(outpath,partfile), 'a+') as file:
                    file.write(str("DNA, part"+str(partnumber)+" = "+str(start)+"-"+str(end)+"\n"))
            else :
                with open(os.path.join(outpath,partfile), 'w') as out:
                    out.write(str("DNA, part"+str(partnumber)+" = "+str(start)+"-"+str(end)+"\n"))
            if os.path.isfile(os.path.join(outpath,partinfo)):
                with open(os.path.join(outpath,partinfo), 'a+') as file:
                    file.write(str(str(start)+"\t"+str(end)+"\t"+str(contig)+"\t"+"part"+str(partnumber)+"\n"))
            else :
                with open(os.path.join(outpath,partinfo), 'w') as out:
                    out.write(str(str(start)+"\t"+str(end)+"\t"+str(contig)+"\t"+"part"+str(partnumber)+"\n"))

## fasta file
missing_data_samples=list()
condsamp=0
condcont=0
if (len(ref_list)-len(bads))==0:
    condcont=condcont+1
    print("[WARN]: all loci were removed according to the parameters used. Please try with other values.")
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[WARN]: all loci were removed according to the parameters used. Please try with other values.")+"\n")
else:
    pass

# par échantillon dans sampling --> pour chaque contig de la liste si dans bads pass sinon si echantillon présent sequence/sinon mettre des gapp

if condsamp==0 and condcont==0:
    print("[INFOS]: computing final output fasta files")
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[INFOS]: computing final output fasta files")+"\n")
    for samp in sampling:
        header= ">"+str(samp)
        concatenate_seq=""
        for contig in ref_list:
            if contig in bads:
                continue
            else:
                if contig in list(clean_catalog_1.keys()):
                    if samp in list(clean_catalog_1[contig].keys()):
                        concatenate_seq=concatenate_seq+str(clean_catalog_1[contig][samp])
                    else:
                        len_contig=ref_len[contig]
                        concatenate_seq=concatenate_seq+str("-"*len_contig)
        flen=len(concatenate_seq)
        Ncount=int(concatenate_seq.count("-"))
        Nratio=float(float(Ncount)/float(len(concatenate_seq)))
        # write a condition for allowed missing data in outgp
        if cond_outgrp > 0:
            if samp in populations_inds["OUT"]:
                if Nratio<lmissing_out:
                    if os.path.isfile(os.path.join(outpath,fname)):
                        with open(os.path.join(outpath,fname), 'a+') as file:
                            old_headers = []
                            end_file=file.tell()
                            file.seek(0)
                            for line in file:
                                if line.startswith(">"):
                                    old_headers.append(line.replace(">",""))
                            if not samp in old_headers:
                                file.seek(end_file)
                                file.write(header+'\n')
                                file.write(str(concatenate_seq)+'\n')
                            else:
                                pass
                    else :
                        with open(os.path.join(outpath,fname), 'w') as out:
                            out.write(header+'\n')
                            out.write(str(concatenate_seq)+'\n')
                    if os.path.isfile(os.path.join(outpath,"concatenated.missing")):
                        with open(os.path.join(outpath,"concatenated.missing"), 'a+') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                    else:
                        with open(os.path.join(outpath,"concatenated.missing"), 'w') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                else:
                    missing_data_samples.append(samp)
                    if os.path.isfile(os.path.join(outpath,"concatenated.missing")):
                        with open(os.path.join(outpath,"concatenated.missing"), 'a+') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                    else:
                        with open(os.path.join(outpath,"concatenated.missing"), 'w') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
            else:
                if Nratio<lmissing:
                    if os.path.isfile(os.path.join(outpath,fname)):
                        with open(os.path.join(outpath,fname), 'a+') as file:
                            old_headers = []
                            end_file=file.tell()
                            file.seek(0)
                            for line in file:
                                if line.startswith(">"):
                                    old_headers.append(line.replace(">",""))
                            if not samp in old_headers:
                                file.seek(end_file)
                                file.write(header+'\n')
                                file.write(str(concatenate_seq)+'\n')
                            else:
                                pass
                    else :
                        with open(os.path.join(outpath,fname), 'w') as out:
                            out.write(header+'\n')
                            out.write(str(concatenate_seq)+'\n')
                    if os.path.isfile(os.path.join(outpath,"concatenated.missing")):
                        with open(os.path.join(outpath,"concatenated.missing"), 'a+') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                    else:
                        with open(os.path.join(outpath,"concatenated.missing"), 'w') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                else:
                    missing_data_samples.append(samp)
                    if os.path.isfile(os.path.join(outpath,"concatenated.missing")):
                        with open(os.path.join(outpath,"concatenated.missing"), 'a+') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                    else:
                        with open(os.path.join(outpath,"concatenated.missing"), 'w') as file:
                            file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
        else:
            if Nratio<lmissing:
                if os.path.isfile(os.path.join(outpath,fname)):
                    with open(os.path.join(outpath,fname), 'a+') as file:
                        old_headers = []
                        end_file=file.tell()
                        file.seek(0)
                        for line in file:
                            if line.startswith(">"):
                                old_headers.append(line.replace(">",""))
                        if not samp in old_headers:
                            file.seek(end_file)
                            file.write(header+'\n')
                            file.write(str(concatenate_seq)+'\n')
                        else:
                            pass
                else :
                    with open(os.path.join(outpath,fname), 'w') as out:
                        out.write(header+'\n')
                        out.write(str(concatenate_seq)+'\n')
                if os.path.isfile(os.path.join(outpath,"concatenated.missing")):
                    with open(os.path.join(outpath,"concatenated.missing"), 'a+') as file:
                        file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                else:
                    with open(os.path.join(outpath,"concatenated.missing"), 'w') as file:
                        file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
            else:
                missing_data_samples.append(samp)
                if os.path.isfile(os.path.join(outpath,"concatenated.missing")):
                    with open(os.path.join(outpath,"concatenated.missing"), 'a+') as file:
                        file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")
                else:
                    with open(os.path.join(outpath,"concatenated.missing"), 'w') as file:
                        file.write(str(samp)+"\t"+str(round(Nratio,2))+"\n")


    if len(missing_data_samples)>0:
        print("%s/%s taxa removed according to the final missing data threshold" % (len(missing_data_samples),len(sampling)))
        with open(os.path.join(outpath,logname), 'a+') as file:
            file.write(str("%s/%s taxa removed according to the final missing data threshold" % (len(missing_data_samples),len(sampling)))+"\n")
    else:
        pass
    fsamp=len(sampling)-len(missing_data_samples)
    print("[INFOS]: final matrix")
    print("loci number: %s" % (len(ref_list)-len(bads)))
    print("samples: %s" % (fsamp))
    print("length (bp): %s" % (flen))
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[INFOS]: final matrix")+"\n")
        file.write(str("loci number: %s" % (len(ref_list)-len(bads)))+"\n")
        file.write(str("samples: %s" % (fsamp))+"\n")
        file.write(str("length (bp): %s" % (flen))+"\n")
        file.close()
else:
    pass
