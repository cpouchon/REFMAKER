#!/usr/bin/env python

import sys
import os, errno
import argparse
import time
from scipy import stats
import numpy as np

start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering consensus sequences from catalog alignments. Script was writen by C. Pouchon (2021).')
parser.add_argument("-c", help="minimal reference overlapp required to process a locus",
                    type=float)
parser.add_argument("-d","--cov", help="path to coverage and depth files",
                    type=str)
parser.add_argument("--df", help="filtering sample depth frequency",
type=float)
parser.add_argument("-f","--fasta", help="path to fasta consensus",
                    type=str)
parser.add_argument("--hf", help="filtering sample heterozygosity frequency",
type=float)
parser.add_argument("-l", help="minimal sequence length",
                    type=int)
parser.add_argument("-M", help="maximal missing loci allowed per sample",
                    type=float)
parser.add_argument("-m", help="maximal missing data allowed per sample",
                    type=float)
parser.add_argument("-m_out", help="maximal missing data allowed for outgroup samples",
                    type=float)
parser.add_argument("-o","--output_dir", help="output directory",
                    type=str)
parser.add_argument("-p","--pop", help="population file",
                    type=str)
parser.add_argument("-r", help="minimum percentage of individuals in a population required to process a locus",
                    type=float)
parser.add_argument("-R", help="minimum percentage of populations required to process a locus",
                    type=float)
parser.add_argument("--vcf", help="merge VCF file",
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



inpath=args.fasta
covpath=args.cov
outpath=args.output_dir
within_clust_pop=args.r #0.5
between_clust_pop=args.R #0.5
popfile=args.pop
moverlap=args.c #0.25
mmissing=args.m #0.8
mout_missing=args.m_out
lmissing=args.M #0.0
mlen=args.l #100
vcf=args.vcf
f_freq=args.df
h_freq=args.hf

mkdir(outpath)
logname="concatenated.log"
open(os.path.join(outpath, logname), 'w').close()

print("[INFOS]: selection of consensus sequences")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: selection of consensus sequences")+"\n")
print("[INFOS]: params set:")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: params set:")+"\n")
print("min. frequency of overlapping contig within population per locus: %s" % (str(within_clust_pop)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("min. frequency of overlapping contig within population per locus: %s" % (str(within_clust_pop)))+"\n")
print("min. frequency of populations sharing a locus: %s" % (str(between_clust_pop)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("min. frequency of populations sharing a locus: %s" % (str(between_clust_pop)))+"\n")
print("min. reference overlapp: %s" % (str(moverlap)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("min. reference overlapp: %s" % (str(moverlap)))+"\n")
print("max. frequency of missing data: %s" % (str(mmissing)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("max. frequency of missing data: %s" % (str(mmissing)))+"\n")
#print("max. frequency of missing loci: %s" % (str(lmissing)))


fname="concatenated.fa"
open(os.path.join(outpath, fname), 'w').close()

old_settings = np.seterr(all='ignore')

populations_inds={}
sampling=list()
with open(popfile) as p:
    poptab = p.readlines()
for line in poptab:
    l=line.rstrip()
    ind=l.split("\t")[0]
    sampling.append(ind)
    pop=l.split("\t")[1]
    if pop in list(populations_inds.keys()):
        populations_inds[pop].append(ind)
    else:
        populations_inds[pop]=[]
        populations_inds[pop].append(ind)


in_files = []
path_to_seek=os.walk(os.path.join(inpath))
for r, d, f in path_to_seek:
    for file in f:
        if file.endswith(".fa"):
            in_files.append(os.path.join(r, file))


#1) parse all files
clean_catalog={}
ref_list=list()
ref_len={}
stat_catalog={}
#ref_pop={}

toprint=str("[INFOS]: parsing consensus fasta files")
print("[INFOS]: parsing consensus fasta files")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(toprint+"\n")

for file in in_files:
    sampleID=os.path.basename(file).replace(str("."+"fa"),"")
    clean_catalog[sampleID]={}
    seqs={}
    cond_init=0
    count=0
    if sampleID in sampling:
        with open(file) as f:
            reftab = f.readlines()
        for line in reftab:
            if line.startswith(">"):
                l=line.rstrip()
                if cond_init==0:
                    cond_init=cond_init+1
                else:
                    if seqID not in list(seqs.keys()):
                        cons="".join(s)
                        seq_overlap=1-(cons.count("N")/len(cons))
                        seq_mapp=len(cons)-cons.count("N")
                        ref_len[seqID]=len(cons)
                        if seq_mapp>=mlen:
                            if float(seq_overlap)>=float(moverlap):
                                cons_f=cons
                                seqs[seqID]=cons_f
                                count=count+1
                                #for pop in populations_inds:
                                #    if sampleID in populations_inds[pop]:
                                #        old_count=ref_pop[seqID][pop]
                                #        new_count=old_count+1
                                #        ref_pop[seqID][pop]=new_count
                            else:
                                cons_f="N"*len(cons)
                                seqs[seqID]=cons_f
                        else:
                            cons_f="N"*len(cons)
                            seqs[seqID]=cons_f
                    else:
                        pass
                    cond_init=cond_init+1
                seqID=l.replace(">","")
                s=[]
                if seqID in ref_list:
                    pass
                else:
                    ref_list.append(seqID)
                    ref_len[seqID]=0
                    #ref_pop[seqID]={}
                    #for pop in populations_inds:
                    #    ref_pop[seqID][pop]=0
            else:
                l=line.rstrip()
                s.append(l)
                if line==reftab[len(reftab)-1]:
                    if seqID not in list(seqs.keys()):
                        cons="".join(s)
                        seq_overlap=1-(cons.count("N")/len(cons))
                        seq_mapp=len(cons)-cons.count("N")
                        ref_len[seqID]=len(cons)
                        if seq_mapp>=mlen:
                            if float(seq_overlap)>=float(moverlap):
                                cons_f=cons
                                seqs[seqID]=cons_f
                                count=count+1
                                #for pop in populations_inds:
                                #    if sampleID in populations_inds[pop]:
                                #        old_count=ref_pop[seqID][pop]
                                #        new_count=old_count+1
                                #        ref_pop[seqID][pop]=new_count
                            else:
                                cons_f="N"*len(cons)
                                seqs[seqID]=cons_f
                        else:
                            cons_f="N"*len(cons)
                            seqs[seqID]=cons_f
                    else:
                        pass
                else:
                    pass
        clean_catalog[sampleID]=seqs
        stat_catalog[sampleID]={}
        stat_catalog[sampleID]["count"]=count
        stat_catalog[sampleID]["fr"]=round(count/len(ref_list),2)
    else:
        continue
    #infos=list()
    #infos.append(sampleID)
    #infos.append(count)
    #infos.append(round(count/len(ref_list),2))
    #stat_catalog[sampleID]=infos

# pop without missing samples

#2) removing samples with high missing loci rate

filt_pop={}
bad_samples=list()
for pop in populations_inds:
    filt_pop[pop]=0
print("[INFOS]: locus coverage per sample")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: locus coverage per sample")+"\n")
print(str("sample")+"\t"+str("c.loci")+"\t"+str("f.loci"))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("sample")+"\t"+str("c.loci")+"\t"+str("f.loci")+"\n")
for samp in stat_catalog:
    for pop in populations_inds:
        if samp in populations_inds[pop]:
            if stat_catalog[samp]["fr"]>=(lmissing):
                old_count=filt_pop[pop]
                new_count=old_count+1
                filt_pop[pop]=new_count
            else:
                bad_samples.append(samp)
            print(str(samp)+"\t"+str(stat_catalog[samp]["count"])+"\t"+str(stat_catalog[samp]["fr"]))
            with open(os.path.join(outpath,logname), 'a+') as file:
                file.write(str(samp)+"\t"+str(stat_catalog[samp]["count"])+"\t"+str(stat_catalog[samp]["fr"])+"\n")

#print("[INFOS]: %s/%s samples removed according to the missing loci threshold at %s" % (len(bad_samples),len(stat_catalog),str(lmissing)))

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
            if "NODE" in l[0]:
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

bad_depth_contigs=list()
for cont in stat_depth:
    bad_ratio=float(stat_depth[cont].count("B")/(stat_depth[cont].count("B")+stat_depth[cont].count("G")))
    if bad_ratio<f_freq:
        pass
    else:
        bad_depth_contigs.append(cont)

print("[INFOS]: %s/%s loci removed according to the depth cutoff" % (len(bad_depth_contigs),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: %s/%s loci removed according to the depth cutoff" % (len(bad_depth_contigs),len(ref_list)))+"\n")

bads=bad_depth_contigs

# # selection on heterozygosity level
# print("[INFOS]: checking for heterozygous loci outliers")
# stat_ht={}
# dict = {}  ## {contig_id : [num_SNP ayant le moins de géno manquants, nb_missgeno]}
# filin = open(vcf, 'r')
# while True:
#     line = filin.readline()
#     if not line: break
#     if line.startswith("#"): continue
#     if not line.startswith("NODE"): continue
#     num_cluster = line.split("\t")[0]
#     if num_cluster not in stat_ht:
#         stat_ht[num_cluster]={}
#         stat_ht[num_cluster]["count"]=0
#         #stat_ht[num_cluster]["ratio"]=0.0
#     cond_SNP=line.split("\t")[4]
#     if cond_SNP == ".": continue
#     l=line.rstrip()
#     genotypes = [i.split(":")[0] for i in l.split("\t")[10:]]
#
#     if genotypes.count("0/1")<1: continue
#     oldc=stat_ht[num_cluster]["count"]
#     newc=oldc+1
#     stat_ht[num_cluster]["count"]=newc
# filin.close()


# selection on heterozygosity level - using pop frequency eg. X% of samples sharing an ht sites counted
print("[INFOS]: checking for outliers in heterozygous loci")
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: checking for outliers in heterozygous loci")+"\n")
stat_ht={}
dict = {}  ## {contig_id : [num_SNP ayant le moins de géno manquants, nb_missgeno]}
filin = open(vcf, 'r')
while True:
    line = filin.readline()
    if not line: break
    #if line.startswith("#CHROM"):
    #    l=line.rstrip()
    #    samp=l.split("\t")[9:]
    if line.startswith("#"): continue
    if not line.startswith("NODE"): continue
    num_cluster = line.split("\t")[0]
    if num_cluster not in stat_ht:
        stat_ht[num_cluster]={}
        stat_ht[num_cluster]["count"]=0
        #stat_ht[num_cluster]["ratio"]=0.0
    cond_SNP=line.split("\t")[4]
    if cond_SNP == ".": continue
    l=line.rstrip()
    genotypes = [i.split(":")[0] for i in l.split("\t")[10:]]

    if (genotypes.count("0/1")/len(genotypes))<h_freq: continue
    oldc=stat_ht[num_cluster]["count"]
    newc=oldc+1
    stat_ht[num_cluster]["count"]=newc
filin.close()

# calc het. ratio
ht_ratios=list()
ht_ids=list()
for num_cluster in stat_ht:
    len_cluster = int(num_cluster.split("_length_")[1].split("_cov")[0])
    covlen = int((max(stat_cov[num_cluster])/100)*len_cluster)
    ht_ratio=round(float(stat_ht[num_cluster]["count"]/len_cluster),3)
    ht_ratios.append(ht_ratio)
    ht_ids.append(num_cluster)
# perform Zscores for ratios
zh = np.abs(stats.zscore(ht_ratios))
threshold = 3
# Position of the outlier
bad_ht_contigs=list()
tmp=np.where(zh > 3)
tmp_id=tmp[0].tolist()
for i in tmp_id:
    cont=ht_ids[i]
    bad_ht_contigs.append(cont)

com_dh=0
com_list=list()
for c in bad_ht_contigs:
    if c in bads:
        com_dh=com_dh+1
        com_list.append(c)
    else:
        bads.append(c)

print("[INFOS]: %s/%s additional loci removed according to the heterozygosity cutoff" % (int(len(bad_ht_contigs)-com_dh),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: %s/%s additional loci removed according to the heterozygosity cutoff" % (int(len(bad_ht_contigs)-com_dh),len(ref_list)))+"\n")

# selection of contigs (ref) according to the population level for the final outputs
# add condition to break if filt_pop empty

cond_outgrp=0
if "OUT" in list(filt_pop.keys()):
    cond_outgrp=cond_outgrp+1
    good_outgp_shared={}
else:
    pass

bad_contigs=list()
ref_pop={}
for contig in ref_list:
    ref_pop[contig]={}
    for pop in populations_inds:
        ref_pop[contig][pop]=0

    for samp in sampling:
        if samp in bad_samples:
            pass
        else:
            seq=clean_catalog[samp][contig]
            if seq.count("N")==len(seq):
                pass
            else:
                for pop in populations_inds:
                    if samp in populations_inds[pop]:
                        old_count=ref_pop[contig][pop]
                        new_count=old_count+1
                        ref_pop[contig][pop]=new_count
                    else:
                        pass
    # count pop
    cond_r=0
    for pop in ref_pop[contig]:
        if ref_pop[contig][pop]/filt_pop[pop]>=within_clust_pop:
            cond_r=cond_r+1
        else:
            pass
    if cond_outgrp==0:
        if cond_r/len(filt_pop.keys())>=between_clust_pop:
            pass
            # infos si outgrp
        else:
            bad_contigs.append(contig)
            if contig in bads:
                pass
            else:
                bads.append(contig)
    else:
        num_pops=len(filt_pop.keys())-1
        if cond_r/num_pops>=between_clust_pop:
            if ref_pop[contig]["OUT"]>=1:
                good_outgp_shared[contig]=ref_pop[contig]["OUT"]
            else:
                pass
        else:
            bad_contigs.append(contig)
            if contig in bads:
                pass
            else:
                bads.append(contig)

print("[INFOS]: %s/%s loci removed according to the population level thresholds" % (len(bad_contigs),len(ref_list)))
with open(os.path.join(outpath,logname), 'a+') as file:
    file.write(str("[INFOS]: %s/%s loci removed according to the population level thresholds" % (len(bad_contigs),len(ref_list)))+"\n")

if cond_outgrp > 0:
    vals=set(list(good_outgp_shared.values()))
    print("[INFOS]: %s/%s loci shared with a least one outgroup" % (len(good_outgp_shared),len(ref_list)-len(bad_contigs)))
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[INFOS]: %s/%s loci shared with a least one outgroup" % (len(good_outgp_shared),len(ref_list)-len(bad_contigs)))+"\n")
    for v in vals:
        lc=list(good_outgp_shared.values())
        count=lc.count(v)
        print("%s/%s loci shared with %s outgroup taxa" % (count,len(ref_list)-len(bad_contigs),v))
        with open(os.path.join(outpath,logname), 'a+') as file:
            file.write(str("%s/%s loci shared with %s outgroup taxa" % (count,len(ref_list)-len(bad_contigs),v))+"\n")

# 4) concat proceed concat

## partition file
end=0
partnumber=0
partfile="concatenated.partitions"
partinfo="concatenated.info"
open(os.path.join(outpath, partfile), 'w').close()
open(os.path.join(outpath, partinfo), 'w').close()

for contig in ref_list:
    if contig in bad_contigs:
        pass
    else:
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
if (len(clean_catalog.keys())-len(bad_samples))==0:
    condsamp=condsamp+1
    print("[WARN]: all samples were removed according to the missing loci threshold. Try with an other value.")
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[WARN]: all samples were removed according to the missing loci threshold. Try with an other value.")+"\n")
else:
    pass
if (len(ref_list)-len(bad_contigs))==0:
    condcont=condcont+1
    print("[WARN]: all loci were removed according to the population thresholds. Try with other values.")
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[WARN]: all loci were removed according to the population thresholds. Try with other values.")+"\n")
else:
    pass

if condsamp==0 and condcont==0:
    print("[INFOS]: computing final output fasta files")
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[INFOS]: computing final output fasta files")+"\n")
    for samp in sampling:
        if samp in bad_samples:
            pass
        else:
            header= ">"+str(samp)
            concatenate_seq=""
            for contig in ref_list:
                if contig in bads:
                    continue
                else:
                    concatenate_seq=concatenate_seq+str(clean_catalog[samp][contig])
            flen=len(concatenate_seq)
            Ncount=int(concatenate_seq.count("N"))
            Nratio=float(float(Ncount)/float(len(concatenate_seq)))
            # write a condition for allowed missing data in outgp
            if samp in populations_inds["OUT"]:
                if Nratio<mout_missing:
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
                if Nratio<mmissing:
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
        print("[INFOS]: %s/%s taxa removed according to the final missing data threshold" % (len(missing_data_samples),len(sampling)))
        with open(os.path.join(outpath,logname), 'a+') as file:
            file.write(str("[INFOS]: %s/%s taxa removed according to the final missing data threshold" % (len(missing_data_samples),len(sampling)))+"\n")
    else:
        pass
    fsamp=len(sampling)-(len(missing_data_samples)+len(bad_samples))
    print("[INFOS]: final loci number: %s" % (len(ref_list)-len(bads)))
    print("[INFOS]: final samples: %s" % (fsamp))
    print("[INFOS]: final matrix length: %s" % (flen))
    with open(os.path.join(outpath,logname), 'a+') as file:
        file.write(str("[INFOS]: final loci number: %s" % (len(ref_list)-len(bads)))+"\n")
        file.write(str("[INFOS]: final samples: %s" % (fsamp))+"\n")
        file.write(str("[INFOS]: final matrix length: %s" % (flen))+"\n")
        file.close()
else:
    pass
