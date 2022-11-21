#!/usr/bin/env python

import sys
import os, errno
import argparse
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='Filtering VCF')
parser.add_argument("-i","--input", help="input VCF file",
                    type=str)
parser.add_argument("-o","--out", help="output VCF file",
                    type=str)
parser.add_argument("-m","--missing", help="max missing genotype to keep a variant",
                    type=float,default=1.0)
parser.add_argument("--maf", help="min minor allele frequency",
                    type=float,default=0.0)
parser.add_argument("-d","--depth", help="min minor allele depth to keep a variant (per sample)",
                    type=int,default=0)
parser.add_argument("-p","--ploidy", help="ploidy set for the calling",
                    type=int,default=2)
parser.add_argument("-t","--type", help="type of variant to keep",
                    type=str,default="all",choices=["snp","indel","all","snp+indel"])
parser.add_argument("--biallelic", help="keep only biallelic variants",
                    action="store_true")
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

# Default parameters
output_file=args.out
input_file=args.input
min_cov=args.depth
max_missing_prop=args.missing
min_maf=args.maf
#mapping_qual=args.mapq
ploidy=args.ploidy
query_type=args.type

keep_biallelic=0
if args.biallelic: keep_biallelic += 1

#filin = open(input_file, 'r')

print("Filtering VCF")
print("# filter parameters: indiv_AD>="+str(min_cov)+", missing<="+str(max_missing_prop)+", variant_type="+str(query_type))
print("# parsing VCF input file")
with open(input_file) as p:
    lines = p.readlines()
print("# done")

#filout = open(output_file, 'w')
open(output_file, 'w').close()
op=os.path.dirname(output_file)
oname=os.path.basename(output_file).split(".")[0]
oend=os.path.basename(output_file).split(".")[1]
filout_var=open(os.path.join(op, str(oname+"_var."+oend)), 'w')
open(os.path.join(op, str(oname+"_var."+oend)), 'w').close()

n_geno = 0
max_missing = 0 # integer
var_count = 0
var_written = 0
na_count = 0
new_na_count = 0
selected_na_count = 0
indel_count=0
snp_bi_count=0
snp_multi_count=0

missformat=list()
for i in range(0,ploidy):
    missformat.append(".")
missgeno="/".join(missformat)

Bar = ProgressBar(len(lines), 60, '\t filtering VCF')
barp=0
#print("# filtering VCF")
for line in lines:
    barp=barp+1
    linest=line.rstrip()
    if line.startswith("#"):
        with open(output_file, 'a+') as filout:
            filout.write(line)
        with open(os.path.join(op, str(oname+"_var."+oend)), 'a+') as filout_var:
            filout_var.write(line)
        #filout.write(line)
        #filout_var.write(line)
        #print(line.rstrip())
        if line.startswith("#CHROM"):  # Used to record the total number of samples
            n_geno=len(linest.split("\t")) - 9
            max_missing = int(max_missing_prop*n_geno)
        continue
        Bar.update(barp)
    else:
        splitted=linest.split("\t")
        ALTinfo=splitted[4]
        if ALTinfo==".":
            #filout.write(linest)
            with open(output_file, 'a+') as filout:
                filout.write(linest+'\n')
            #print(linest.rstrip())
            Bar.update(barp)
        else:
            var_count += 1
            #define the type of variant (snp/other)
            REFinfo=splitted[3]
            ALT=ALTinfo.split(",")
            num_alleles=1
            complex=0
            if len(REFinfo)>1:
                # cond 1: reference has multiple nt
                for i in ALT:
                    if len(i)>1:
                        complex += 1
                    else:
                        num_alleles += 1
                if len(ALT)==1 and complex==0:
                    var_type="indel"
                    indel_count += 1
                else:
                    var_type="complex"
            else:
                # cond 2: reference has a single nt
                for i in ALT:
                    if len(i)>1:
                        complex += 1
                    else:
                        num_alleles += 1
                if len(ALT)==1 and complex==0:
                    var_type="snp"
                    snp_bi_count += 1
                elif len(ALT)==1 and complex>0:
                    var_type="indel"
                    indel_count += 1
                elif len(ALT)>1 and complex==0:
                    var_type="snp"
                    snp_multi_count += 1
                else:
                    var_type="complex"
            #conditions to keep the querried variant type
            if query_type!="all":
                if query_type=="indel" and var_type!="indel":
                    continue
                elif query_type=="snp+indel" and var_type=="complex":
                    continue
                elif query_type=="snp" and var_type!="snp":
                    continue
            if keep_biallelic==1 and num_alleles>2:
                    continue

            # condition for MAPPING_QUAL
            # condition for missing genotype
            # condition for depth
            start_geno = 9
            line_towrite = "\t".join(splitted[:start_geno])
            missing_count = 0
            ref_count = 0
            alt_count = 0
            FORMAT={}
            for f_part in splitted[7].split(";"):
                if "=" in f_part:
                    k=f_part.split("=")[0]
                    v=f_part.split("=")[1]
                    FORMAT[k]=v
                else:
                    continue
            INFOS=splitted[8].split(":")
            AD_INFOS=INFOS.index("AD")

            #if "MQM" in list(FORMAT.keys()):
            #    if float(FORMAT["MQM"]) < mapping_qual: continue
            #elif "MQ" in list(FORMAT.keys()):
            #    if float(FORMAT["MQ"]) < mapping_qual: continue

            for geno in splitted[start_geno:]:
                geno_info = geno.split(":")
                genotype = geno_info[0]
                geno_ad=geno_info[AD_INFOS]
                if genotype==missgeno:
                    missing_count += 1
                    na_count += 1
                    geno_new=genotype+":"+":".join(geno_info[1:])
                else:
                    tmp_ad=geno_ad.split(",")
                    tmp_geno=genotype.split("/")
                    cond_ad=0
                    # condition ad for alt allele(s) -> consider as missing if not filled
                    for i in range(1,len(tmp_ad)):
                        if tmp_ad[i] != ".":
                            if int(tmp_ad[i]) < min_cov:
                                cond_ad += 1
                    if cond_ad>0:
                        genotype=missgeno
                        missing_count += 1
                        na_count += 1
                        tmp=list()
                        for ele in geno_info[1:]:
                            if '\n' in ele:
                                tmp.append(".\n")
                            else:
                                tmp.append(".")
                        new_info=":".join(tmp)
                        geno_new = genotype+":"+":".join(tmp)

                    else:
                        geno_new = genotype+":"+":".join(geno_info[1:])
                        #count ref/alt alleles
                        for a in tmp_geno:
                            if a=="0":
                                ref_count += 1
                            else:
                                alt_count += 1
                #rewrite geno info
                line_towrite += "\t"+geno_new

            maf = 0
            if ref_count+alt_count>0:
                if alt_count>0:
                    Bar.update(barp)
                    maf = (float) (min(ref_count,alt_count))/(ref_count+alt_count)
                    ## Filter on missing count and maf
                    if missing_count <= max_missing and maf >= min_maf:
                        with open(output_file, 'a+') as filout:
                            filout.write(line_towrite+'\n')
                        with open(os.path.join(op, str(oname+"_var."+oend)), 'a+') as filout_var:
                            filout_var.write(line_towrite+'\n')
                        #filout.write(line_towrite)
                        #filout_var.write(line_towrite)
                        var_written += 1
                else:
                    Bar.update(barp)
                    split_line_towrite=line_towrite.split("\t")
                    split_line_towrite[4]="."
                    new_line_towrite="\t".join(split_line_towrite)
                    with open(output_file, 'a+') as filout:
                        filout.write(new_line_towrite+'\n')
                    #filout.write(new_line_towrite)

#filin.close()
#filout.close()
#filout_var.close()
print("")
print("# total variants: "+str(var_count))
print("  indels: "+str(indel_count)+", snps: "+str((snp_multi_count+snp_bi_count))+" (biallelic: "+str(snp_bi_count)+")")
print("# filtered variants: "+str(var_written))
print("done")
