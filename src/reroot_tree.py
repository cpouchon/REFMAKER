#!/usr/bin/env python
import sys
from ete3 import Tree
import os, errno
import argparse
import time


start_time = time.time()

parser = argparse.ArgumentParser(description='Reroot tree(s) according to outgroups taxa')
parser.add_argument("-i","--input", help="input tree(s)",
                    type=str)
parser.add_argument("--outgroups", help="list of outgroup taxa",
                    type=str)
parser.add_argument("--partitions", help="list of partitions",
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

intrees=args.input
outgroup_names = [str(item) for item in args.outgroups.split(',')]
names = args.partitions

#outgroup_names =["Grias_cauliflora_79571_SRR8599511","Gustavia_augusta_372736_SRR8599522","Gustavia_hexapetala_372738_TOU004262_BGN_FNM","Gustavia_hexapetala_372738_TOU010328_BGN_PFL","Gustavia_serrata_1421040_SRR8599531","Gustavia_sp_79572_513"]
#print(outgroup_names)
mkdir("./locus_trees")

listnames=list()
with open(names) as tab:
    for part in tab:
        l = part.rstrip()
        name = l.split(" ")[1]
        listnames.append(name)
lnum=0
with open(intrees) as trees:
    for line in trees:
        t = Tree(line.rstrip())
        #print(t)
        outgroups_in_tree = list(set(t.get_leaf_names()).intersection(set(outgroup_names)))
        #print(outgroups_in_tree)
        tname=listnames[lnum]
        #print(tname)
        lnum=lnum+1
        ingp=list(set(t.get_leaf_names()).difference(set(outgroup_names)))
        t.set_outgroup(ingp[0])
        if len(outgroups_in_tree)>1:
            ancestor = t.get_common_ancestor(outgroups_in_tree)
            t.set_outgroup(ancestor)
            t.write(format=1, outfile=os.path.join("./locus_trees/", str(str(tname)+".tre")))
        elif len(outgroups_in_tree) == 1:
            t.set_outgroup(outgroups_in_tree[0])
            t.write(format=1, outfile=os.path.join("./locus_trees/", str(str(tname)+".tre")))
        else:
            continue
