#!/usr/bin/env python

import glob
import argparse
import sys
import os, errno
import networkx as nx
import matplotlib.pyplot as plt
import markov_clustering as mc
import random
from scipy.sparse import csr_matrix, isspmatrix, issparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Filtering of reference metacontigs after blast')
parser.add_argument("-b","--blast", help="output from blast of contigs into database",
                    type=str)
parser.add_argument("-i","--infile", help="reference fasta",
                    type=str)
parser.add_argument("-l","--overlap", help="min. overlap",
                    type=float)
parser.add_argument("-s","--similarlity", help="min. similarlity",
                    type=float)
parser.add_argument("-o","--outpath", help='path to outpout results', type=str)


if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

try:
    from matplotlib.pylab import show, cm, axis
except ImportError:
    sys.stderr.write("Matplotlib not present\n")
    raise

def draw_graph2(matrix, clusters, **kwargs):
    """
    Visualize the clustering

    :param matrix: The unprocessed adjacency matrix
    :param clusters: list of tuples containing clusters as returned
                     by 'get_clusters'
    :param kwargs: Additional keyword arguments to be passed to
                   networkx.draw_networkx
    """
    # make a networkx graph from the adjacency matrix
    graph = nx.Graph(matrix)

    # map node to cluster id for colors
    cluster_map = {node: i for i, cluster in enumerate(clusters) for node in cluster}
    colors = [cluster_map[i] for i in range(len(graph.nodes()))]

    # if colormap not specified in kwargs, use a default
    if not kwargs.get("cmap", False):
        kwargs["cmap"] = cm.tab20

    # draw
    nx.draw_networkx(graph, node_color=colors, **kwargs)
    #axis("off")
    #show(block=False)

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


blastout=args.blast
outpath=args.outpath
infile=args.infile

contigs={}
contigs_align={}
contigs_spos={}
contigs_qpos={}
toremove=list()
cont_to_remove=list()
accessions={}
access_dict={}
bitscore_dict={}
best_hit={}

# blastn ref vs ref

seqs=to_dict_remove_dups(SeqIO.parse(infile, "fasta"))


blastf=open(blastout, "r")

print("Parsing blast catalog")
Bar = ProgressBar(len(blastf.readlines()), 60, '\t ')
barp=0
blastf=open(blastout, "r")
for l in blastf.readlines():
    barp=barp+1
    if l.startswith("#"):
        continue
    else:
        contigid=l.rstrip().split("\t")[0]
        bitscore=float(l.rstrip().split("\t")[11])
        sim=float(l.rstrip().split("\t")[2])
        aln=int(l.rstrip().split("\t")[3])
        sstart=int(l.rstrip().split("\t")[8])
        send=int(l.rstrip().split("\t")[9])
        qstart=int(l.rstrip().split("\t")[6])
        qend=int(l.rstrip().split("\t")[7])
        refid=l.rstrip().split("\t")[1]
        if sstart>send:
            sstart=int(l.rstrip().split("\t")[9])
            send=int(l.rstrip().split("\t")[8])
        if qstart>qend:
            qstart=int(l.rstrip().split("\t")[7])
            qend=int(l.rstrip().split("\t")[6])
    #print(contigid,refid,start,end)
    if contigid in list(contigs.keys()):
        if refid in list(contigs[contigid].keys()):
            tmp=contigs_spos[contigid][refid]
            cond_intersect=0
            for p in contigs_spos[contigid][refid]:
                min=int(p.split("-")[0])
                max=int(p.split("-")[1])
                #print(min,max)
                rhit=set(range(min,max))
                rc=set(range(sstart,send))
                intersection=rhit.intersection(rc)
                if len(intersection)==0:
                    pass
                else:
                    cond_intersect=cond_intersect+1
                    continue
            if cond_intersect==0:
                #old_score=float(contigs[contigid][refid])
                #new_score=float(old_score)+float(bitscore)
                old_aln=int(contigs_align[contigid][refid])
                new_aln=int(old_aln)+int(aln)
                #contigs[contigid][refid]=new_score
                contigs[contigid][refid].append(sim)
                contigs_align[contigid][refid]=new_aln
                contigs_spos[contigid][refid].append(str(sstart)+"-"+str(send))
                contigs_qpos[contigid][refid].append(str(qstart)+"-"+str(qend))
        else:
            contigs[contigid][refid]=[sim]
            contigs_align[contigid][refid]=aln
            contigs_spos[contigid][refid]=[]
            contigs_spos[contigid][refid].append(str(sstart)+"-"+str(send))
            contigs_qpos[contigid][refid]=[]
            contigs_qpos[contigid][refid].append(str(qstart)+"-"+str(qend))
    else:
        #print(contigid, refid)
        contigs.setdefault(contigid, {})
        contigs_align.setdefault(contigid, {})
        contigs_spos.setdefault(contigid, {})
        contigs_qpos.setdefault(contigid, {})
        #contigs[contigid][refid]=bitscore
        contigs[contigid][refid]=[sim]
        contigs_align[contigid][refid]=aln
        contigs_spos[contigid][refid]=[]
        contigs_spos[contigid][refid].append(str(sstart)+"-"+str(send))
        contigs_qpos[contigid][refid]=[]
        contigs_qpos[contigid][refid].append(str(qstart)+"-"+str(qend))
    Bar.update(barp)
print("done")
print("")

cond_l=args.overlap
cond_s=args.similarlity

bad=list()
bad_subject=dict()
bad_querries=dict()


connect_red=list()
connect_dict=dict()
connect_edges=list()

print("Finding edges nodes")
Bar = ProgressBar(len(contigs_align), 60, '\t ')
barp=0

for cont in contigs_align:
    barp=barp+1
    for sub in contigs_align[cont]:
        tmp=list()
        rm=list()
        rmlen=0
        if cont==sub:
            pass
        else:
            for p in contigs_qpos[cont][sub]:
                    if len(tmp)==0:
                        tmp.append(p)
                    else:
                        if p in tmp:
                            rm.append(p)
                        else:
                            tmp.append(p)
            if len(rm)>0:
                for s in rm:
                    a=int(p.split("-")[0])
                    b=int(p.split("-")[1])
                    rmlen=rmlen+(b+a)
        sublen=int(sub.split("_length_")[1].split("_cov_")[0])
        if (contigs_align[cont][sub]-rmlen)/sublen >= cond_l:
            if sum(contigs[cont][sub])/len(contigs[cont][sub]) >= cond_s*100:
                if cont not in connect_dict:
                    connect_dict[cont]=list()
                    connect_dict[cont].append(sub)
                else:
                    connect_dict[cont].append(sub)
                connect_edges.append((cont,sub))
                if cont!=sub:
                    connect_red.append((cont,sub))
                #print(cont,sub,round((contigs_align[cont][sub]-rmlen)/sublen,2),round(sum(contigs[cont][sub])/len(contigs[cont][sub]),2))
    Bar.update(barp)
print("done")
print("")

print("Computing graph and adjacency matrix")
G = nx.Graph()
#Gr=nx.Graph()
G.add_edges_from(connect_edges)
#Gr.add_edges_from(connect_red)
try:
    matrix = nx.to_scipy_sparse_array(G)
except:
    matrix = nx.to_scipy_sparse_matrix(G)
#matrix_numpy=matrix.toarray()
if isspmatrix(matrix):
    matrix2=matrix
else:
    matrix2=csr_matrix(matrix)

print("Clustering using MCL algorithm")
print("1) selection of the best inflation value giving the highest modularity score")
infval=[i / 10 for i in range(11, 36)]
bestQ=0
bestV=0
condinit=0
Bar = ProgressBar(len(infval), 60, '\t ')
barp=0
for inflation in infval:
    barp=barp+1
    result = mc.run_mcl(matrix2, inflation=inflation)
    clusters = mc.get_clusters(result)
    Q = mc.modularity(matrix=result, clusters=clusters)
    #print("inflation:", inflation, "modularity:", Q)
    Bar.update(barp)
    if condinit==0:
        bestQ=Q
        bestV=inflation
        condinit+=1
    else:
        if Q >= bestQ:
            bestQ=Q
            bestV=inflation
        else:
            continue
print("done")

print("2) finding clusters using the optimized cluster inflation value")
result = mc.run_mcl(matrix2, inflation=bestV)
clusters = mc.get_clusters(result)
print("done")
#plot of clusters

print("Save cluster plot")
draw_graph2(matrix, clusters, node_size=10, with_labels=False, edge_color="silver",width=0.01)
plt.savefig(os.path.join(outpath, str("clusters.pdf")),format='pdf')
plt.savefig(os.path.join(outpath, str("clusters.svg")),format='svg')
print("done")
print("")
# now select the best contigs from each clust
# condition
print("Selection of representative contigs from clusters")
Bar = ProgressBar(len(clusters), 60, '\t ')
barp=0
selected=list()
for clust in clusters:
    barp=barp+1
    unique=list()
    if len(clust)>1:
        for cont in clust:
            cfound = list(filter(lambda x: x.count(cont) > 0, clusters))
            cond_uniq=0
            for v in cfound:
                if v==clust:
                    continue
                else:
                    cond_uniq+=1
            if cond_uniq<1:
                unique.append(list(G.nodes)[cont])
    else:
        cont=clust[0]
        cfound = list(filter(lambda x: x.count(cont) > 0, clusters))
        cond_uniq=0
        for v in cfound:
            if v==clust:
                continue
            else:
                cond_uniq+=1
        if cond_uniq<1:
            unique.append(list(G.nodes)[cont])
    if len(unique)>0:
        selcont="NA"
        selscore=0
        for contig in unique:
            contig_lgth=int(contig.split("_length")[1].split("_cov")[0].replace("_",""))
            contig_cov=float(contig.split("_cov")[1].split("_")[1])
            val=contig_lgth*contig_cov
            if val >= selscore:
                selcont=contig
                selscore=val
        selected.append(selcont)
    Bar.update(barp)
print("done")

open(os.path.join(outpath, "selected_clean_metassemblies.fa"), 'w').close()

print("Exporting selected metacontigs")
Bar = ProgressBar(len(selected), 60, '\t ')
barp=0
for c in selected:
    barp=barp+1
    if c in list(seqs.keys()):
        header=">"+seqs[c].id
        seq=str(seqs[c].seq)
        if os.path.isfile(os.path.join(outpath,"selected_clean_metassemblies.fa")):
            with open(os.path.join(outpath,"selected_clean_metassemblies.fa"), 'a+') as file:
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
            with open(os.path.join(outpath,"selected_clean_metassemblies.fa"), 'w') as out:
                out.write(header+'\n')
                out.write(str(seq)+'\n')
    Bar.update(barp)
print("done")
