from scipy import stats
import numpy as np


file="cov_catalog.txt"

cov_dict={}
covs=list()
ids=list()

with open(file) as f:
    reftab = f.readlines()
for line in reftab:
    l=line.rstrip().split("\t")
    if "NODE" in l[0]:
        cont=l[0]
        cov=float(l[5])
        depth=float(l[6])
        cov_dict[cont]=depth
        covs.append(depth)
        ids.append(cont)
    else:
        continue

z = np.abs(stats.zscore(covs))
threshold = 3
# Position of the outlier
print(np.where(z > 3))

# voir pour prendre que ceux qui sont en excès sup. pas les plus faible.. ?

# mkdir Calling/covs -->
# mettre la fonction en ouvrant le catalogue par indivs
# faire un dict global avec NODE list 0,1,NA si excès (1) ou non (0), manquants
# ensuite si (1)/(0+1) > 50% --> identifie en bad_depth

# voir le cov mean sur la ref pour l'ensemble --> si <50% --> identifie en bad_cov / mettre des 0.0 en fonction de la longueur
# file number = indivs

# parser le merge catalog VCF
# sites ht / contigs
