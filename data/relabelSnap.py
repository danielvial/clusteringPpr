# ensures edge list from SNAP has consecutive node labels 0,1,...,n-1

import numpy as np
import sys

# import edge list
edges = np.loadtxt(sys.argv[1]+'.txt',int)

# create node map
nodeMap = np.unique(edges)

# relabel nodes
n = nodeMap.shape[0]
key = np.arange(n)
index = np.digitize(edges.ravel(),nodeMap,right=True)
newEdges = key[index].reshape(edges.shape)

# save cleaned file and number nodes
np.savetxt(sys.argv[1]+'-r.txt',newEdges,delimiter='\t',fmt='%d')
