import multiprocessing
import scipy
import math
import sys
import numpy as np

data_file = sys.argv[1]
size = int(sys.argv[2])
gen = int(sys.argv[3])

d = np.loadtxt(data_file)
#N = len(d[:,1])
N = gen*size

x = []
y = []
for i in range(0, N):
    iDom = 0
    for j in range(0, N):
        if (i != j):
            if (d[j,1] < d[i,1]) and (d[j,2] < d[i,2]):
                iDom = iDom + 1
    if (iDom == 0):
        print(i, d[i,1], d[i,2], d[i,3], d[i,4], d[i,5])
