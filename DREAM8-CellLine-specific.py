# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 11:33:32 2013

@author: xinghualu
"""
import cPickle

from PyGibbCAMP import PyGibbCAMP
import os
import numpy as np

dataDir = 'ProcessedData/cellline.specific.tables/UACC812/'
nodeFile = "ProcessedData/name.matching.csv"


dataMatrix =  dataDir   + "data.matrix.csv"
perturbMatrix =  dataDir   + "perturbation.table.csv"
if os.path.exists(dataDir + "missDataMatrix.csv"): 
    missDataMatrix =  dataDir + "missDataMatrix.csv"
else:
    missDataMatrix = None

nParents = 3
for alpha in np.arange(1., .4, -0.1, dtype=np.float):
    for i in range(5):
        pickleFileName = dataMatrix + ".chain" + str(i) + ".nParents" + str(nParents) +  ".alpha-" + str(alpha) + ".pickle"
        net = PyGibbCAMP(nodeFile = nodeFile, dataMatrixFile = dataMatrix, perturbMatrix = perturbMatrix, missingDataMatrix= missDataMatrix)
        net.trainGibbsEM(pickleDumpFile = pickleFileName, nParents = nParents, nChains = 1, alpha=alpha, maxIter = 300)

