# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 11:33:32 2013

@author: xinghualu
"""


from PyGibbCAMP import PyGibbCAMP

import numpy as np

dataDir = 'ProcessedData/insilico/'

nodeFile = "ProcessedData/insilico/name.matching.csv"
dataMatrix =  dataDir  + "data.matrix.insilico.csv"
perturbMatrix =  dataDir  + "perturbation.table.insilico.csv"
missDataMatrix = None
net = PyGibbCAMP(nodeFile = nodeFile, dataMatrixFile = dataMatrix, perturbMatrix = perturbMatrix, missingDataMatrix= missDataMatrix)

nParents = 3  
for alpha in np.arange(0.1, 1, 0.1, dtype=np.float):
    for i in range(5):
        pickleFileName = dataMatrix + ".chain" + str(i) + ".nParents" + str(nParents) +  ".alpha-" + str(alpha) + ".pickle"
        net.trainGibbsEM(pickleDumpFile = pickleFileName, nChains = 1, alpha=alpha, maxIter = 500, nParents = nParents)

