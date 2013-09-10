# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 11:33:32 2013

@author: xinghualu
"""
import cPickle
import numpy as np
from PyGibbCAMP import PyGibbCAMP

dataDir = 'ProcessedData/'

nodeFile = dataDir + "name.matching.csv"
dataMatrix = dataDir + "data.matrix.normalized.csv"
missDataMatrix = dataDir + "missDataMatrix.csv"
perturbMatrix = dataDir + "perturbation.table.csv"

nParents = 2
for alpha in np.arange(1., .5, -0.1, dtype=np.float):
    for i in range(5):
        pickleFileName = dataMatrix + ".chain" + str(i) + ".nParents" + str(nParents) +  ".alpha-" + str(alpha) + ".pickle"
        model = PyGibbCAMP(nodeFile = nodeFile, dataMatrixFile = dataMatrix, perturbMatrix = perturbMatrix, missingDataMatrix= missDataMatrix)
        model.trainGibbsEM(pickleDumpFile = pickleFileName, nParents = nParents, nChains = 1, alpha= alpha, maxIter = 500)

#cPickle.dump(model, open("final-model-09-06-13-alpha.05.nParent.4.pickle", 'wb'))
