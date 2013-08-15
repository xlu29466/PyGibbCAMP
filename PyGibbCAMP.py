# -*- coding: utf-8 -*-
"""
@ PyCAMP  Python Causal Modeling of Pathways, a python implmentation for modeling
causal relationship bewtween cellular signaling proteins, particularly phosphorylated
proteins based on reverse phase protein array (RPPA) data.


Created on Wed Aug 14 19:16:25 2013

@author: Xinghua  Lu
"""

import networkx as nx
import numpy as np
from numpy import matlib
from rpy2 import robjects 
import rpy2.robjects.numpy2ri
import os


## Class PyCAMP
#  This class represent the regulatory network among signaling proteins, in 
#  particular protein kinases and their target proteins.  Given observed
#  RPPA data, which reflect the concentration of total or phosphorylated
#  proteins of interest in a cell population, the CAMP model represent the 
#  observed total fluoresecent signals of POIs, the phosphorylation states of
#  POIs, and the activation states of POIs.  The goal is to learn if a causal
#  relationship between the activation state of a POI and the phosphorylation
#  state of POIs.


class PyGibbCAMP:
    
    ## Constructor
    #  @param nodeFile  A string of pathname of file containing nodes.  The 
    #                   name, type, measured
    #  @param edgeFile  A string of pathname of file containing edges
    #  @param 
    def __init__(self, nodeFile = None, edgeFile = None, dataMatrixFile = None, perturbMatrix = None):
        """
        """
        self.network = None
        self.data = None
        
        if nodeFile and edgeFile:
            self.initNetwork(nodeFile, edgeFile)
            
        if self.network and dataMatrixFile:
            self.assocData(dataMatrixFile)            

    
    
    def initNetwork(self, nodeFile, edgeFile):
        if not nodeFile or not edgeFile:
            raise Exception("Calling 'intiNetwork' without empty file names")
        if self.network:
            self.network.clear()
            
        self.network = nx.DiGraph()
        
        # add nodes and edges
            
    def assocData(self, dataMatrixFile):
        if not dataMatrixFile:
            raise Exception ("Calling 'assocData' with an empty file name")
    
    ## Add an edge between two nodes.  Can be called to perform graph search
    #  @param source  The source node of th edge
    #  @param sink    The sink node of the edge
    def addEdge(self, source, sink):
        self.network.add_edge(source, sink)

    ## Remove an edge between two nodes.  Can be called to perform graph search
    #  @param source  The source node of th edge
    #  @param sink    The sink node of the edge        
    def deleteEdge(self, source, sink):
        self.network.remove_edge(source, sink)
        
    ## Calculate the marginal probability of observing the measured data by
    #  integrating out all possible setting of latent variable states and 
    #  model parameters.
    def calcObservedMarginal(self):
        # this can be easily achieved by taking expectation of observed 
        # phosphorylation states 
        pass
    
    ## Perform Gibbs sampling to perform EM inference of network model
    def trainModelByGibbs(self):
        pass
    
    
    ## Gibbs sampling to update latent variables
    def _updateLatentByGibbs(self):
        pass
    
    ## Inference of model parameters through logistic regression 
    def _updateParams(self):
        pass
    
    
    def toXGMML(self):
    
    
    
    
    
