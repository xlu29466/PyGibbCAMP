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
    #  @param dataMatrixFile  A string to data
    def __init__(self, edgeList, dataMatrixFile = None, perturbMatrix = None):
        self.network = None
        self.data = None
        
        if nodeList:
            self.initNetwork(nodeList)
            
        if self.network and dataMatrixFile:
            self.assocData(dataMatrixFile)            

    
    ## Create a network
    #  @param  edgeList  A list of PySigNetNode objects
    #
    #  This function take a list of tuples consisting of a source and a sink to
    #  create a network consisting of the nodes.  Both source and sink should
    #  should be instances of a PySigNetNode object.  This function
    #  create a "networkx" instance, with node indexed by the name of the nodes
    #  and each network node is indexed by the name of the node
    def initNetwork(self, edgeList):
        if not nodeList:
            raise Exception("Calling 'intiNetwork' without empty object list")
        if self.network:
            self.network.clear()
            
        self.network = nx.DiGraph()
        
        for source, sink in nodeList:
            self.network.add_node(source.name, 'data' = source)
            self.network.add_node(sink.name, 'data' = sink)
            self.network.add_edge(source.name, sink.name)
            
            
    ## Return a node object by name
    # @param  name  The name of the node
    def getNodeByName(self, name):
        if name in self.network:
            return self.network[name]['data']
        else:
            return None
            
            
    def assocData(self, dataFileName):
        if not dataMatrixFile:
            raise Exception ("Calling 'assocData' with an empty file name")
            
        try:
            f = open(dataFileName)
            lines = f.readlines()
        except IOError:
            print "Fail to read  file " + dataFileName
        
        colnames = None
        if len(lines) == 1:  # Mac version csv, with "\r" as return
            lines = lines[0].split("\r")
            colnames = lines.pop(0).split(',') # split  header and extract colnames
            map(lambda x: x.rstrip(), lines)  # remove the "\r"
            lines = "\n".join(lines)  # use "\n" to join lines
        else:
            colnames = lines.pop(0).split(',')
            lines = "".join(lines)
            
        # read in data and generate a numpy data matrix
        self.data = np.genfromtxt(StringIO(lines), delimiter = ",", usecol=tuple(range(1, len(colnames))))
            
        #check in which column the data for a node in graph locates
        self.dictNode2MatrixIndx = dict()  
        for node in self.network:
            nodeIndex = colnames.index(node)
            if not nodeIndex:  # there is no data for the node
                raise Exception("The data for node " + node + " is missing.  Quit!")
            self.dictNode2MatrixIndx[node] = nodeIndex
        
        # find column indices for the predecessors of a node
        self.dictParentOfNodeToMatrixIndx = dict()  # a diction of list
        for node in self.network:
            self.dictParentOfNodeToMatrixIndx[node] = list() 
            preds = self.network.predecessors(node)

            for p in preds:                
                self.dictParentOfNodeToMatrixIndx[node].append(colnames.index[p])             
                

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
    
    
    
    
    
