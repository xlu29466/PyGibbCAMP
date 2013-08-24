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
import os, math, pickle
from SigNetNode import SigNetNode

import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()   # enable directly pass numpy arrary or matrix as arguments to rpy object
R = robjects.r                      # load R instance
R.library("glmnet")
glmnet = R('glmnet')                # make glmnet from R a callable python object

R.library("mixtools")
normalmixEM = R('normalmixEM')


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
    #  @param edgeFile  A list of tuples, each containing a source and sink node 
    #                   of an edge
    #  @param dataMatrixFile  A string to data
    def __init__(self, nodeFile = None, edgeFile = None, dataMatrixFile = None, perturbMatrix = None):
        self.network = None
        self.data = None
        self.nChains = 1
        
        if nodeList:
            self.initNetwork(nodeFile, edgeFile)
            
        if self.network and dataMatrixFile and perturbMatrix:
            self.assocData(dataMatrixFile)            

    
    ## Create a network
    #  @param  edgeList  A list of PySigNetNode objects
    #
    #  This function take a list of tuples consisting of a source and a sink to
    #  create a network consisting of the nodes.  Both source and sink should
    #  should be instances of a PySigNetNode object.  This function
    #  create a "networkx" instance, with node indexed by the name of the nodes
    #  and each network node is indexed by the name of the node
    def initNetwork(self, nodeFile, edgeFile):
        if not nodeFile:
            raise Exception("Calling 'intiNetwork' with empty nodeFile name")
        if not edgeFile:
            raise Exception("Calling 'initNetwork' with empty edge file" )
            
        if self.network:
            self.network.clear()
            
        self.network = nx.DiGraph()
        
        try:
            nf = open(nodeFileName, "r")
            nodeLines = nf.readlines()
            nf.close()
        except IOError:
            print "Failed to open the file containing nodes"
            return
            
        if len(nodeLines) == 1:  # Mac files end a line with \r instead of \n
            nodeLines = nodeLines[0].split("\r")
           
        try:
            ef = open(edgeFileName, 'r')
            edgeLines = ef.readlines()
            ef.close()
        except IOError:
            print "Fail to open file " + nodeFileName
            return
        if len(edgeLines) == 1:
            edgeLines = edgeLines[0].split("\r")  

        # creat a networkx graph 
        self.network = nx.DiGraph()
        print "Creating network"

        # parse nodes
        for line in nodeLines:
            #print line
            protein, antibody = line.split(',')
            fluo = antibody + 'F'
            if protein not in self.network:
                self.network.add_node(protein, nodeObj = SigNetNode(protein, 'activityState', False))
            self.network.add_node(antibody, nodeObj= SigNetNode(antibody, 'phosState', False))
            self.network.add_node(fluo, nodeObj = SigNetNode(fluo, 'fluorescence', True))
            self.network.add_edge(antibody, protein)
            self.network.add_edge(antibody, fluo)
        print "Added " + str(len(nodeLines)) + "antibody nodes"

        # parse edges
        for line in edgeLines:
            #print "Current edge line: " + line
            source, sink = line.split(',')
            if source not in self.network or sink not in self.network:
                print "Cannot add edges between nodes not in the network"
                continue
            self.network.add_edge(source, sink)
        print "Added " + str(len(edgeLines)) + " edges.  Done with creating network"        
            
            
    def assocData(self, dataFileName, perturbMatrix):
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
    def removeEdge(self, source, sink):
        self.network.remove_edge(source, sink)
        
        
    ## Calculate the marginal probability of observing the measured data by
    #  integrating out all possible setting of latent variable states and 
    #  model parameters.
    def calcObservedMarginal(self):
        # this can be easily achieved by taking expectation of observed 
        # phosphorylation states 
        pass
    
    
    def initParams(self):
        """Initialize the parameter vector associated with each node 
           with a random vector sampled from standard Gaussian.
           
           The parameters associated with each node for each MCMC chain is 
           stored in two dimensional numpy array with nChain * nPredecessors 
        """
        print "Initialize parameters associated with each node in each MCMC chain"
        for nodeId in self.network: 
            if self.network.node[nodeId]['nodeObj'].type == 'fluorescence':                
                # To do, init the mean and sd of fluorescence signal
                mixGuassians = normalmixEM(self.data[:,self.dictNode2MatrixIndx[nodeId]])
                nodeObj.mus = matlib.repmat(mixGuassian[2], self.nChain, 1)
                nodeObj.sigmas = matlib.repmat(mixGaussian[3], self.nChains, 1)               
            else:
                parents = self.network.predecessors(nodeId)
                if len(parents) > 0:
                    self.dictNodeParams[nodeId] = np.random.randn(self.nChains, len(parents) + 1)
                else:
                    self.dictNodeParams[nodeId]  = None
                
        

    ## Perform Gibbs sampling to perform EM inference of network model
    def trainModelByGibbs(self):
        pass




    def _calcNodeMarginal(self, nodeId, c):
        """
        Calculate the marginal probability of a node's state set to "1" 
        
        args:
             nodeId   A string id of the node of interest
             c        An integer indicate the chain from which the parameter 
                         vector to be used  
        """
        nCases, nVariables = np.shape(self.data)
        
        # collect the state of the predecessors of the node
        predIndices = self.dictParentOfNodeToMatrixIndx[nodeId]  
        logProbOneCondOnParents = 0
        logProbZeroCondOnParents = 0
        if len(predIndices) > 0:  # if the node has parents  
            # calculate p(node = 1 | parents);   
            nodeParams = self.dictNodeParams[nodeId][c,:] 
            predStates =  np.column_stack((np.ones(nCases), self.nodeStates[c][:, predIndices])) 
            pOneCondOnParents = 1 / (1 + np.exp( - np.dot(predStates, nodeParams)))  
            logProbOneCondOnParents  = np.log(pOneCondOnParents)
            logProbZeroCondOnParents = np.log(1 - pOneCondOnParents)

        # collect  evidence from all children 
        logProbChildCondOne = 0  # the prob of child conditioning on current node == 1
        logProdOfChildCondZeros = 0
        children = self.network.successors(nodeId)
        curNodeIndx = self.dictNode2MatrixIndx[nodeId]

        if len(children) > 0:
            for child in children:   
                nodeObj = self.network.node[child]['nodeObj']
                childNodeIndx = self.dictNode2MatrixIndx[child]
                if nodeObj.type == "fluorescence":
                    logProbChildCondOne += self.nodeStates[:,curNodeIndx] * 0.5 * nodeObj.sigma[1] * np.square(self.data[:,childNodeIndx] - nodeObj.mus[1]) 
                    logProdOfChildCondZeros +=  (1 - self.nodeStates[:,curNodeIndx] ) * 0.5 * nodeObj.sigma[0] * np.square(self.data[:,childNodeIndx] - nodeObj.mus[0])
                else:  
                    # collect data and parameters associated with the node
                    childNodeParams = self.dictNodeParams[child][c,:] 
                    curChildStates = self.nodeStates[c][:, self.dictNode2MatrixIndx[child]]                    
                        
                        # find out predecessors and the position of curNode in the list
                        childPredIndices = self.dictParentOfNodeToMatrixIndx[child] 
                        curNodePosInPredList = childPredIndices.index(curNodeIndx)  
        
                    # Collect states of the predecessors of the child
                    # Set the state of current node to ones 
                    childPredStatesWithCurSetOne = self.nodeStates[c][:, childPredIndices]    
                    childPredStatesWithCurSetOne[:, curNodePosInPredList] = np.ones(nCases)  
                    childPredStatesWithCurSetOne = np.column_stack((np.ones(nCases), childPredStatesWithCurSetOne)) # padding data with a column ones as bias
                    pChildCondCurNodeOnes = 1 / (1 + np.exp(-np.dot(childPredStatesWithCurSetOne, childNodeParams)))
                    logProbChildCondOne += np.log (curChildStates * pChildCondCurNodeOnes + (1 - curChildStates) * (1 - pChildCondCurNodeOnes))
        
                    # set the state of the current node (nodeId) to zeros 
                    childPredStatesWithCurSetZero = self.nodeStates[c][:,childPredIndices]
                    childPredStatesWithCurSetZero [:, curNodePosInPredList] = np.zeros(nCases)
                    childPredStatesWithCurSetZero = np.column_stack((np.ones(nCases), childPredStatesWithCurSetZero))
                    pChildCondCurNodeZeros = 1 / (1 + np.exp(- np.dot(childPredStatesWithCurSetZero, childNodeParams))) 
                    logProdOfChildCondZeros += np.log(curChildStates * pChildCondCurNodeZeros + (1 - curChildStates) * (1 - pChildCondCurNodeZeros))

        # now we can calculate the marginal probability of current node 
        curNodeMarginal = 1 / (1 + np.exp(logProbZeroCondOnParents + logProdOfChildCondZeros - logProbOneCondOnParents - logProbChildCondOne))
        return curNodeMarginal
    

    def parseGlmnetCoef(self, glmnet_res):        
        """ Parse the 'beta' matrix returned by calling glmnet through RPy2.
            Return the first column of 'beta' matrix of the glmnet object 
            with all none zero values returned by the glmnet
            """
                    
        # Read in lines of beta matrix txt, which is a nVariables * nLambda.
        # Since we call glmnet by padding x with a column of 1s, we only work
        # with the 'beta' matrix returned by fit
        betaLines = StringIO(str(glmnet_res.rx('beta'))).readlines()
        dimStr = re.search("\d+\s+x\s+\d+", betaLines[1]).group(0)
        if not dimStr:
            raise Exception("'parse_glmnet_res' could not determine the dims of beta")
        nVariables , nLambda = map(int, dimStr.split(' x ')) 
        betaMatrix = np.zeros( (nVariables, nLambda), dtype=np.float)
        
        # glmnet print beta matrix in mulitple blocks with 
        # nVariable * blockSize
        blockSize = len(betaLines[4].split()) - 1
        curBlockColStart = - blockSize
        for line in betaLines:  #read in blocks
            m = re.search('^V\d+', line)
            if not m:  # only find the lines begins with 'V\d'
                continue
            else:
                rowIndx = int(m.group(0)[1:len(m.group(0))]) 
            if rowIndx == 1:
                curBlockColStart += blockSize
                
            # set 'rowIndx' as start from 0
            rowIndx -= 1

            fields = line.rstrip().split()
            fields.pop(0)
            if len(fields) != blockSize:
                blockSize = len(fields)
            for j in range(blockSize):
                if fields[j] == '.':
                    continue
                else:
                    betaMatrix[rowIndx, curBlockColStart + j] = float(fields[j])                 
                            
        # scan through the beta matrix and return the first all none zero 
        # or the last column                    
        for j in range(nLambda):
            if not np.any(betaMatrix[:, j] == 0.):
                break
        
        return betaMatrix[:,j]        
      
        
    def _updteParams(self, alpha = 0.05, keepRes = False):
        # Update the parameter associated with each node, p(n | Pa(n)) using logistic regression,
        # using expected states of precessors as X and current node states acrss samples as y
        nCases, nVariables = np.shape(self.data)
        for nodeId in self.network:
            predIndices = self.dictParentOfNodeToMatrixIndx[nodeId]
            if len(predIndices) == 0:
                continue
            
            nodeIdx = self.dictNode2MatrixIndx[nodeId]
            nodeObj = self.network.node[nodeId]['nodeObj']
            
            if keepRes:
                self.network.node[nodeId]['nodeObj'].fitResults = []

            for c in range(self.nChains): 
                expectedPredState = self.expectedStates[c][:, predIndices]
                curNodeData = self.data[:,nodeIdx]
                if self.network.node[nodeId]['nodeObj'].type == "fluorescence":
                    nodeObj.mus[c,0] = np.mean ((1-expectedPredState) * curNodeData)
                    nodeObj.mus[c, 0] = np.std ((1-expectedPredState) * curNodeData)
                    nodeObj.sigmas[c, 1] = np.mean (expectedPredState * curNodeData)
                    nodeObj.sigmas[c, 1] = np.mean (expectedPredState * curNodeData)
                
                else:   
                    x = np.column_stack((np.ones(nCases), expectedPredState))
                    y = self.nodeStates[c][:, nodeIdx]
                    
                    #check if all x and y are of same value, which will lead to problem for glmnet
                    rIndx = map(lambda z: int(math.floor(z)), np.random.rand(3) * nCases)
                    if sum(y) == nCases:                        
                        y[rIndx] = 0                        
                    elif sum( map(lambda x: 1 - x, y)) == nCases:
                        y[rIndx] = 1        
                    y = robjects.vectors.IntVector(y)
                
                    allOnes = np.where(np.sum(x[:, 1:nVariables],0) == nCases)
                    for c in allOnes[0]:
                        rIndx = map(lambda z: int(math.floor(z)), np.random.rand(3) * nCases)
                        x[rIndx, c+1] = 0 
                    allZeros = np.where(np.sum(np.ones(np.shape(x)) - x, 0) == nCases) 
                    for c in allZeros[0]:
                        rIndx = map(lambda z: int(math.floor(z)), np.random.rand(3) * nCases)
                        x[rIndx, c] = 1
                    
                
                    # call logistic regression using glmnet from Rpy
                    fit = glmnet (x, y, alpha = .05, family = "binomial", intercept = 0)
                    # extract coefficients from Rpy2 vector object
                    self.dictNodeParams[nodeId][c,:] = self.parseGlmnetCoef(fit) 
                    
                    if keepRes:
                        self.network.node[nodeId]['nodeObj'].fitResults.append(fit)
        
                
                
    
    ## Gibbs sampling to update latent variables
    def _updateLatentByGibbs(self):
        pass
    
    ## Inference of model parameters through logistic regression 
    def _updateParams(self):
        pass
    
    
    def toXGMML(self):
    
    
    
    
    
