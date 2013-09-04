# -*- coding: utf-8 -*-
"""
@ PyCAMP  Python Causal Modeling of Pathways, a python implmentation for modeling
causal relationship bewtween cellular signaling proteins, particularly phosphorylated
proteins based on reverse phase protein array (RPPA) data.  This model is designed
to model the signal transduction through series of protein phosphorylation cascade,
in which phosphorylation of a protein often activate the protein, which in turn
lead to phosphorylation of other proteins.  This model represent
the phosphorylation state(s) and activation state of a protein separately such that the model
is capable of capture the fact that, at certain time, phosphorylation of a protein
can be decoupled by drug and inhibitors. 


Created on Wed Aug 14 19:16:25 2013

@author: Xinghua  Lu
"""

import networkx as nx
import numpy as np
from numpy import matlib
from rpy2 import robjects 
import  math, cPickle, re
from SigNetNode import SigNetNode
from StringIO import StringIO
from NamedMatrix import NamedMatrix

import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()   # enable directly pass numpy arrary or matrix as arguments to rpy object
R = robjects.r                      # load R instance
R.library("glmnet")
glmnet = R('glmnet')                # make glmnet from R a callable python object

R.library("mixtools")
normalmixEM = R('normalmixEM')




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
        
        self.dictPerturbEffect = { 'TR.GSK690693.Inhibitors' : ('AKT1',	0), \
        'TR.GSK690693_GSK1120212.Inhibitors' : [('AKT1', 0), ('MAP2K1', 0)],\
        'TR.EGF.Stimuli' : ('EGFR', 1), 'TR.FGF1.Stimuli': ('EGFR', 1)}

        
        if nodeFile and edgeFile:
            self.initNetwork(nodeFile, edgeFile)
            
        if self.network and dataMatrixFile and perturbMatrix:
            self.assocData(dataMatrixFile, perturbMatrix)            

    
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
            nf = open(nodeFile, "r")
            nodeLines = nf.readlines()
            nf.close()
        except IOError:
            raise Exception( "Failed to open the file containing nodes")
            return
            
        if len(nodeLines) == 1:  # Mac files end a line with \r instead of \n
            nodeLines = nodeLines[0].split("\r")
           
        try:
            ef = open(edgeFile, 'r')
            edgeLines = ef.readlines()
            ef.close()
        except IOError:
            raise Exception ("Fail to open file " + edgeFile)
            return
        if len(edgeLines) == 1:
            edgeLines = edgeLines[0].split("\r")  

        # creat a networkx graph 
        self.network = nx.DiGraph()
        print "Creating network"

        self.dictProteinToAntibody = dict()
        self.dictAntibodyToProtein = dict()
        # parse nodes
        for line in nodeLines:
            #print line
            protein, antibody = line.rstrip().split(',')
            if protein not in self.dictProteinToAntibody:
                self.dictProteinToAntibody[protein] = []
            self.dictProteinToAntibody[protein].append(antibody)
            self.dictAntibodyToProtein[antibody] = protein
            
            fluo = antibody + 'F'
            if protein not in self.network:
                self.network.add_node(protein, nodeObj = SigNetNode(protein, 'activeState', False))
            self.network.add_node(antibody, nodeObj= SigNetNode(antibody, 'phosState', False))
            self.network.add_node(fluo, nodeObj = SigNetNode(fluo, 'fluorescence', True))
            self.network.add_edge(antibody, protein)
            self.network.add_edge(antibody, fluo)
        print "Added " + str(len(nodeLines)) + "antibody nodes"

        # parse edges
        for line in edgeLines:
            #sorce should be protien name, sink can be protein name or an antibody
            source, sink = line.rstrip().split(',')
            if source not in self.dictProteinToAntibody or sink not in self.network:
                print "Cannot add edges between nodes not in the network"
                continue
            
            if self.network.node[sink]['nodeObj'].type == "phosState":
                self.network.add_edge(source, sink)
            else: # both source and sink are proteins, 
                for ab in self.dictProteinToAntibody[sink]: # look up antibodies of sink
                    self.network.add_edge(source, ab)
                    
        cycles = nx.simple_cycles(self.network)
        for cyc in cycles:
            if len(cyc) == 4:  #a cycle involve 2, cut one protein to antibody edge  
                for n in cyc:
                    if self.network.node[n]['nodeObj'].type == 'activeState':
                        for child in self.network.successors(n):
                            if child in cyc:
                                self.network.remove_edge(n, child)
                                break
                
        print "Added " + str(len(self.network.edges())) + " edges. "  
        
            
    ## Parsing data files
    #  @param dataFileName  File name for RPPA data.  The format should be  
    #                       a nCases * nAntibody CVS matrix file
    #  @param perturbMatrix  File containing description of perturbation 
    #        conditions.  The file format: each row corresponds
    #        to a perturbation condition containing three fields:
    #        1) canse id, which should be in the case in the data matrix 
    #        2) perturbed protein, which should be in the node file
    #        3) value of perturbed protein, '1' --> activation, '0' --> inactivation
    def assocData(self, dataFileName, perturbMatrix, missingDataMatrix):
        # parse data mastrix by calling NamedMatrix class
        self.data = NamedMatrix(dataFileName)
        self.dataFileName = dataFileName
        
        self.perturbData = NamedMatrix(perturbMatrix)
        self.missingDataMatrix = NamedMatrix(missingDataMatrix)
           

        
    ## Init parameters of the model
    #  In Bayesian network setting, the joint probability is calculated
    #  through the product of a series conditional probability.  The parameters
    #  of the PyCAMP model defines p(x | Pa(X)).  For observed fluorescent node
    #  the conditional probability is a mixture of two Gaussian distribution.  
    #  therefore, the parameters are two pairs of mu and sigma.  For
    #  the hidden variables representing phosphorylation states and activation
    #  states of proteins, the conditional probability is defined by a logistic
    #  regression. Therefore, the parameters associated with such a node is a 
    #  vector of real numbers.
    # 
    def _initParams(self):
        print "Initialize parameters associated with each node in each MCMC chain"
        for nodeId in self.network: 
            self._initNodeParams(nodeId)
            
    def _initNodeParams(self, nodeId):
        nodeObj = self.network.node[nodeId]['nodeObj']
        if nodeObj.type == 'fluorescence':                
            # Estimate mean and sd of fluo signal using mixture model
            if nodeId in self.missingDataMatrix.getColnames():
                nodeData = self.data.getValuesByCol(self.missingDataMatrix[:,nodeId] == 0, nodeId)
            else:
                nodeData = self.data.getValuesByCol(nodeId)
            mixGuassians = normalmixEM(nodeData)
            # mus and sigmas are represented as nChain x 2 matrices
            nodeObj.mus = matlib.repmat(np.array(mixGuassians[2]), self.nChain, 1)
            nodeObj.sigmas = matlib.repmat(np.array(mixGuassians[3]), self.nChains, 1)               
        else:
            preds = self.network.predecessors(nodeId)
            if len(preds) > 0:
                nodeObj.paramNames = preds
                nodeObj.params = np.random.randn(self.nChains, len(preds) + 1)
            else:
                nodeObj.params  = None
                
    
    ## Initialize latent variables
    #    
    #
    def _initHiddenStates(self):
        hiddenNodes = [n for n in self.network if not self.network.node[n]['nodeObj'].bMeasured]
        nCases, nAntibody = self.data.shape()
        caseNames = self.data.getRownames()
        
        self.hiddenNodeStates = list()
        for c in range(self.nChains):
            tmp = np.zeros((nCases, len(hiddenNodes)))
            tmp[np.random.rand(nCases, len(hiddenNodes)) < 0.25] = 1
            self.hiddenNodeStates.append(NamedMatrix(npMatrix = tmp, colnames = hiddenNodes, rownames = caseNames))
        
        
    ## Calculate the marginal probability of observing the measured data by
    #  integrating out all possible setting of latent variable states and 
    #  model parameters.
    def calcEvidenceLikelihood(self):
        # this can be easily achieved by taking expectation of observed 
        # phosphorylation states 
        loglikelihood = 0        
        obsNodes = [n for n in self.network if self.network.node[n]['nodeObj'].bMeasured]
        for c in range(self.nChains):
            for nodeId in obsNodes:
                curNodeData = self.data.getValuesByCol(nodeId)
                pred = self.network.predecessors(nodeId)
                predStates = self.hiddenNodeStates[c].getValuesByCol(pred)
                nodeObj = self.network.node[nodeId]['nodeObj'] # where parameters are saved

                loglikelihood += sum(- predStates * (np.log(nodeObj.sigma[c, 1])\
                + 0.5 * np.square(curNodeData - nodeObj.mus[c, 1]) / np.square(nodeObj.sigma[c, 1])) \
                + (1 - predStates) *  (- np.log(nodeObj.sigma[c, 0])\
                + 0.5 * np.square(curNodeData - nodeObj.mus[c, 0]) / np.square(nodeObj.sigma[c, 0])))
        
        loglikelihood /= self.nChains
            
        
        
    ## Perform graph search
    def train(self, nChains = 10, alpha = 0.2, nParents = 4, nSamples = 10, pickleDumpFile = None, maxIter = 1000):
        self.nChains = nChains
        self.alpha = alpha  
        self.likelihood = list()
        self.nSamples = nSamples
        self.nParents = nParents
        if pickleDumpFile:
            self.pickleDumpFile = pickleDumpFile
        else:
            self.pickleDumpFile = self.dataFileName + "alpha" + str(self.alpha) +  ".pickle"        
        
        
        # Starting EM set up Markov chains  to train a model purely based on prior knowledge
        self._initHiddenStates()
        self._initParams()
        self._gibbsEM(maxIter)
        self.trimEdgeByConsensus(.9)
        
        # search for possible structures by modifying cannonical network
        nCases, nAntibodies = np.shape(self.data)
        coefMatrix = np.corrcoef(np.transpose(self.data))
        for i in range(nAntibodies):
            antibodyI = self.data.getColnames()[i]
            proteinI = self.dictAntibodyToProtein[antibodyI]

            corrI = coefMatrix[:,i].tolist() # corr of I to other
            corrIsorted = reversed(sorted(corrI))
            corrIsorted.pop(0)  #correlation with self 
            # foreach node, we consider 3 top correlated as potential parents
            
            for j in range(5): 
                cycleEdge = None
                currLikelihood = np.mean(self.liklihood[-5:-1])
                indxJ = corrI.index(corrIsorted[j])
                antibodyJ = self.data.getColnames()[indxJ]
                proteinJ = self.dictAntibodyToProtein[antibodyJ]
                if proteinI == proteinJ: # two antibodies again a common protein
                    continue                    
                if self.network.has_edge(proteinJ, antibodyI):
                    continue
                if self.network.has_edge(proteinI, antibodyJ):
                    cycleEdge = (proteinI, antibodyJ)

                self._initNodeParams(antibodyI)
                self._gibbsEM(100)
                newLikelihood = np.mean(self.liklihood[-5:-1])
                if (newLikelihood - currLikelihood)  < 0:
                    self.network.remove_edge(proteinJ, antibodyI)
                    if cycleEdge:
                        self.network.add_edge(cycleEdge[0], cycleEdge[1])
        
        # now try to delete edges that does contribute to evidence
        self.trimEdgeByConsensus(.9)
        return self  
              
        
    ## Perform Gibbs sampling to perform EM inference of network model
    #  
    #   @param nChains  The number of Markov chains in the model
    #   @alpha  The 'alpha' parameter, which will be used to train glmnet model for logistic regression
    #   @nSamples  The number of samples to collect before EM update parameters
    #   @pickleDumpFile The file name to save the best model while training
    #   @maxIter  The number of maximal interation to train the model
    #
    #  This is the main function of the class, which handles all training aspects
    #  
    def _gibbsEM(self, maxIter):
        self.likelihood = list()        
        
        self.expectedStates = list()
        nCases, nAntibodies = np.shape(self.data)
        for c in range(self.nChains):                  
            # each chain collect expected statistics of nodes from samples along the chain
            self.expectedStates.append(np.zeros(np.shape(self.data)))

        print "Starting EM: alpha = " + str(self.alpha) + "; nChains = " + str(self.nChains) + "; nSamples = " + str (self.nSamples)
        optLikelihood = float("-inf")
        bConverged = False
        sampleCount = 0
        self._updteParams(self.alpha)
        for nIter in range(maxIter):
            likelihood = self.calcEvidenceLikelihood()
            print "nIter: " + str(nIter ) + "; log likelihood of evidence: " + str(likelihood)
            self.likelihood.append(likelihood)
            
            # E-step of EM
            self._updateStates()            
            if  (nIter+1) % 2 == 0: # we collect sample every other iteration
                sampleCount += 1
                for c in range(self.nChains):
                    self.expectedStates[c] = self.expectedStates[c] + self.hiddenNodeStates[c]                
                
            # M-step of EM.  We only update parameters after a collecting a certain number of samples
            if sampleCount >= self.nSamples:
                sampleCount = 0
                 # take expectation of sample states
                self.expectedStates = map(lambda x: x / self.nSamples, self.expectedStates)
                self._updteParams(self.alpha)
                
                likelihood = self.calcEvidenceLikelihood()
                self.likelihood.append(likelihood)   
                print "nIter: " + str(nIter + 1) + "; log likelihood of evidence: " + str(likelihood)                    

                # collect the current best fit models
                if likelihood > optLikelihood:
                    optLikelihood = likelihood
                    try:
                        cPickle.dump(self, open(self.pickleDumpFile, 'wb'))
                    except: 
                        raise Exception("Cannot create pickle dumpfile " + self.pickleDumpFile)

                bConverged = self._checkConvergence()
                if bConverged:
                    print "EM converged!"
                    break
                
                for c in range(self.nChains):  # clear expectedStates
                    self.expectedStates[c] = np.zeros(np.shape(self.hiddenNodeStates))


            
    def _checkConvergence(self):
        # To do, add convergence checking code
        if len(self.likelihood) < 20:
            return False
            
        ml = np.mean(self.likelihood[-5:-1])
        ratio = abs(self.likelihood[-1] - ml ) / abs(ml)        
        return ratio <= 0.001

        
    def _updateStates(self):
        nCases, antibody = np.shape(self.data)
        nCases, nHiddenNodes = np.shape(self.hiddenNodeStates[0])
        
        # interate through all nodes. 
        for c in range(self.nChains):
            for nodeId in self.network:
                # skip observed nodes
                if self.network.node[nodeId]['nodeObj'].bMeasured:
                    #check if observed variables has missing values
                    continue
                
                curNodeMarginal = self.calcNodeCondProb(nodeId, c)
                
                # sample states of current node based on the prob, and update 
                sampleState = np.zeros(nCases)
                sampleState[curNodeMarginal >= np.random.rand(nCases)] = 1.
                curNodeIndx = self.hiddenNodeStates[c].findColIndices(nodeId)
                self.hiddenNodeStates[c].data[:, curNodeIndx] = sampleState
                
                # clamp the activationState of perturbed nodes to a fix value
                if nodeId in self.perturbData:
                    for caseName, value in self.perturbData[nodeId]:
                        self.hiddenNodeStates[c].setCellValue(caseName, nodeId, value)
                
            # To do:  need to represent perturbation as nodes in network, 

            
    def calcNodeCondProb(self, nodeId, c):
        """
        Calculate the marginal probability of a node's state set to "1" conditioning 
        on all evidence.
        
        args:
             nodeId   A string id of the node of interest
             c        An integer indicate the chain from which the parameter 
                         vector to be used  
        """
        nodeObj = self.network.node[nodeId]['nodeObj']
        if nodeObj.type == "fluorescence":
            raise Exception("Call _caclNodeMarginalProb on an observed variable " + nodeId)

        nCases, nAntibody = np.shape(self.data)        

        # collect the state of the predecessors of the node
        preds = self.network.predecessors(nodeId)        
        logProbOneCondOnParents = 0
        logProbZeroCondOnParents = 0
        if len(preds) > 0:  # if the node has parents  
            # calculate p(curNode = 1 | parents);                 
            nodeParams = nodeObj.params[c,:] 
            predStates =  np.column_stack((np.ones(nCases), self.hiddenNodeStates.getValuesByCol(preds))) 
            pOneCondOnParents = 1 / (1 + np.exp( - np.dot(predStates, nodeParams)))
            pOneCondOnParents[pOneCondOnParents == 1] -= np.finfo(np.float).eps
            pOneCondOnParents[pOneCondOnParents == 0] += np.finfo(np.float).eps
            logProbOneCondOnParents  = np.log(pOneCondOnParents)
            logProbZeroCondOnParents = np.log(1 - pOneCondOnParents)

        # collect  evidence from  children 
        logProbChildCondOne = 0  # the prob of child conditioning on current node == 1
        logProdOfChildCondZeros = 0
        
        children = self.network.successors(nodeId)
        if len(children) > 0:
            for child in children:   
                nodeObj = self.network.node[child]['nodeObj']
                if nodeObj.type == "fluorescence":
                    curChildData = self.data.getValuesByCol(child)  
                    # calculate the probability using mixture Gaussian
                    logProbChildCondOne +=  (- np.log( nodeObj.sigma[c, 1]) + 0.5 * np.square(curChildData - nodeObj.mus[c, 1]) / np.square(nodeObj.sigma[c,1])) 
                    logProdOfChildCondZeros += (- np.log( nodeObj.sigma[c, 0]) + 0.5 * np.square(curChildData - nodeObj.mus[c, 0]) / np.square(nodeObj.sigma[c,0])) 
                    
                else:  # current child is a latent variable
                    # collect data and parameters associated with the node
                    curChildStates = self.hiddenStates[c].getValuesByCol(child)                    
                    
                    # Collect states of the predecessors of the child
                    childPreds = self.network.predecessors(child)
                    childNodeParams = nodeObj.params[c,:]
                    childPredStates = self.hiddenNodeStates[c].getValuesByCol(childPreds)
                    childPredStates = np.column_stack((np.ones(nCases), childPredStates)) # padding data with a column ones as bias
                    
                    # Set the state of current node to ones 
                    curNodePosInPredList = childPreds.index(nodeId) + 1 # offset by 1 because padding 
                    childPredStates[:, curNodePosInPredList] = np.ones(nCases)                
                    pChildCondCurNodeOnes = 1 / (1 + np.exp(-np.dot(childPredStates, childNodeParams)))
                    pChildCondCurNodeOnes[pChildCondCurNodeOnes==1] -= np.finfo(np.float).eps
                    pChildCondCurNodeOnes[pChildCondCurNodeOnes==0] += np.finfo(np.float).eps
                    logProbChildCondOne += np.log (curChildStates * pChildCondCurNodeOnes + (1 - curChildStates) * (1 - pChildCondCurNodeOnes))
                    
                    # set the state of the current node (nodeId) to zeros 
                    childPredStates [:, curNodePosInPredList] = np.zeros(nCases)
                    pChildCondCurNodeZeros = 1 / (1 + np.exp(- np.dot(childPredStates, childNodeParams))) 
                    pChildCondCurNodeZeros[pChildCondCurNodeZeros==1]  -= np.finfo(np.float).eps
                    pChildCondCurNodeZeros[pChildCondCurNodeZeros==0]  += np.finfo(np.float).eps
                    logProdOfChildCondZeros += np.log(curChildStates * pChildCondCurNodeZeros + (1 - curChildStates) * (1 - pChildCondCurNodeZeros))

        # now we can calculate the marginal probability of current node 
        curNodeMarginal = 1 / (1 + np.exp(logProbZeroCondOnParents + logProdOfChildCondZeros - logProbOneCondOnParents - logProbChildCondOne))
        return curNodeMarginal
    

    def parseGlmnetCoef(self, glmnet_res):        
        """ Parse the 'beta' matrix returned by calling glmnet through RPy2.
            Return the first column of 'beta' matrix of the glmnet object 
            with 3 or more non-zero values 
            """
        # read in intercept; a vector of length of nLambda
        a0 = np.array(glmnet_res.rx('a0'))[0]
        
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
                            
        return a0, betaMatrix       
      

      
        
    def _updteParams(self, alpha = 0.1):
        # Update the parameter associated with each node, p(n | Pa(n)) using logistic regression,
        # using expected states of precessors as X and current node states acrss samples as y
        nCases, nVariables = np.shape(self.data)
        for nodeId in self.network:
            preds = self.network.predecessors(nodeId)
            predIndices = self.hiddenNodeStates[0].findColIndices(preds)
            if len(preds) == 0:
                continue
            
            nodeObj = self.network.node[nodeId]['nodeObj']            
            nodeObj.fitResults = []
            for c in range(self.nChains): 
                expectedPredState = self.expectedStates[c][:, predIndices]
                curNodeData = self.data.getValuesByCol(nodeId)
                if self.network.node[nodeId]['nodeObj'].type == "fluorescence":
                    nodeObj.mus[c,0] = np.mean ((1-expectedPredState) * curNodeData)
                    nodeObj.sigma[c, 0] = np.std ((1-expectedPredState) * curNodeData)
                    nodeObj.mus[c, 1] = np.mean (expectedPredState * curNodeData)
                    nodeObj.sigmas[c, 1] = np.std (expectedPredState * curNodeData)                
                else:   
                    #x = np.column_stack((np.ones(nCases), expectedPredState))
                    x =  expectedPredState
                    y = self.hiddenNodeStates.getValuesByCol(nodeId) 
                    
                    #check if all x and y are of same value, which will lead to problem for glmnet
                    rIndx = map(lambda z: int(math.floor(z)), np.random.rand(3) * nCases)
                    if sum(y) == nCases:  # if every y == 1                      
                        y[rIndx] = 0                        
                    elif sum( map(lambda x: 1 - x, y)) == nCases:
                        y[rIndx] = 1        
                    y = robjects.vectors.IntVector(y)
                
                    allRwoSumOnes = np.where(np.sum(x, 0) == nCases)[0]
                    for c in allRwoSumOnes:
                        rIndx = map(lambda z: int(math.floor(z)), np.random.rand(3) * nCases)
                        x[rIndx, c] = 0 
                    allZeros = np.where(np.sum(np.ones(np.shape(x)) - x, 0) == nCases) 
                    for c in allZeros[0]:
                        rIndx = map(lambda z: int(math.floor(z)), np.random.rand(3) * nCases)
                        x[rIndx, c] = 1
                    
                    # call logistic regression using glmnet from Rpy
                    fit = glmnet (x, y, alpha = alpha, family = "binomial")
                    
                    # extract coefficients from Rpy2 vector object
                    a0, betaMatrix = self.parseGlmnetCoef(fit)
                    for j in range(np.shape(betaMatrix)[1]):
                        if sum(betaMatrix[:, j] != 0) > self.nParents:
                            break
                    if j >= len(a0):
                        j = len(a0) - 1
                    nodeObj.params[c,:] = betaMatrix[:, j]
                    nodeObj.fitResults.append(fit)


  
    def toXGMML(self):
        pass
    
    
    def trimEdgeByConsensus(self, percent):
        antibodies = [n for n in self.network if self.network[n]['nodeObj'].type == 'phosState']
        for node in antibodies:            
            preds = self.network.predecessors(node)
            if len(preds) > 0:
                nodeParams = self.network.node[node]['nodeObj'].params 
                nChain, nParams = np.shape(nodeParams)
            
                for i in range(1, nParams):
                    if sum(nodeParams[:,i]==0) > math.floor(nChain * percent):
                        self.network.remove_edge(preds[i-1], node)
            

    
    def trimLeaves(self):
        while True: 
            nodesToRemove = list()
            for node in self.network:
                children = self.network.successors(node)
                if len(children) == 0 and not self.network.node[node]['nodeObj'].bMeasured:
                    nodesToRemove.append(node)

            if len(nodesToRemove)  == 0:
                break

            for node in nodesToRemove:
                self.network.remove_node(node)
                
    
    
    
