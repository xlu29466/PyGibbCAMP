# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 08:29:35 2013

@author: xinghua
"""

# this script process the trained PySBVSigNet models and identify a best model
# for rat and human respectively.  
import os, re, cPickle
from PyGibbCAMP import PyGibbCAMP
import networkx as nx
import os, sys
from graphNew2011_xgmml import *


def draw_xgmml(p_graph, p_graphName):
   g = Graph(name=p_graphName,directed=True)
   for edge in p_graph.edges():
      edgeColor = ''
      if p_graph[edge[0]][edge[1]]['effect']=='+':
         edgeColor += 'red'
      else:
         edgeColor += 'blue'
      n1 = Node(description='descreption for node 1', id=edge[0], label=edge[0], color='102,153,0')
      n2 = Node(description='descreption for node 1', id=edge[1], label=edge[1], color='102,153,0')
      edge = Edge(n1.id, n2.id, description='edge description', color=edgeColor,style = 'solid')

      g.add_node(n1)
      g.add_node(n2)
      g.add_edge(edge)
   g.drawGraph()


#pickledir = "/Users/xinghualu/Dropbox/Dream8/PyGibbCAMP/ProcessedData/cellline.specific.tables/BT20/"

cellLine = "UACC812"
pickledir = "/Users/xinghualu/Dropbox/Dream8/PyGibbCAMP/ProcessedData/cellline.specific.tables/" + cellLine + "/"
#pickledir = "/Users/xinghualu/Dropbox/Dream8/PyGibbCAMP/ProcessedData/insilico/"
os.chdir(pickledir)
files = os.listdir(pickledir)
bestModels = dict()
nKeeps = 5
n = 0
for f in files:    
    if re.search('pickle$', f):
        n += 1
        print "Reading pickle file " + f
        m = cPickle.load(open(f, 'rb'))

        if n <= nKeeps: 
            bestModels[f] = m
        else:
            for model in bestModels:
                if m.likelihood[-1] > bestModels[model].likelihood[-1]:                    
                    del (bestModels[model])
                    bestModels[f] = m
                    break
        
print "Best models: " + str(bestModels)

f = open("likelihoods.txt", 'w')
for m in bestModels:
    f.write(m + " likelihood: " + str(bestModels[m].likelihood[-1]) + "\n")
f.close()

stimuli = ['EGF', 'FGF1', 	'HGF',	'IGF1',	'Insulin',	'NRG1',	'PBS',	'Serum']
#stimuli = ['EGF']#, 	'FGF1', 	'HGF',	'IGF1',	'Insulin',	'NRG1',	'PBS',	'Serum']
#stimuli = ['loLIG1',	'hiLIG1',	'loLIG2',	'hiLIG2']

pthreshold = 0.2
for s in stimuli:
    outNetwork = nx.DiGraph()
    print "Extract network for " + s
    for model in bestModels:
        print "Processing model: " + model
        tmpNet = bestModels[model].getStimuliSpecificNet(s, pthreshold)
        f = open (model + "." + s + ".p" + str(pthreshold) + ".edges.txt", 'w')
        for u, v in tmpNet.edges():
            f.write( u + "\t" + tmpNet.edge[u][v]['effect'] + "\t"  + v + "\n")
            if not outNetwork.has_edge(u, v):
                outNetwork.add_edge(u, v, effect = [tmpNet.edge[u][v]['effect']])
            else:
                outNetwork.edge[u][v]['effect'].append(tmpNet.edge[u][v]['effect'])
        f.close()
        draw_xgmml(tmpNet, model + s + ".p" + str(pthreshold))
        
    fileName1 = pickledir + cellLine + "." + s + ".combinedNet.edges.txt"
    fileName2 = pickledir + cellLine + "." + s + ".combinedNet.edges.prob.txt"
    f1 = open (fileName1, 'w')
    f2 = open (fileName2, 'w')
    
    for u, v in outNetwork.edges():
        effect = 0
        prob = 0
        for e in outNetwork.edge[u][v]['effect']:
            if e == '+':
                effect += 1
            else:
                effect -=1
        avgEffect = effect / len(outNetwork.edge[u][v]['effect'])
        if avgEffect > 0: 
            f1.write(u + "\t1\t"  + v + "\n")
            outNetwork.edge[u][v]['effect']= '+'
        else:
            f1.write(u + "\t-1\t"  + v + "\n")
            outNetwork.edge[u][v]['effect']= '-'
        prob = float(len(outNetwork.edge[u][v]['effect'])) / len(bestModels)
        f2.write(u + "\t" + str(prob) + "\t"  + v + "\n")
    f1.close()
    f2.close()
    draw_xgmml(outNetwork, pickledir + cellLine + "." + s + ".combinedNet")
    

        