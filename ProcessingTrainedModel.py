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
import networkx, os, sys
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


pickledir = "/Users/xinghualu/Dropbox/Dream8/PyGibbCAMP/ProcessedData/cellline.specific.tables/MCF7/"
#pickledir = "/Users/xinghualu/Dropbox/Dream8/PyGibbCAMP/ProcessedData/insilico/"
os.chdir(pickledir)
files = os.listdir(pickledir)
bestModel = None
bestLikelihood = float("-inf")
bestModelName = None
for f in files:    
    if re.search('pickle$', f):
        print "Reading pickle file " + f
        m = cPickle.load(open(f, 'rb'))

        if  m.likelihood[-1] > bestLikelihood:
            bestLikelihood = m.likelihood[-1]
            bestModel = m
            bestModelName = f
        
#f = pickledir + "data.matrix.csvalpha-1.0.pickle"        
print "Best human model: " + bestModelName

stimuli = ['EGF', 	'FGF1', 	'HGF',	'IGF1',	'Insulin',	'NRG1',	'PBS',	'Serum']
#stimuli = ['loLIG1',	'hiLIG1',	'loLIG2',	'hiLIG2']
for s in stimuli:
    print "Extract network for " + s
    tmpNet = bestModel.getStimuliSpecificNet(s)
    print "Total number of nodes: " + str(len(tmpNet.nodes())) + "; total number of edges: " + str(len(tmpNet.edges()))
    graphmlFile = bestModelName + "." + s + "edges.txt"
    f = open (graphmlFile, 'w')
    for u, v in tmpNet.edges():
        f.write(u + "\t" + tmpNet.edge[u][v]['effect'] + "\t" + v + "\n")
    f.close()
    
    draw_xgmml(tmpNet, bestModelName + "." + s)

        