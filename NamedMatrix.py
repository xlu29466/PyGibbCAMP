# -*- coding: utf-8 -*-
"""
A wrapper class enable access data matrix elements by col and row names

Created on Sun Aug 25 08:40:33 2013



@author: xinghualu
"""

import numpy as np
from StringIO import StringIO

class NamedMatrix:
    ## Constructor
    #  @param  filename=None  A string point to a text matrix file
    #  @param delimiter=','  A string indicate the delimiter separating fields in txt
    #  @param npMatrix=None  A reference to a numpy matrix
    #  @colnames  A string array of column names
    #  @rownames  A string array of rownames

    def __init__(self, filename = None, delimiter = ',', npMatrix = None, colnames = None, rownames = None):
        
        if filename and npMatrix:  
            raise Exception ("Cannot create a NamedMatrix with both 'npMatrix' and 'filename' arguments set")
        if not filename and not npMatrix:
            raise Exception ("Attempt to create a NameMatrix without 'filename' or an 'npMatrix'")
        
        if filename:
            try:
                f = open(filename, 'r')
                lines = f.readlines()
            except IOError:
                print "Fail to read  file " + filename
            
            colnames = None
            if len(lines) == 1:  # Mac version csv, with "\r" as return
                lines = lines[0].split("\r")
                colnames = lines.pop(0).split(',') # split  header and extract colnames
                map(lambda x: x.rstrip(), lines)  # remove the "\r"
                lines = "\n".join(lines)  # use "\n" to join lines
            else:
                colnames = lines.pop(0).split(',')
                lines = "".join(lines)
                
            # extract condition name
            self.rownames = list()
            for l in lines():
                rownames.append(l.split(',')[0])         
                            
            # read in data and generate a numpy data matrix
            self.data = np.genfromtxt(StringIO(lines), delimiter = ",", usecol=tuple(range(1, len(colnames))))

            
        if npMatrix:
            self.data = npMatrix
            nrow, ncol = np.shape(self.data)
            if colnames:
                if len(colnames) == ncol:
                    self.colnames = colnames
                else:
                    raise Exception("Dimensions of input colnames and matrix do not agree")
            else:
                self.colnames = list()
                for c in range(ncol):
                    self.colnames.append('c' + str(c))
            if rownames:
                if len(rownames) == nrow:
                    self.rownames = rownames
                else:
                    raise Exception("Dimensions of input rownames and matrix do not agree")
            else:
                self.rownames = list()
                for r in range(nrow):
                    self.rownames.append('r' + str(r))
                    
    def setColnames(self, colnames):
        if len(colnames) == len(self.colnames):
            self.colnames = colnames
        else:
            raise Exception("New colnames vector has differnt dimension as the original colnames")
            
    def getColnames(self):
        return self.colnames
            
    def setRownames(self, rownames):
        if len(rownames) == len(self.rownames):
            self.rownames = rownames
        else:
            raise Exception("New rownames vector has differnt dimension as the original colnames")
    
    def getRownames(self):
        return self.rownames
            
    def getValuesByCol(self, colnames):
        if set(colnames) - set(self.colnames):
            raise Exception("Try to access nonexisting column")
        colIndx = [lambda x: self.colnames.index(x) for x in colnames if x in self.colnames]
        return self.data[:, colIndx]
        
    def setValuesByColName(self, values, col):      
        self.data[:,self.colnames.index(col)] = values
        
    def getValuesByRow(self, rownames):
        if set(rownames) - set(self.rownames):
            raise Exception("Try to access non-existing rows")
        rowIndx = [lambda x: self.colnames.index(x) for x in rownames if x in self.rownames]
        return self.data[rowIndx, :]
        
     
    def shape(self):
        if self.data:
            return np.shape(self.data)
            
        else:
            return None
            
    ## Return the position indices of colnames  
    def findColIndices(self, colnames):
        if set(colnames) - set(self.colnames):
            raise Exception("Unknown column name is used to query index")
            
        return [lambda x: self.colnames.index(x) for x in colnames]
        
    ## Return the position indices of rownames 
    def findRowIndices(self, rownames):
        if set(rownames) - set(self.rownames):
            raise Exception("Unknown column name is used to query index")
            
        return [lambda x: self.rownames.index(x) for x in rownames]
        
    
    def setCellValue(self, rowname, colname, value):
        value = np.float(value) # force it into a np.float
        self.data[self.rownames.index(rowname), self.colnames.index(colname)] = value
            
        
    
    
    
        
        
            
            
            
            
                        
  
                    
