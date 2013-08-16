

class SigNetNode:
    """ An object represetation of nodes in a signaling network"""

    ## constructor
    # @param    name  String representation of the node
    # @param    nodeType  String represetation of the type of a node.  Possible
    # @param    bMeasured 
    def __init__(self, name, nodeType, bMeasured):
        self.name = name
        self.nodeType = nodeType
        self.bMeasured = bMeasured

    def isMeasured(self):
        return bMeasured

