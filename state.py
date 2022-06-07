import numpy as np

class State:
    st = False
    ln = None
    rn = None
    
    def __init__(self, ud: bool, l: 'State'=None, r: 'State'=None):
        self.st = ud
        self.ln = l
        self.rn = r
    
    def __repr__(self):
        return str(self.st)
    
    def getNum(self):
        if self.st: return 1
        else: return -1
    
    def flip(self):
        self.st = ~self.st
    
    def setLN(self, l: 'State'):
        self.ln = l
    
    def setRN(self, r: 'State'):
        self.rn = r
    
    def setNs(self, l: 'State', r: 'State'):
        self.setLN(l)
        self.setRN(r)