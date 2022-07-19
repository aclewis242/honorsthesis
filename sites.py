import numpy as np

class Site: # This class represents individual sites in the lattice
    st = False
    ln = None
    rn = None
    
    def __init__(self, ud: bool, l: 'Site'=None, r: 'Site'=None):
        self.st = ud
        self.ln = l # Left neighbor
        self.rn = r # Right neighbor
    
    def __repr__(self):
        return str(self.st)
    
    def __int__(self):
        if self.st: return 1
        else: return -1
    
    def __add__(self, s2):
        return int(self) + int(s2)
    
    __radd__ = __add__
    
    def __mul__(self, s2):
        return int(self) * int(s2)
    
    __rmul__ = __mul__
    
    def flip(self):
        self.st = not self.st
    
    def setLN(self, l: 'Site'):
        self.ln = l
    
    def setRN(self, r: 'Site'):
        self.rn = r
    
    def setNs(self, l: 'Site', r: 'Site'):
        self.setLN(l)
        self.setRN(r)