import numpy as np

class Site: # This class represents individual sites in the lattice
    st = False
    ln = None
    rn = None
    rb = None
    
    def __init__(self, ud: bool, l: 'Site'=None, r: 'Site'=None, b: bool=True):
        self.st = ud
        self.ln = l # Left neighbor
        self.rn = r # Right neighbor
        self.rb = b # Bond with right neighbor
    
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
    
    def mb(self):
        self.rb = not self.rb