from bonds import Bond

class Site: # This class represents individual sites in the lattice
    st = False
    ln = None
    lb = None
    rn = None
    rb = None
    
    def __init__(self, ud: bool, ln: 'Site'=None, lb: Bond=None, rn: 'Site'=None, rb: Bond=None):
        self.st = ud
        self.ln = ln # Left neighbor
        self.lb = lb # Left bond
        self.rn = rn # Right neighbor
        self.rb = rb # Right bond
    
    def flip(self):
        self.st = not self.st
    
    def setNs(self, l: 'Site', r: 'Site'):
        self.ln = l
        self.rn = r
    
    def mkall(self):
        self.lb.mk()
        self.rb.mk()
    
    def brkall(self):
        self.lb.brk()
        self.rb.brk()
    
    def __repr__(self):
        return str(self)
    
    def __int__(self):
        if self.st: return 1
        else: return -1
    
    def __add__(self, s2):
        return int(self) + int(s2)
    
    __radd__ = __add__
    
    def __mul__(self, s2):
        return int(self) * int(s2)
    
    __rmul__ = __mul__

    def __str__(self):
        if self.st: return "↑"
        else: return "↓"