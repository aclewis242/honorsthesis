from bonds import Bond

class Site: # This class represents individual sites in the lattice
    st = False
    ln = None
    lb = None
    rn = None
    rb = None
    ind = None
    J = 1
    B = 0
    
    def __init__(self, ud: bool, ln: 'Site'=None, lb: Bond=None, rn: 'Site'=None, rb: Bond=None, ind: int=None, J: float=1, B: float=0):
        self.st = ud # Direction
        self.ln = ln # Left neighbor
        self.lb = lb # Left bond
        self.rn = rn # Right neighbor
        self.rb = rb # Right bond
        self.ind = ind # Index of site
        self.J = J # Interaction energy
        self.B = B # Magnetic field strength
    
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
    
    def E(self):
        return -(self.J*self*(self.rn*self.rb + self.ln*self.lb) + self.B*self)
    
    def copy(self):
        return Site(ud=self.st, ln=self.ln, lb=self.lb, rn=self.rn, rb=self.rb, ind=self.ind, J=self.J, B=self.B)

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
        if self.st: return '↑'
        else: return '↓'