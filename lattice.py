import numpy as np
from numpy import *
from sites import Site
from bonds import Bond

class Lattice:
    lat = np.empty(0, dtype=Site)
    J = 1
    B = 0
    E = 0
    dE = 0
    size = 0
    temp = 1
    K = 1

    def __init__(self, size: int=100, intEn: float=1, demonEn: float=2, magF: float=0, dir: bool=None): # Initialises lattice
        self.J = intEn # Interaction energy (keep this at 1)
        self.B = magF # Magnetic field
        self.K = self.J/self.temp # "temp" is technically actually kT in this case; merged for convenience
        self.size = size
        self.lat = np.empty(size, dtype=Site)
        altSpin = True
        for s in range(self.size):
            spin = True
            if dir is not None: spin = dir
            else: spin = altSpin
            if s%2 == 1: altSpin = not altSpin
            self.lat[s] = Site(ud=spin, ind=s, J=self.J, B=self.B)
            self.lat[s].rb = Bond()
        for s in range(self.size):
            self.lat[s].setNs(self.lat[s-1], self.lat[(s+1)%self.size])
            self.lat[s].lb = self.lat[s-1].rb
        self.E = self.energy()
        self.dE = demonEn
    
    def energy(self): # Returns instantaneous energy. Not for use inside loops (since multi-layer loops are slow)
        tot = 0.0
        for s in self.lat:
            tot += -self.J*s*s.rn*s.rb - self.B*s
        return tot
    
    def spin(self): # Returns overall spin direction
        return np.mean(self.lat)
    
    def metropolis(self, ind: int=1, temp: float=1): # Implementation of the Metropolis algorithm
        b = 1/temp # Measured relative to the Boltzmann constant, i.e. k = 1
        site = self.lat[ind]
        diffE = 2*(self.J*site*(site.rn + site.ln) + self.B*site) # Calculates difference in pre- and post-flip energy
        # To find probability of the flip occurring:
        doFlip = 1
        if diffE > 0:
            probOfFlip = np.e**(-b*(diffE))
            doFlip = np.random.choice([0, 1], p=[1-probOfFlip, probOfFlip])
        if doFlip:
            site.flip()
            self.E += diffE
        return self.lat
    
    def metroDemon(self, site: Site, brk: bool=True, rev: bool=False):
        order = [self.flip, self.bonds]
        if not brk: order.pop()
        if rev: order = np.flip(order)
        for f in order: f(site)
        return self.lat

    def flip(self, site: Site):
        diffE = -2*site.E()
        if diffE <= self.dE:
            self.trans(diffE)
            site.flip()
    
    def bonds(self, site: Site):
        initE = site.E()
        site.rb.flip()
        diffE = site.E() - initE
        if diffE <= self.dE: self.trans(diffE)
        else: site.rb.flip()
    
    def trans(self, E):
        self.E += E
        self.dE -= E
    
    def enSp(self): # Theoretical spin energy
        return -(self.size-1)*self.J*tanh(self.K)
    
    def entrSp(self): # Theoretical spin entropy
        return (self.enSp()/self.temp + (self.size*log(2) + (self.size - 1)*log(cosh(self.K))))/self.size
    
    def enDe(self): # Theoretical demon energy
        return self.J/(exp(self.K) - 1)
    
    def entrDe(self): # Theoretical demon entropy
        return self.enDe()/self.temp - log(1 - exp(-self.K))
    
    def entr(self, dist): # Calculate entropy of given distribution
        dist = dist[dist != 0]/sum(dist)
        ent = -np.sum(dist*log(dist))
        return dist, ent

    def __repr__(self):
        rv = ""
        for s in self.lat:
            rv += str(s)
            if s.rb:
                rv += " "
            else:
                rv += "|"
        return rv