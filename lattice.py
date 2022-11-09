import numpy as np
import sites as st
from bonds import Bond
from math import exp
from mpmath import nsum, inf

class Lattice:
    lat = np.empty(0, dtype=st.Site)
    J = 1
    B = 0
    E = 0
    dE = 0
    size = 0
    temp = 1

    def __init__(self, size: int=100, intEn: float=1, demonEn: float=2, magF: float=0, dir: bool=None): # Initialises lattice
        self.J = intEn # Interaction energy (keep this at 1)
        self.B = magF # Magnetic field
        self.size = size
        self.lat = np.empty(size, dtype=st.Site)
        altSpin = True
        for s in range(self.size):
            spin = True
            if dir is not None: spin = dir
            else: spin = altSpin
            if s%2 == 1: altSpin = not altSpin
            self.lat[s] = st.Site(ud=spin)
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

    def metroDemon(self, ind: int=1): # Implementation of the Metropolis algorithm, demon ver. (Creutz)
        site = self.lat[ind]
        diffE = (site.rb + site.lb)*self.J*site*(site.rn + site.ln) + self.B*site # Calculates difference in pre- and post-flip energy
        if diffE < 0 or self.dE - diffE >= -2:
            site.flip()
            if self.dE - diffE == -2:
                print("diffE: {}".format(diffE))
                diffE/=2
                site.rb.brk()
            else:
                site.mkall()
            self.dE -= diffE
            self.E += diffE
            # if self.E != self.energy():
            #     print("bruh")
            #     print("self.E: {}".format(self.E))
            #     print(site.rb)
            #     print("self.energy(): {}".format(self.energy()))
            #     print(self.lat)
        return self.lat
    
    def bf(self, n): # Boltzmann
        return exp(-n*self.J/self.temp) # k is treated as 1 for our purposes
    
    def entr(self): # Theoretical entropy
        E = -nsum(lambda n: n*self.J*self.bf(n), [0, inf])/nsum(lambda n: self.bf(n), [0, inf])
        A = -np.log(float(nsum(lambda n: self.bf(n), [0, inf])))
        return (E - A)/self.temp
    
    def __repr__(self):
        rv = ""
        for s in self.lat:
            rv += s
            if s.rb:
                rv += " "
            else:
                rv += "|"
        return rv