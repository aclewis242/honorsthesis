import numpy as np
import sites as st

def randState():
    return np.random.choice(np.array([True, False]))

class Lattice:
    lat = np.empty(0, dtype=st.Site)
    J = 1
    B = 0
    E = 0
    dE = 0
    size = 0

    def __init__(self, size: int=20, intEn: float=1, demonEn: float=5, magF: float=0, dir: bool=None): # Initialises lattice
        self.J = intEn # Interaction energy (keep this at 1)
        self.B = magF # Magnetic field
        self.lat = np.empty(size, dtype=st.Site)
        self.size = size
        for s in range(self.size):
            spin = True
            if dir is not None: spin = dir
            else: spin = randState()
            self.lat[s] = st.Site(ud=spin)
        for s in range(self.size):
            self.lat[s].setNs(self.lat[s-1], self.lat[(s+1)%self.size])
        self.E = self.energy()
        self.dE = demonEn
    
    def energy(self): # Returns instantaneous energy. Not for use inside loops (since multi-layer loops are slow)
        tot = 0.0
        for s in self.lat:
            tot += -self.J*s*s.rn - self.B*s
        return tot
    
    def spin(self): # Returns overall spin direction
        return np.mean(self.lat)
    
    def metropolis(self, temp: float=1): # Implementation of the Metropolis algorithm
        b = 1/temp # Measured relative to the Boltzmann constant, i.e. k = 1
        siteInd = np.random.choice(range(self.size)) # Selects random site
        site = self.lat[siteInd]
        diffE = 2*(self.J*site*(site.rn + site.ln) + self.B*site) # Calculates difference in pre- and post-flip energy
        # To find probability of the flip occurring:
        doFlip = 1
        if diffE > 0:
            probOfFlip = np.e**(-b*(diffE))
            doFlip = np.random.choice([0, 1], p=[1-probOfFlip, probOfFlip])
        if doFlip:
            self.lat[siteInd].flip()
            self.E += diffE
        return self.lat

    def metroDemon(self): # Implementation of the Metropolis algorithm, demon ver. (Kreutz)
        siteInd = np.random.choice(range(self.size)) # Selects random site
        site = self.lat[siteInd]
        diffE = 2*(self.J*site*(site.rn + site.ln) + self.B*site) # Calculates difference in pre- and post-flip energy
        newDE = self.dE - diffE
        if newDE > 0:
            self.lat[siteInd].flip()
            self.dE = newDE
            self.E += diffE
        return self.lat
    
    def sweep(self, alg=metropolis):
        for i in range(self.size):
            self.lat = alg()
        return self.E/self.size