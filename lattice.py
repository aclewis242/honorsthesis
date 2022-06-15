import numpy as np
import sites as st

def randState():
    return np.random.choice(np.array([True, False]))

class Lattice:
    lat = np.empty(0, dtype=st.Site)
    J = 1
    B = 0

    def __init__(self, size: int=500, intEn: float=1, magF: float=0, dir: bool=None):
        self.J = intEn
        self.B = magF
        self.lat = np.empty(size, dtype=st.Site)
        for s in range(self.size()):
            spin = True
            if dir is not None: spin = dir
            else: spin = randState()
            self.lat[s] = st.Site(ud=spin)
        for s in range(self.size()):
            self.lat[s].setNs(self.lat[s-1], self.lat[(s+1)%self.size()])
    
    def size(self):
        return self.lat.size
    
    def energy(self):
        tot = 0.0
        for s in self.lat:
            tot += -self.J*s*s.rn - self.B*s
        return tot
    
    def metropolis(self, temp: float=200):
        b = 1/temp # Measured relative to the Boltzmann constant, i.e. k = 1
        siteInd = np.random.choice(range(self.size()))
        site = self.lat[siteInd]
        currE = self.energy()
        newE = currE + 2*(self.J*site*(site.rn + site.ln) + self.B*site)
        doFlip = 1
        if newE > currE:
            probOfFlip = np.e**(-b*(newE - currE))
            doFlip = np.random.choice([0, 1], p=[1-probOfFlip, probOfFlip])
        if doFlip:
            self.lat[siteInd].flip()
        return self.lat
    
    def sweep(self):
        for i in range(self.size()):
            self.lat = self.metropolis()
        return self.energy()/self.size()