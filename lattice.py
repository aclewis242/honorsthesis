import numpy as np
import sites as st

def randState():
    return np.random.choice(np.array([True, False]))

class Lattice:
    lat = np.empty(0, dtype=st.Site)
    J = 1
    B = 0
    E = 0
    size = 0

    def __init__(self, size: int=20, intEn: float=1, magF: float=0, dir: bool=None):
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
    
    def energy(self):
        tot = 0.0
        for s in self.lat:
            tot += -self.J*s*s.rn - self.B*s
        return tot
    
    def spin(self):
        return np.mean(self.lat)
    
    def metropolis(self, energy, temp: float=1):
        b = 1/temp # Measured relative to the Boltzmann constant, i.e. k = 1
        siteInd = np.random.choice(range(self.size))
        site = self.lat[siteInd]
        newE = energy + 2*(self.J*site*(site.rn + site.ln) + self.B*site)
        doFlip = 1
        if newE > energy:
            probOfFlip = np.e**(-b*(newE - energy))
            doFlip = np.random.choice([0, 1], p=[1-probOfFlip, probOfFlip])
        if doFlip:
            self.lat[siteInd].flip()
            energy = newE
        return self.lat, energy
    
    def sweep(self):
        for i in range(self.size):
            self.lat, self.E = self.metropolis(self.E)
        return self.E/self.size