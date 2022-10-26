import numpy as np
import sites as st

class Lattice:
    lat = np.empty(0, dtype=st.Site)
    J = 1
    B = 0
    E = 0
    dE = 0
    size = 0
    # order = np.empty(0)

    def __init__(self, size: int=100, intEn: float=1, demonEn: float=1, magF: float=0, dir: bool=None): # Initialises lattice
        self.J = intEn # Interaction energy (keep this at 1)
        self.B = magF # Magnetic field
        self.size = size
        self.lat = np.empty(size, dtype=st.Site)
        # self.order = np.array([14, 6, 9, 0, 8, 19, 2, 16, 12, 11, 3, 18, 10, 15, 7, 4, 1, 5, 17, 13])
        # self.order = np.linspace(0, size-1, size)
        # rand.shuffle(self.order)
        altSpin = True
        for s in range(self.size):
            spin = True
            if dir is not None: spin = dir
            else: spin = altSpin
            altSpin = not altSpin
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
            self.lat[ind].flip()
            self.E += diffE
        return self.lat

    def metroDemon(self, ind: int=1): # Implementation of the Metropolis algorithm, demon ver. (Creutz)
        site = self.lat[ind]
        diffE = 2*(self.J*site*(site.rn + site.ln) + self.B*site) # Calculates difference in pre- and post-flip energy
        newDE = self.dE - diffE
        if newDE > 1: # and (site*site.rn<0 or site*site.ln<0):
            self.lat[ind].flip()
            self.dE = newDE
            self.E += diffE
        return self.lat
    
    # def sweep(self, order, alg=metroDemon):
    #     for i in order:
    #         self.lat = alg(ind=int(i))
    #     return self.E/self.size
    
    # def printIt(self, alg=metroDemon):
    #     steps = range(self.size)
    #     Es = np.zeros_like(steps)
    #     for i in self.order:
    #         self.lat = alg(ind=int(i))
    #         Es[np.where(self.order==i)] = self.E
    #     plt.plot(steps, Es)
    #     plt.savefig("it.png")
    #     plt.close()