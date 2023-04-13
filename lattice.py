import numpy as np
from sites import Site
from bonds import Bond

class Lattice:
    '''
    The (1D) lattice of spins. This is the core of the model.

    ### Attributes
    lat: The array of sites
    J: Interaction energy
    B: Magnetic field strength
    E: System energy
    dE: Demon energy
    size: Lattice size, i.e. length in this case
    
    (NB: J and B are not meant to be changed from their default values!)
    '''
    lat = np.empty(0, dtype=Site)
    J = 1
    B = 0
    E = 0
    dE = 0
    size = 0

    def __init__(self, size: int=500, intEn: float=1, demonEn: float=0, magF: float=0, dir: bool=None, align: float=None):
        '''
        Initialises the lattice.

        ### Parameters
        size: Lattice size
        intEn: Interaction energy
        demonEn: Initial demon energy
        magF: Magnetic field strength
        dir: Initial direction of all spins (`True` is up, `False` is down)
        align: Ratio of up- to down-spins
        '''
        self.J = intEn
        self.B = magF
        self.size = size
        self.lat = np.empty(size, dtype=Site)
        altSpin = True
        for s in range(self.size): # Creates the initial set of spins & bonds
            spin = True
            if dir is not None: spin = dir
            else: spin = altSpin
            if s%2 == 1: altSpin = not altSpin
            self.lat[s] = Site(ud=spin, ind=s, J=self.J, B=self.B)
            self.lat[s].rb = Bond()
        for s in range(self.size): # Joins the spins and bonds together properly
            self.lat[s].setNs(self.lat[s-1], self.lat[(s+1)%self.size])
            self.lat[s].lb = self.lat[s-1].rb
        if align is not None: # Generate directions based on 'align'
            align = align%1
            for s in range(self.size): self.lat[s].st = np.random.choice([False, True], p=[1-align, align])
        self.E = self.energy
        self.dE = demonEn
    
    @property
    def energy(self) -> int:
        '''
        Calculates the instantaneous energy of the system.
        (NB: Avoid using inside loops! The lattice can get quite large, so running this frequently will make the simulation take much longer.
        Store and update the energy dynamically instead!)

        ### Returns
        tot: Total energy of system
        '''
        tot = 0
        for s in self.lat: tot += -self.J*s*s.rn*s.rb - self.B*s
        return tot
    
    def metropolis(self, site: Site, temp: float=5, *args, **kwargs):
        '''
        Decides whether or not to update a particular site according to the Metropolis algorithm.

        ### Parameters
        site: The site to be checked
        '''
        # 'temp' is measured relative to the Boltzmann constant, i.e. kT = temp
        b = 1/temp
        diffE = 2*(self.J*site*(site.rn + site.ln) + self.B*site) # Calculates difference in pre- and post-flip energy
        doFlip = True
        if diffE > 0:
            probOfFlip = np.e**(-b*(diffE))
            doFlip = np.random.choice([False, True], p=[1-probOfFlip, probOfFlip])
        if doFlip:
            # site.flip()
            self.trans(diffE)
            if abs(self.E) >= self.size: self.trans(-diffE)
    
    def demon(self, site: Site, brk: bool=True, rev: bool=False, **kwargs):
        '''
        Decides whether or not to update a particular site according to the Creutz algorithm.

        ### Parameters
        site: The site to be checked
        brk: Allow/disallow bond breaking/making
        rev: Whether or not the system is currently reversing
        '''
        order = [self.flip, self.bonds]
        if not brk: order.pop()
        if rev: order = np.flip(order)
        [f(site) for f in order]

    def flip(self, site: Site):
        '''
        Decides whether or not to flip the site. Will always occur if the demons have enough energy for it.

        ### Parameters
        site: The site to be checked
        '''
        diffE = -2*site.E
        if diffE <= self.dE:
            self.trans(diffE)
            site.flip()
    
    def bonds(self, site: Site):
        '''
        Decides whether or not to break the site's (right) bond. Will always occur if the demons have enough energy for it.

        ### Parameters
        site: The site to be checked
        '''
        initE = site.E
        site.rb.flip()
        diffE = site.E - initE
        if diffE <= self.dE: self.trans(diffE)
        else: site.rb.flip()
    
    def trans(self, E: int):
        '''
        Transfers energy from the demons to the system.

        ### Parameters
        E: The quantity of energy to be transferred over
        '''
        self.E += E
        self.dE -= E

    def __repr__(self):
        '''
        Represents the lattice as a string. '|' indicates a broken bond.
        (NB: Used only when printing to the console - this is not a general string method!)

        ### Returns
        rv: The string representation of the lattice
        '''
        rv = ''
        for s in self.lat:
            rv += str(s)
            if s.rb: rv += ' '
            else: rv += '|'
        return rv