from bonds import Bond

class Site:
    '''
    Represents an individual site within the lattice.

    ### Attributes
    st: The direction the site points. `True` is up, `False` is down
    ln: The site to its left ('left neighbor')
    lb: The bond connecting it to its left neighbor
    rn: The site to its right ('right neighbor')
    rb: The bond connecting it to its right neighbor
    ind: The index of the site within the lattice
    J: Interaction energy
    B: Magnetic field strength

    (NB: J and B are not meant to be changed from their default values!)
    '''
    st = False
    ln = None
    lb = None
    rn = None
    rb = None
    ind = None
    J = 1
    B = 0
    
    def __init__(self, ud: bool, ln: 'Site'=None, lb: Bond=None, rn: 'Site'=None, rb: Bond=None, ind: int=None, J: float=1, B: float=0):
        '''
        Initialises the site.

        ### Parameters
        ud: The direction the site points. `True` is up, `False` is down
        ln: The site to its left ('left neighbor')
        lb: The bond connecting it to its left neighbor
        rn: The site to its right ('right neighbor')
        rb: The bond connecting it to its right neighbor
        ind: The index of the site within the lattice
        J: Interaction energy
        B: Magnetic field strength
        '''
        self.st = ud
        self.ln = ln
        self.lb = lb
        self.rn = rn
        self.rb = rb
        self.ind = ind
        self.J = J
        self.B = B
    
    def flip(self):
        '''
        Flips the orientation of the site.
        '''
        self.st = not self.st
    
    def setNs(self, l: 'Site', r: 'Site'):
        '''
        Sets the site's neighbors.

        ### Parameters
        l: The site's future left neighbor
        r: The site's future right neighbor
        '''
        self.ln = l
        self.rn = r
    
    @property
    def E(self) -> int:
        '''
        Returns the energy of the site.
        '''
        return -(self.J*self*(self.rn*self.rb + self.ln*self.lb) + self.B*self)

    def __repr__(self):
        '''
        Represents the site as a string, specifically for console printing purposes.
        '''
        return str(self)
    
    def __int__(self) -> int:
        '''
        Represents the site as a number. Positive is up, negative is down.
        '''
        if self.st: return 1
        else: return -1
    
    def __add__(self, s2):
        '''
        Allows the '+' operator to be used with site objects.

        ### Parameters
        s2: The subject of the addition. Can be any type capable of being cast as an integer
        '''
        return int(self) + int(s2)
    
    __radd__ = __add__ # (Allows addition from the other direction, i.e. 5 + s instead of s + 5)
    
    def __mul__(self, s2):
        '''
        Allows the '*' operator to be used with site objects.

        ### Parameters
        s2: The subject of the multiplication. Can be any type capable of being cast as an integer
        '''
        return int(self) * int(s2)
    
    __rmul__ = __mul__ # (Allows multiplication from the other direction, i.e. 5*s instead of s*5)

    def __str__(self) -> str:
        '''
        Represents the site's direction with a single text character.
        '''
        if self.st: return '↑'
        else: return '↓'