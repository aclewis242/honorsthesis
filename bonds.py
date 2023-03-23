class Bond:
    '''
    Represents the bond between two sites in the lattice.

    ### Attributes
    bd: Whether or not the bond exists
    '''
    bd = None

    def __init__(self, bu: bool=True):
        '''
        Initialises the bond.

        ### Parameters
        bu: Whether or not the bond exists
        '''
        self.bd = bu
    
    def mk(self):
        '''
        Makes the bond.
        '''
        self.bd = True
    
    def brk(self):
        '''
        Breaks the bond.
        '''
        self.bd = False
    
    def flip(self):
        '''
        Makes a broken bond or breaks a made bond.
        '''
        self.bd = not self.bd
    
    def __int__(self):
        '''
        Represents the bond as an integer.
        '''
        return int(self.bd)
    
    def __add__(self, b2):
        '''
        Allows the '+' operator to be used with bond objects.

        ### Parameters
        b2: The subject of the addition. Can be any type capable of being cast as an integer
        '''
        return int(self) + int(b2)
    
    __radd__ = __add__ # (Allows addition from the other direction, i.e. 5 + b instead of b + 5)
    
    def __mul__(self, b2):
        '''
        Allows the '*' operator to be used with bond objects.

        ### Parameters
        b2: The subject of the multiplication. Can be any type capable of being cast as an integer
        '''
        return int(self) * int(b2)
    
    __rmul__ = __mul__ # (Allows multiplication from the other direction, i.e. 5*b instead of b*5)

    def __bool__(self) -> bool:
        '''
        Represents the bond as a boolean.
        '''
        return self.bd
    
    def __repr__(self):
        '''
        Represents the bond as a string, specifically for console printing purposes.
        '''
        return str(self.bd)