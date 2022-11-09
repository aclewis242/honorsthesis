class Bond:
    bd = None

    def __init__(self, bu: bool=True):
        self.bd = bu
    
    def mk(self):
        self.bd = True
    
    def brk(self):
        self.bd = False
    
    def __int__(self):
        return self.bd
    
    def __add__(self, b2):
        return int(self) + int(b2)
    
    __radd__ = __add__
    
    def __mul__(self, b2):
        return int(self) * int(b2)
    
    __rmul__ = __mul__

    def __bool__(self):
        return self.bd
    
    def __repr__(self):
        return str(self.bd)