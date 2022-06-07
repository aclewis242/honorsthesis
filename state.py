import numpy as np

class State:
    st = False
    
    def __init__(self, ud: bool=np.random.choice(np.array([True, False]))):
        self.st = ud
    
    def __repr__(self):
        return str(self.st)
    
    def getNum(self):
        if self.st: return 1
        else: return -1
    
    def flip(self):
        self.st = ~self.st