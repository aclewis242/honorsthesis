class state:
    st = False
    
    def __init__(self, ud: bool):
        self.st = ud
    
    def getNum(self):
        if self.st: return 1
        else: return -1
    
    def flip(self):
        self.st = ~self.st