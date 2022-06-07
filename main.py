import numpy as np
import state

def randState():
    return np.random.choice(np.array([True, False]))

def buildLattice(size: int=20):
    lattice = np.empty(size, dtype=bool)
    for s in range(lattice.size):
        lattice[s] = randState()
    return lattice

if __name__ == "__main__":
    lat = buildLattice()
    print(lat)
    flip(lat[0])
    print(lat)