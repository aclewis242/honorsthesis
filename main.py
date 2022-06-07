import numpy as np
import state as st

def randState():
    return np.random.choice(np.array([True, False]))

def buildLattice(size: int=20):
    lattice = np.empty(size, dtype=st.State)
    for s in range(lattice.size):
        lattice[s] = st.State(ud=randState())
    for s in range(lattice.size):
        lattice[s].setNs(lattice[s-1], lattice[(s+1)%20])
    return lattice

if __name__ == "__main__":
    lat = buildLattice()