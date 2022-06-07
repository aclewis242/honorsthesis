import numpy as np
import state as st

def buildLattice(size: int=20):
    lattice = np.empty(size, dtype=st.State)
    for s in range(lattice.size):
        lattice[s] = st.State()
    return lattice

if __name__ == "__main__":
    lat = buildLattice()