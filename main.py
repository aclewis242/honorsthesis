import numpy as np
import sites as st
import lattice as lat

if __name__ == "__main__":
    l = lat.Lattice()
    print(l.lat)
    l.metropolis()