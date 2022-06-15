import numpy as np
import lattice as lat
import matplotlib.pyplot as plt

if __name__ == "__main__":
    l = lat.Lattice()
    numSteps = 500
    steps = range(numSteps)
    energies = np.zeros_like(steps, dtype=float)
    energies = np.append(energies, 0.0)
    energies[0] = l.energy()/l.size()
    print(energies[0])
    for i in steps:
        energies[i+1] = l.sweep()
    steps = range(numSteps+1)
    
    plt.plot(steps, energies)
    plt.ylabel("Energy")
    plt.xlabel("Timestep")
    plt.savefig("energy.png")
    plt.show()