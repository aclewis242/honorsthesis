import numpy as np
import lattice as lat
import matplotlib.pyplot as plt

if __name__ == "__main__":
    l = lat.Lattice()
    numSteps = 500 # "Time"
    steps = range(numSteps)
    energies = np.zeros_like(steps, dtype=float)
    energies = np.append(energies, 0.0)
    demonEnergies = np.zeros_like(energies, dtype=float)
    spins = np.zeros_like(energies, dtype=float)
    energies[0] = l.energy()/l.size
    demonEnergies[0] = l.dE/l.size
    spins[0] = l.spin()
    for i in steps:
        energies[i+1] = l.sweep(alg=l.metroDemon)
        demonEnergies[i+1] = l.dE/l.size
        spins[i+1] = l.spin()
    steps = range(numSteps+1)
    
    # ENERGIES
    plt.plot(steps, energies, label='System energies')
    if demonEnergies[0] != demonEnergies[-1] or demonEnergies[0] != demonEnergies[1]:
        plt.plot(steps, demonEnergies, label='Demon energies')
        plt.plot(steps, demonEnergies + energies, label='Total energy')
        plt.legend()
    plt.ylabel("Energy")
    plt.xlabel("Timestep")
    plt.savefig("energy.png")
    plt.show()

    # SPINS
    plt.plot(steps, spins)
    plt.ylabel("Spin direction")
    plt.xlabel("Timestep")
    plt.savefig("spin.png")
    plt.show()