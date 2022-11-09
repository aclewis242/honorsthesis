import numpy as np
import lattice as lat
import matplotlib.pyplot as plt
import random as rand

if __name__ == "__main__":
    ### SPECIFY INPUT PARAMETERS HERE ###
    size = 100
    useDemon = True
    numSteps = 5000 # "Time"

    l = lat.Lattice(size=size)
    steps = range(numSteps)
    energies = np.zeros_like(steps, dtype=float)
    energies = np.append(energies, 0.0)
    demonEnergies = np.zeros_like(energies, dtype=float)
    spins = np.zeros_like(energies, dtype=float)
    dEdist = np.zeros_like(energies, dtype=int)
    energies[0] = l.energy()/l.size
    demonEnergies[0] = l.dE/l.size
    spins[0] = l.spin()
    name = "ising"
    algr = l.metropolis
    if useDemon:
        name = "demon"
        algr = l.metroDemon
    order = np.zeros_like(steps)
    for i in steps:
        if i < numSteps/2:
            order[i] = int(rand.random()*size)
        else:
            order[i] = order[steps[-1]-i]
    for i in range(order.size):
        l.lat = algr(ind=order[i])
        energies[i+1] = l.E/l.size
        demonEnergies[i+1] = l.dE/l.size
        spins[i+1] = l.spin()
        dEdist[int(l.dE/4)] += 1
    dEdist = np.trim_zeros(dEdist)/np.sum(dEdist)
    entr = -np.sum(dEdist*np.log10(dEdist))
    print(entr)
    print(l.entr())
    steps = range(numSteps+1)
    
    # ENERGIES
    plt.plot(steps, energies, label='System energy')
    if useDemon:
        plt.plot(steps, demonEnergies, label='Demon energy')
        plt.plot(steps, demonEnergies + energies, label='Total energy')
        plt.legend()
    plt.ylabel("Energy")
    plt.xlabel("Timestep")
    plt.savefig("energy-{}.png".format(name))
    plt.show()

    # ENERGY DISTRIBUTION
    if useDemon:
        plt.bar(range(dEdist.size), dEdist)
        plt.ylabel("P")
        plt.xlabel("Demon energy")
        plt.savefig("dedist.png")
        plt.show()

    # SPINS
    plt.plot(steps, spins)
    plt.ylabel("Spin direction")
    plt.xlabel("Timestep")
    plt.savefig("spin-{}.png".format(name))
    plt.show()