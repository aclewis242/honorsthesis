import numpy as np
import lattice as lat
import matplotlib.pyplot as plt
import random as rand

K = 1.380649e-23 # Boltzmann constant

if __name__ == "__main__":
    ### SPECIFY INPUT PARAMETERS HERE ###
    size = 100
    useDemon = True
    numSteps = 5000 # "Time"
    bkt = 1000 # Size of buckets for spin distribution

    l = lat.Lattice(size=size)
    steps = range(numSteps)
    energies = np.zeros_like(steps, dtype=float)
    energies = np.append(energies, 0.0)
    demonEnergies = np.zeros_like(energies, dtype=float)
    spins = np.zeros_like(energies, dtype=float)
    dEdist = np.zeros_like(energies, dtype=int)
    spdist = np.linspace(-bkt, bkt, 2*bkt+1)/bkt
    spdistEm = np.zeros_like(spdist, dtype=int)
    energies[0] = l.energy()/l.size
    demonEnergies[0] = l.dE/l.size
    spins[0] = l.spin()
    name = "ising"
    algr = l.metropolis
    if useDemon:
        name = "demon"
        algr = l.metroDemon2
    order = np.zeros_like(steps)
    for i in steps:
        if i < numSteps/2:
            order[i] = int(rand.random()*size)
        else:
            order[i] = order[steps[-1]-i]
    rev = False
    for i in range(order.size):
        if i >= numSteps/2: rev = True
        l.lat = algr(l.lat[order[i]], rev=rev)
        # print(l)
        energies[i+1] = l.E/l.size
        demonEnergies[i+1] = l.dE/l.size
        spins[i+1] = l.spin()
        dEdist[l.dE] += 1
        spdistEm[np.where(spdist == int(l.spin()*bkt)/bkt)[0]] += 1

    dEdist, dEentr = l.entr(dEdist)
    spdist, spentr = l.entr(spdistEm)
    print(dEentr)
    print(spentr)
    entr = dEentr + spentr
    entrTheo = l.entrDe() + l.entrSp()
    print(entr)
    print(entrTheo)
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

    # SPIN DISTRIBUTION
    plt.bar(range(spdist.size), spdist)
    plt.ylabel("P")
    plt.xlabel("Spin")
    plt.savefig("spdist.png")
    plt.show()