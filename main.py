from lib import *

if __name__ == "__main__":
    ### SPECIFY INPUT PARAMETERS HERE ###
    size = 500          # Sets lattice size
    doDemons = True     # Decides whether to use the Creutz or the Metropolis algorithm. (Creutz is 'True')
    doFreeze = True     # Turn bond freezing on/off for the middle third
    doSum = True        # Decides whether to compute the log. differences of the total or instantaneous distributions.
    splitEvens = False  # Decides whether to compute the log. differences normally or odd-to-odd, even-to-even only.
    timeExtension = 1   # How many extra times you want the simulation to run. Only meaningful for Creutz with freezing on
    blockCount = 27     # Keep this at a multiple of 3 so it syncs well with the bond making/breaking
    time = 50           # Change this and only this if you want to run it for more/less time

    ### INITIAL SETUP ###
    l = lat.Lattice(size=size, align=0.65)
    numSteps = 6*blockCount*time
    steps = range(numSteps)
    name = 'ising'
    algr = l.metropolis
    if doDemons:
        name = 'demon'
        algr = l.metroDemon
    else:
        timeExtension = 0
        splitEvens = False
    if not doFreeze: timeExtension = 0
    order = genRandList(int(numSteps/2), size)
    rev = False
    eDiffBlock = 6*time
    eDiffsUnp = []
    if timeExtension: order = np.array([np.append(order, genRandList(numSteps, size, rev=False)) for i in range(timeExtension)]).flatten()
    energies = np.zeros_like(order, dtype=float)
    energies = np.append(energies, 0.0)
    demonEnergies = np.zeros_like(energies, dtype=float)
    dEdist = np.zeros_like(energies, dtype=int)
    energies[0] = l.energy
    demonEnergies[0] = l.dE
    dEblock = np.zeros_like(dEdist, dtype=int)

    ### SIMULATION LOOP ###
    for i in range(order.size):
        if i >= numSteps/2: rev = True
        if i >= numSteps/2 - numSteps/6 and i < numSteps/2 + numSteps/6: brk = not doFreeze
        else: brk = True
        algr(l.lat[order[i]], brk=brk, rev=rev)
        energies[i+1] = l.E
        demonEnergies[i+1] = l.dE
        if doDemons: dEblock[int(l.dE)] += 1
        else: dEblock[abs(int(l.E/4))] += 1
        if i%eDiffBlock == 0 and i != 0:
            eDiffsUnp.append(dEblock.copy())
            dEblock_trim = np.trim_zeros(dEblock, 'b')
            ind = int(i/eDiffBlock)
            plt.bar(range(dEblock_trim.size), norm(dEblock_trim))
            plt.xlabel('Energy')
            plt.ylabel('Probability')
            plt.title(f'Distribution at time block {ind}')
            plt.savefig(f'blocks-{name}/{ind}.png')
            plt.close()
            if not doSum: dEblock = np.zeros_like(dEdist, dtype=int)
    
    ### ENERGY LOGARITHM DIFFERENCE PROCESSING ###
    eDiffsUnp.append(dEblock)
    eDiffsUnp = np.array(eDiffsUnp)
    if doSum: dEdist = dEblock
    steps = range(numSteps+1)
    eDiffsNZ = []
    lgt = longest(eDiffsUnp)
    for e in eDiffsUnp:
        dEdist += e
        e[lgt] = '7'
        eDiffsNZ.append(np.trim_zeros(e, 'b')[:lgt])
    eDiffsNZ = np.array(eDiffsNZ)
    dEdist = np.trim_zeros(dEdist)
    eDiffsLn = np.log(eDiffsNZ)
    eDiffs = []
    for e in eDiffsLn:
        diff = np.zeros(e.size)
        for i in np.arange(1, e.size):
            diff[i] = e[i] - e[i-1-splitEvens]
            if abs(diff[i]) == np.inf or i-1-splitEvens < 0: diff[i] = None
        eDiffs.append(diff)
    eDiffs = np.array(eDiffs)
    
    ### RESULTS PLOTTING
    # Energy over time
    plt.plot(energies, label='System energy')
    if doDemons:
        plt.plot(demonEnergies, label='Demon energy')
        plt.plot(demonEnergies + energies, label='Total energy')
        plt.legend()
    plt.ylabel('Energy')
    plt.xlabel('Timestep')
    plt.title('Energy in time')
    plt.savefig(f'energy-{name}.png')
    plt.show()

    # Final energy distribution
    plt.bar(range(dEdist.size), norm(dEdist))
    plt.ylabel('Probability')
    plt.xlabel('Energy')
    plt.title('Final energy distribution')
    plt.savefig(f'edist-{name}.png')
    plt.show()
    
    # Energy differences
    tblocks = np.arange(1, (2*timeExtension + (not timeExtension))*blockCount+1)*eDiffBlock
    eDiffsT = np.transpose(eDiffs)
    [plt.plot(tblocks, eDiffsT[e], color=getColor(e/len(eDiffsT))) for e in range(len(eDiffsT))]
    plt.xlabel('Timestep')
    plt.ylabel(r'$\Delta\ln(E)$')
    plt.title('Difference of energy natural logarithms')
    plt.savefig(f'ediffs-{name}.png')
    plt.show()