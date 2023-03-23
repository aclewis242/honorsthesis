from lib import *

if __name__ == "__main__":
    ### SPECIFY INPUT PARAMETERS HERE ###
    size = 500
    doDemons = True
    doFreeze = True
    doSum = True
    splitEvens = False
    timeExtension = 1
    blockCount = 27 # Keep this at a multiple of 3 so it syncs well with the bond making/breaking
    time = 50 # Change this and only this if you want to run it for more/less time
    # bkt = 1000 # Size of buckets for spin distribution

    l = lat.Lattice(size=size, align=0.65)
    numSteps = 6*blockCount*time
    steps = range(numSteps)
    # spins[0] = l.spin()
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
    # spins = np.zeros_like(energies, dtype=float)
    dEdist = np.zeros_like(energies, dtype=int)
    # spdist = np.linspace(-bkt, bkt, 2*bkt+1)/bkt
    # spdistEm = np.zeros_like(spdist, dtype=int)
    energies[0] = l.energy()/l.size
    demonEnergies[0] = l.dE/l.size
    dEblock = np.zeros_like(dEdist, dtype=int)
    for i in range(order.size):
        if i >= numSteps/2: rev = True
        if i >= numSteps/2 - numSteps/6 and i < numSteps/2 + numSteps/6: brk = not doFreeze
        else: brk = True
        l.lat = algr(l.lat[order[i]], brk=brk, rev=rev)
        energies[i+1] = l.E/l.size
        demonEnergies[i+1] = l.dE/l.size
        # spins[i+1] = l.spin()
        if doDemons: dEblock[int(l.dE)] += 1
        else: dEblock[abs(int(l.E/4))] += 1
        if i%eDiffBlock == 0 and i != 0:
            eDiffsUnp.append(dEblock.copy())
            dEblock_trim = np.trim_zeros(dEblock, 'b')
            plt.bar(range(dEblock_trim.size), dEblock_trim)
            plt.savefig(f'blocks-{name}/{int(i/eDiffBlock)}.png')
            plt.close()
            if not doSum: dEblock = np.zeros_like(dEdist, dtype=int)
        # spdistEm[np.where(spdist == int(l.spin()*bkt)/bkt)[0]] += 1
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
    print(f'eDiffsLn shape: {eDiffsLn.shape}')
    eDiffs = []
    for e in eDiffsLn:
        diff = np.zeros(e.size)
        for i in np.arange(1, e.size):
            diff[i] = e[i] - e[i-1-splitEvens]
            if abs(diff[i]) == np.inf or i-1-splitEvens < 0: diff[i] = None
        # print(f'diffs size: {diff.size}')
        eDiffs.append(diff)
    eDiffs = np.array(eDiffs)
    print(f'eDiffs shape: {eDiffs.shape}')
    
    # ENERGIES
    plt.plot(energies, label='System energy')
    if doDemons:
        plt.plot(demonEnergies, label='Demon energy')
        plt.plot(demonEnergies + energies, label='Total energy')
        plt.legend()
    plt.ylabel('Energy')
    plt.xlabel('Timestep')
    plt.savefig(f'energy-{name}.png')
    plt.show()

    # ENERGY DISTRIBUTION
    plt.bar(range(dEdist.size), dEdist)
    plt.ylabel('P')
    plt.xlabel('Energy')
    plt.savefig(f'edist-{name}.png')
    plt.show()

    # ENERGY LN
    # if doDemons:
    #     lndE = np.log(demonEnergies) - np.log(demonEnergies[0])
    #     plt.plot(lndE)
    #     plt.show()
    
    # ENERGY DIFFS
    tblocks = np.arange(1, (2*timeExtension + (not timeExtension))*blockCount+1)*eDiffBlock
    eDiffsT = np.transpose(eDiffs)
    [plt.plot(tblocks, eDiffsT[e], color=getColor(e/len(eDiffsT))) for e in range(len(eDiffsT))]
    plt.title('Difference of energy natural logarithms')
    plt.savefig(f'ediffs-{name}.png')
    plt.show()

    # # SPINS
    # plt.plot(steps, spins)
    # plt.ylabel("Spin direction")
    # plt.xlabel("Timestep")
    # plt.savefig("spin-{}.png".format(name))
    # plt.show()

    # # SPIN DISTRIBUTION
    # plt.bar(range(spdist.size), spdist)
    # plt.ylabel("P")
    # plt.xlabel("Spin")
    # plt.savefig("spdist.png")
    # plt.show()