from lib import *

def run(temp=2.5):      # 'temp' is used only in the Metropolis algorithm
    ### SPECIFY INPUT PARAMETERS HERE ###
    size = 500          # Sets lattice size
    doDemons = True     # Decides whether to use the Creutz or the Metropolis algorithm (Creutz is 'True')
    doFreeze = True     # Turn bond freezing on/off for the middle third
    doBrks = True       # Decides whether or not to allow bond formation/breaking
    doSum = True        # Decides whether to compute the log. differences of the total or instantaneous distributions
    splitEvens = False  # Decides whether to compute the log. differences normally or odd-to-odd, even-to-even only
    timeExtension = 1   # How many extra times you want the simulation to run. Only meaningful for Creutz with freezing on
    blockCount = 81     # Keep this at a multiple of 3 so it syncs well with the bond making/breaking
    time = 20           # Change this and only this if you want to run it for more/less time
    eLim = 50           # Maximum energy to show on plots

    ### INITIAL SETUP ###
    genLats(s=size)
    l = lat.Lattice(size=size, align=0.65)
    # with open('lattice.lat', 'rb') as f: l = pkl.load(f)
    l = load(l)
    if not doDemons: l.dE = -size - l.energy
    initL = l.__repr__()
    numSteps = 6*blockCount*time
    name = 'metro'
    algr = l.metropolis
    if doDemons:
        name = 'creutz'
        algr = l.demon
    else:
        timeExtension = 0
        splitEvens = False
        doBrks = False
    if not doFreeze: timeExtension = 0
    if not doBrks: doFreeze = True
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
    sys.setrecursionlimit(10000)
    factor = 4
    if doBrks: factor = 1
    with open(f'init-{name}.lat', 'wb') as f: pkl.dump(l, f)

    ### SIMULATION LOOP ###
    for i in range(order.size):
        if i >= numSteps/2: rev = True
        if i >= numSteps/2 - numSteps/6 and i < numSteps/2 + numSteps/6: brk = not doFreeze
        else: brk = doBrks
        algr(l.lat[order[i]], brk=brk, rev=rev, temp=temp)
        energies[i+1] = l.E
        demonEnergies[i+1] = l.dE
        dEblock[abs(int(l.dE/factor))] += 1
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
    with open(f'fin-{name}.lat', 'wb') as f: pkl.dump(l, f)
    if doDemons and l.__repr__() == initL: print('Reversibility condition upheld!')
    elif doDemons and not timeExtension: print('Reversibility condition failed!')
    
    ### ENERGY LOGARITHM DIFFERENCE PROCESSING ###
    eDiffsUnp.append(dEblock)
    eDiffsUnp = np.array(eDiffsUnp)
    if doSum: dEdist = dEblock
    eDiffsNZ = []
    lgt = longest(eDiffsUnp)
    if lgt > eLim: lgt = eLim
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
    plt.plot(demonEnergies, label='Demon energy')
    plt.plot(demonEnergies + energies, label='Total energy')
    plt.legend()
    plt.ylabel('Energy')
    plt.xlabel('Timestep')
    plt.title('Energy in time')
    plt.savefig(f'energy-{name}.png')
    plt.close()

    # Final energy distribution
    dEdist = dEdist[:eLim]
    colors = [getColor(e/dEdist.size) for e in range(dEdist.size)]
    plt.bar(range(dEdist.size), norm(dEdist), color=colors)
    plt.ylabel('Probability')
    plt.xlabel(f'Energy ({factor}x)')
    plt.title('Final energy distribution')
    plt.savefig(f'edist-{name}.png')
    plt.close()
    
    # Energy differences
    tblocks = np.arange(1, (2*timeExtension + (not timeExtension))*blockCount+1)*eDiffBlock
    eDiffsT = np.transpose(eDiffs)
    [plt.plot(tblocks, eDiffsT[e], color=getColor(e/len(eDiffsT))) for e in range(len(eDiffsT))]
    plt.xlabel('Timestep')
    plt.ylabel(r'$\Delta\ln(E)$')
    plt.title('Difference of energy natural logarithms')
    plt.savefig(f'ediffs-{name}.png')
    plt.close()

    # Energy differences (filtered)
    tblocks = np.arange(1, (2*timeExtension + (not timeExtension))*blockCount+1)*eDiffBlock
    eDiffsT = np.transpose(eDiffs)
    t1 = 0
    for e in range(int(len(eDiffsT)/2)):
        plt.plot(tblocks, eDiffsT[e], color=getColor(2*e/len(eDiffsT)), label=f'E={factor*e}')
        if not timeExtension:
            temp = -factor/np.mean((eDiffsT[e])[int(len(eDiffsT[e])/2):])
            if e == 1: t1 = temp
            print(f'Temperature (E={factor*e}):\t{temp}')
    plt.xlabel('Timestep')
    plt.ylabel(r'$\Delta\ln(E)$')
    plt.title('Difference of energy natural logarithms (filtered)')
    plt.savefig(f'ediffs-{name}-filtered-noleg.png')
    plt.legend()
    plt.savefig(f'ediffs-{name}-filtered.png')
    plt.close()

    return t1

def tempDist():
    temps = np.array(np.linspace(3.5, 8, 20))
    diffs = []
    for t in temps:
        dEdiff = run(t)
        diffs.append(dEdiff)
    # temps.tofile('temps.csv')
    # np.array(diffs).tofile('diffs.csv')
    plt.plot(temps, diffs, linestyle='None', marker='o')
    plt.title('Energy divergence')
    plt.ylabel(r'$\Delta E$')
    plt.xlabel('Temperature')
    # plt.savefig('ediv.png')
    # plt.show()

if __name__ == '__main__':
    run()