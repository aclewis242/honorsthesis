import numpy as np
import lattice as lat
import matplotlib.pyplot as plt
import random as rand
from numpy import cos, pi

K = 1.380649e-23 # Boltzmann constant

def longest(l: list):
    maxLen = -1
    for a in l:
        a = np.trim_zeros(a, 'b')
        if a.size > maxLen: maxLen = a.size
    return maxLen

def memAdd(o):
    print(hex(id(o)))

def genRandList(l: int, size: int, rev: bool=True):
    rL = []
    [rL.append(int(rand.random()*size)) for i in range(l)]
    rL = np.array(rL)
    if rev: rL = np.append(rL, np.flip(rL))
    return rL

def procHex(*args): return '#'+''.join([f'0{hex(s)[2:]}'[-2:] for s in args])

def getColor(i):
    r = int(-120*cos(pi*i) + 120)
    g = int(-120*cos(0.5*pi*i) + 120)
    b = int(70*cos(0.75*pi*(0.5*i + 0.25)) + 70)
    return procHex(r,g,b)

def dispGrad():
    grad = np.linspace(0, 1, 100)
    for x in grad:
        plt.plot([x for z in np.zeros_like(grad)], grad, 'o', color=getColor(x+.0))
    plt.savefig('grad.png')
    plt.show()

def printDist(dist: np.ndarray, xax: str='t'):
    name = 'Time block'
    if xax == 'e': name = 'Energy'
    [print(f'{name} {a}: {dist[a]}') for a in range(dist.shape[0])]