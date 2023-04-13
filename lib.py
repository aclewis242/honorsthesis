'''
A general-purpose collection of functions used in processing the simulation's outputs.
'''

import sys
import os
import numpy as np
import lattice as lat
import matplotlib.pyplot as plt
import random as rand
import pickle as pkl
from numpy import cos, pi

def longest(l: list) -> int:
    '''
    Finds the length of the longest element in a list (excluding trailing zeros).

    ### Parameters
    l: The list of iterables

    ### Returns
    maxLen: The longest length in l
    '''
    maxLen = -1
    for a in l:
        a = np.trim_zeros(a, 'b')
        if a.size > maxLen: maxLen = a.size
    return maxLen

def memAdd(o):
    '''
    Prints the memory address of a given object.

    ### Parameters
    o: The object
    '''
    print(hex(id(o)))

def genRandList(l: int, size: int, rev: bool=True) -> np.ndarray:
    '''
    Generates a list of random positions for the simulation to run through.

    ### Parameters
    l: The length of the desired list
    size: The size of the lattice
    rev: Whether or not the system is supposed to reverse

    ### Returns
    rL: The randomly-generated list
    '''
    rL = np.array([int(rand.random()*size) for i in range(l)])
    if rev: rL = np.append(rL, np.flip(rL))
    return rL

def procHex(*args):
    '''
    Turns a sequence of numbers (between 0 and 256) into a hexadecimal color code.
    Typically, only 3 inputs should be used (RGB), but a fourth is supported (RGBA).

    ### Parameters
    *args: The sequence of numbers
    '''
    return '#'+''.join([f'0{hex(s)[2:]}'[-2:] for s in args])

def getColor(i: float):
    '''
    Calculates RGB color based on position in a custom gradient and returns it as a hexadecimal color code.

    ### Parameters
    i: A decimal between 0 and 1 representing position within the gradient.
    '''
    r = int(-120*cos(pi*i) + 120)
    g = int(-120*cos(0.5*pi*i) + 120)
    b = int(70*cos(0.75*pi*(0.5*i + 0.25)) + 70)
    return procHex(r,g,b)

def dispGrad():
    '''
    Displays the gradient described in `getColor()` as an image.
    '''
    grad = np.linspace(0, 1, 100)
    for x in grad:
        plt.plot([x for z in np.zeros_like(grad)], grad, 'o', color=getColor(x+.0))
    plt.savefig('grad.png')
    plt.show()

def printDist(dist: np.ndarray, xax: str='t'):
    '''
    Prints a given distribution to the console.

    ### Parameters
    dist: The distribution to be printed
    xax: The variable on the x-axis of the distribution. Default is time, set to 'e' for energy
    '''
    name = 'Time block'
    if xax == 'e': name = 'Energy'
    [print(f'{name} {a}: {dist[a]}') for a in range(dist.shape[0])]

def norm(dist: np.ndarray):
    '''
    Normalises a given distribution.

    ### Parameters
    dist: The distribution to be normalised
    '''
    return dist/sum(dist)

def genLats(n: int=10, s: int=500, a: float=0.7):
    '''
    Generates several .lat files containing lattice objects.

    ### Parameters
    n: The number of lattices to generate
    s: The size of the lattices
    a: The alignment fraction to use (0-1)
    '''
    checkDir('lats')
    sys.setrecursionlimit(100000)
    for i in range(n):
        l = lat.Lattice(size=s, align=a)
        with open(f'lats/{l.energy}.lat', 'wb') as f: pkl.dump(l, f)

def load(l: lat.Lattice):
    '''
    'Load a lattice from a file,' sort of. It doesn't actually 'load' anything, per se, but it's useful for typing purposes (i.e.
    telling Python that the object loaded from the file is, in fact, a 'Lattice').

    ### Parameters
    l: The lattice object to load.
    '''
    return l

def checkDir(*args):
    '''
    Make sure the directory in question does, in fact, exist. Used in file saving.

    ### Parameters
    args: The directory(ies) to check. Works within the directory of the project.
    '''
    for dir in args:
        if not os.path.exists(dir): os.mkdir(dir)