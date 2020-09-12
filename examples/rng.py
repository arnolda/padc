from matplotlib.pyplot import *
from random import random

class R250:
    m = 2**31 - 1

    def __init__(self, init_rng):
        # Zustand mit anderem Generator initialisieren
        self.state = []
        for i in range(250):
            self.state.append(init_rng())

    def __call__(self):
        newval = (self.state[0] + self.state[250 - 103]) % self.m
        # neuen Wert anhaengen
        self.state.append(newval)
        # ... und aeltesten rauswerfen
        del self.state[0]
        return self.state[-1]

class Rand:
    a = 1103515245
    b = 12345
    m = 2**31

    def __init__(self, seed):
        self.state = seed

    def __call__(self):
        self.state = (self.state*self.a + self.b) % self.m
        return self.state

class Randu:
    m = 1 << 31 - 1
    a = 65539
    def __init__(self, seed):
        self.state = seed
    def __call__(self):
        self.state = (self.state*self.a) % (self.m + 1)
        return self.state - 1
    
class FloatRNG:
    def __init__(self, int_rng):
        self.rng = int_rng

    def __call__(self):
        return float(self.rng())/self.rng.m
        
for rng in FloatRNG(Randu(123)), FloatRNG(Rand(123)), FloatRNG(R250(Rand(123))), random:
    figure()

    series = [rng() for i in xrange(1000)]
    
    for start in xrange(len(series)):
        id = start
        row = [ series[id] ]
        finished = False
        for i in xrange(1,11):
            id = id + i
            if id >= len(series):
                finished = True
                break
            row.append(series[id])
        if not finished:
            plot(range(0,len(row)), row)

show()
