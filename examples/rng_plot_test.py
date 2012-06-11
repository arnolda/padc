# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
#
# Pseudozufallszahlentests
#
############################################
from scipy import *
from numpy.random import *
import math
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D

class Randmine:
    def __init__(self):
        self.state = 123
    def next(self):
        m = 1 << 32
        a = 49
        b = 1975
        self.state = (a * self.state + b) % m
        return float(self.state & ((1 << 31) - 1))/((1 << 31) - 1)

class Rand:
    def __init__(self):
        self.state = 123
    def next(self):
        m = 1 << 32
        a = 1103515245
        b = 12345
        self.state = (a * self.state + b) % m
        return float(self.state & ((1 << 31) - 1))/((1 << 31) - 1)

class Randu:
    def __init__(self):
        self.state = 123
    def next(self):
        m = 1 << 31
        a = 65539
        self.state = (self.state*a) % m
        return float(self.state - 1)/float(m - 1)

class Minstd:
    def __init__(self):
        self.state = 123
    def next(self):
        m = (1 << 31) - 1
        a = 16807
        self.state = (self.state*a) % m
        return float(self.state - 1)/float(m - 1)

##########################################

rng_type = Randu

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.15,wspace=0.3, left=0.1,right=0.95)

##########################################

rng = rng_type()

data = [ rng.next() for x in range(10000) ]

graph = figure.add_subplot(121)

graph.plot(data[:-1], data[1:], ".", markersize=1)

##########################################

rng = rng_type()
data = [ rng.next() for x in range(10000) ]

rect = figure.add_subplot(1, 2, 2).get_position()
graph = Axes3D(figure, rect)

graph.scatter(data[:-2], data[1:-1], data[2:], s=1, marker="o",
              edgecolors="blue", facecolors="blue")

pyplot.show()
