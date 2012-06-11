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
