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
# Bisektion
#
############################################
from numpy import *
import math
import matplotlib.pyplot as pyplot

def f(r, phi0):
    return exp(-r)/r - phi0

def graphical_bisection(a0, b0, f, n):
    "Bisektion"

    an = a0
    bn = b0

    res = []

    for r in range(0, n):
        mn = (an + bn)/2.0
        res.append((an, bn, mn))

        if f(mn)*f(an) < 0:
            bn = mn
        else:
            an = mn

    return res

#######################################################

phi0 = 2
n = 7

a0 = 0.1
b0 = 1

rmin=0.01
rmax=1.01

path = graphical_bisection(a0, b0, lambda r: f(r, phi0), n)
        
figure = pyplot.figure(figsize=(4,4))

graph = figure.add_subplot(111)

x = linspace(rmin, rmax, 200)

graph.plot(x, f(x, phi0), "b-")
graph.plot(x, zeros(x.shape), "b:",linewidth=0.5)

for an, bn, mn in reversed(path):
    h = f(mn, phi0)
    graph.plot((mn, mn), (0, h), "k--", dashes=(1,1), linewidth=1)
    graph.plot((an, bn), (h, h), "r-", linewidth=1)
    graph.plot((an, an), (h-0.01, h+0.01), "r-", linewidth=1)
    graph.plot((bn, bn), (h-0.01, h+0.01), "r-", linewidth=1)

print ("%.4f %.4f %.4f" % path[-1])

graph.axis([0,rmax,-1.5,1])

figure.savefig("bisektion.pdf")
