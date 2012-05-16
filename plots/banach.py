# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Sukzessive Substitution
#
############################################
from numpy import *
import math
import matplotlib.pyplot as pyplot

def g(r, phi0):
    return exp(-r)/phi0

def graphical_subst(x0, g, n):
    "graphische sukzessive Substitution durch Spiegelung an der Winkelhalbierenden"
    
    # Schnecke
    ssx = []
    ssy = []
    # Achsenmarker
    achse = []

    xn = x0

    for r in range(0, n):
        ssx.append(xn)
        achse.append((xn, g(xn)))
        xn = g(xn)
        ssy.append(xn)
        ssx.append(xn)
        ssy.append(xn)

    achse.append((xn, xn))

    return ssx, ssy, achse

#######################################################

phi0 = 2
n = 7
x0 = 0.05

rmax=0.5

ssx, ssy, achse = graphical_subst(x0, lambda r: g(r, phi0), n)
        
figure = pyplot.figure(figsize=(8,4))

graph = figure.add_subplot(121)

x = linspace(0, rmax, 200)

graph.plot(x, g(x, phi0), "b-")
graph.plot((0,rmax), (0,rmax), "k-")

graph.plot(ssx, ssy, "r-")
graph.plot(ssx[1::2], ssy[1::2], "ro", markersize=3)

for n in range(len(achse)):
    xn, gxn = achse[n]
    graph.plot((xn, xn), (0, gxn), "k:", linewidth=1)
    if n < 3:
        graph.text(xn+.01, 0, "$x_%d$" % n)

graph.axis([0,rmax,0,rmax])

#######################################################

phi0 = 0.25
x0 = 1
n=10
rmin=0
rmax=4

ssx, ssy, achse = graphical_subst(x0, lambda r: g(r, phi0), n)

graph = figure.add_subplot(122)

x = linspace(rmin, rmax, 200)

graph.plot(x, g(x, phi0), "b-")
graph.plot((rmin,rmax), (rmin,rmax), "k-")

graph.plot(ssx, ssy, "r-")
graph.plot(ssx[1::2], ssy[1::2], "ro", markersize=3)

for n in range(len(achse)):
    xn, gxn = achse[n]
    graph.plot((xn, xn), (rmin, gxn), "k:", linewidth=1)
    if n == 0 or n > 8:
        graph.text(xn+.01, 0, "$x_{%d}$" % n)

graph.axis([rmin,rmax,rmin,rmax])

figure.savefig("banach.pdf")
