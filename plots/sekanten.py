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
# Regula falsi
#
############################################
from numpy import *
import math
import matplotlib.pyplot as pyplot

def f(r, phi0):
    return exp(-r)/r - phi0

def graphical_regulafalsi(x0, x1, f, n):
    "Regula falsi"

    xnn = x0
    xn  = x1

    res = [x0, x1]

    for r in range(0, n):
        fprime = (f(xn) - f(xnn))/(xn - xnn)
        xnn = xn
        xn = xn - f(xn)/fprime
        res.append(xn)

    return res

#######################################################

phi0 = 2
n = 7
x0 = 0.05
x1 = 0.1

rmin=0.01
rmax=0.5

path = graphical_regulafalsi(x0, x1, lambda r: f(r, phi0), n)
        
figure = pyplot.figure(figsize=(8,4))

graph = figure.add_subplot(121)

x = linspace(rmin, rmax, 200)

graph.plot(x, f(x, phi0), "b-")
graph.plot(x, zeros(x.shape), "b:",linewidth=0.5)

for n in range(2, len(path)):
    xnnn = path[n-2]
    xnn = path[n-1]
    xn = path[n]

    graph.plot((xnnn, xnnn), (0, f(xnnn, phi0)), "k:", linewidth=1)
    graph.plot((xnnn, xnn), (f(xnnn, phi0), f(xnn, phi0)), "r--", linewidth=2)
    graph.plot((xnnn, xnn), (f(xnnn, phi0), f(xnn, phi0)), "ko", markersize=4)
    graph.plot((xnn, xn), (f(xnn, phi0), 0), "r-", linewidth=2)
    if n < 7:
        graph.text(xnn, -1.3, "$x_%d$" % (n-1),
                   horizontalalignment='center')

print abs(xn - xnn)

graph.axis([0,rmax,-1.5,20])

#######################################################

phi0 = 0.25
n = 7
x0 = 1.8
x1 = 1.9

rmin=0.5
rmax=2

path = graphical_regulafalsi(x0, x1, lambda r: f(r, phi0), n)
        
graph = figure.add_subplot(122)

x = linspace(0.01, rmax, 200)

graph.plot(x, f(x, phi0), "b-")
graph.plot(x, zeros(x.shape), "b:",linewidth=0.5)

for n in range(2, len(path)):
    xnnn = path[n-2]
    xnn = path[n-1]
    xn = path[n]

    graph.plot((xnnn, xnnn), (0, f(xnnn, phi0)), "k:", linewidth=1)
    graph.plot((xnnn, xnn), (f(xnnn, phi0), f(xnn, phi0)), "r--", linewidth=2)
    graph.plot((xnnn, xnn), (f(xnnn, phi0), f(xnn, phi0)), "ko", markersize=4)
    graph.plot((xnn, xn), (f(xnn, phi0), 0), "r-", linewidth=2)
    if n < 6:
        graph.text(xnn, -0.06, "$x_%d$" % (n-1),
                   horizontalalignment='center')

graph.axis([rmin, rmax,-0.25,0.75])
graph.set_xticks(arange(rmin, rmax+0.01, 0.25))

figure.savefig("regulafalsi.pdf")
