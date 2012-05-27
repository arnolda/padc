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

def graphical_regulafalsi(a0, b0, f, n):
    "Regula falsi"

    an = a0
    bn = b0

    res = []

    for r in range(0, n):
        mn = (f(an)*bn - f(bn)*an)/(f(an) - f(bn))
        res.append((an, bn, mn))

        if f(mn)*f(an) < 0:
            bn = mn
        else:
            an = mn

    return res

#######################################################

figure = pyplot.figure(figsize=(8,4))

rmin = 0.01
rmax = 1.01
step = 0.05

phi0 = 2
n = 7

a0 = 0.25
b0 = 1

#######################################################

graph = figure.add_subplot(121)

path = graphical_bisection(a0, b0, lambda r: f(r, phi0), n)
        
x = linspace(rmin, rmax, 200)

graph.plot(x, f(x, phi0), "b-")
graph.plot(x, zeros(x.shape), "b:",linewidth=0.5)

h = -step
for an, bn, mn in reversed(path):
    fx = f(mn, phi0)
    graph.plot((mn, mn), (0, fx), "k--", dashes=(1,1), linewidth=1)
    graph.plot((an, bn), (h, h), "g-", linewidth=1)
    graph.plot((an, an), (h-0.01, h+0.01), "g-", linewidth=1)
    graph.plot((bn, bn), (h-0.01, h+0.01), "g-", linewidth=1)
    h -= step

print ("%.4f %.4f %.4f" % path[-1])

graph.axis([0,rmax,-1.5,1])

#######################################################

graph = figure.add_subplot(122)

path = graphical_regulafalsi(a0, b0, lambda r: f(r, phi0), n)
        
x = linspace(rmin, rmax, 200)

graph.plot(x, f(x, phi0), "b-")
graph.plot(x, zeros(x.shape), "b:",linewidth=0.5)

h = -step
for an, bn, mn in reversed(path):
    fx = f(mn, phi0)
    graph.plot((mn, mn), (0, fx), "k--", dashes=(1,1), linewidth=1)
    graph.plot((an, bn), (h, h),  "g-", linewidth=1)
    graph.plot((an, an), (h-0.01, h+0.01), "g-", linewidth=1)
    graph.plot((bn, bn), (h-0.01, h+0.01), "g-", linewidth=1)
    graph.plot((an, bn), (f(an, phi0), f(bn, phi0)), "r--", linewidth=1)
    graph.plot((an, bn), (f(an, phi0), f(bn, phi0)), "ko", markersize=4)
    h -= step

print ("%.4f %.4f %.4f" % path[-1])

graph.axis([0,rmax,-1.5,1])

figure.savefig("bisektion.pdf")
