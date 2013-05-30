# -*- coding: utf-8 -*-

# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Binning-Analyse nach Jahnke
############################################

from numpy import *
import math
import numpy.random as rand
import matplotlib.pyplot as pyplot

# Beispieldaten aus der Geschwindigkeitskomponente
# eines Teilchens
vac = loadtxt("v_0.5.data.gz")[:,42]
N=len(vac)
kmax = 1000
Delta = 0.01
# Variance of the data
sigma=1.0

vars = []
errs = []
ks = arange(1, kmax)

for k in ks:
    Nb = N/k
    Ntot = Nb*k
    data = vac[:Ntot]

    print k, Nb, data.shape

    blocked = array([ mean(data[block*k:(block+1)*k]) for block in range(Nb) ])
    v = var(blocked, ddof=1)

    vars.append(v)
    errs.append(v/Nb)

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.15,bottom=0.15,wspace=0.3,right=0.95)

# left, error estimate
###################################
graph = figure.add_subplot(121)
graph.set_xscale("log")
graph.set_yscale("log")
graph.set_xlabel("$k$")
graph.set_ylabel(u"$\epsilon^2(k)$")
graph.plot(ks, errs , "b.", markeredgecolor="b")
head = array([ k for k in ks if k < 20 ])
graph.plot(head, 0.25e-4*head , "k--")

# left, error estimate
###################################
graph = figure.add_subplot(122)
graph.set_xscale("log")
graph.set_yscale("log")
graph.set_xlabel("$k$")
graph.set_ylabel(u"$\\tau_{int}$")

print "expected variance", sigma, "measured", var(vac, ddof=1)
graph.plot(ks, Delta*array(vars)*ks/2.0/sigma , "r.", markeredgecolor="r")

figure.savefig("binning.pdf")
