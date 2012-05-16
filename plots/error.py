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

# Fehlerabschätzung korrelierter Daten
############################################

from numpy import *
import math
import numpy.random as rand
import matplotlib.pyplot as pyplot

# 150 unabhängige Beispieldaten aus den Geschwindigkeitskomponenten
# der Teilchen
vac = loadtxt("v_0.5.data.gz").transpose()[1:]

N = 200
deltas = arange(1,30,1)

vars    = zeros(deltas.shape)
varvars = zeros(deltas.shape)

for k in range(len(deltas)):
    delta = deltas[k]
    sumvar  = 0
    sum2var = 0
    S = 0
    for sample in vac:
        for offset in range(0,delta,10):
            data = sample[offset:offset + N*delta:delta]
            if len(data) != N:
                print "das geht schief!!", len(data), N
            mean = sum(data)/N
            var = (sum(data*data) - N*mean**2)/(N-1)
            # Varianz vom Datensatz
            sumvar  += var
            sum2var += var*var
            S += 1
    # Mittelwert und Varianz der Varianzen
    meanvar = sumvar/S
    varvar = (sum2var - S*meanvar**2)/(S - 1)
    vars[k]    = meanvar
    varvars[k] = varvar

figure = pyplot.figure(figsize=(4,4))
figure.subplots_adjust(left=0.2)

t = 0.01*deltas
graph = figure.add_subplot(111)
graph.set_xlabel(u"Δ")
graph.set_ylabel(u"σ²-Schätzer")
graph.errorbar(t, vars, yerr=varvars,  linestyle="o")
t = linspace(0.005,0.3,100)
graph.plot(t, 1 - 2*0.17/t/N , "k-")
graph.axis((0,0.305,0.5,1.2))
figure.savefig("error.pdf")
