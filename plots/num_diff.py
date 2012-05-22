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
# Numerische Differentiation
#
############################################
from numpy import *
import math
import matplotlib.pyplot as pyplot

a=100

def f(x):
    return sin(a*x**2)

def fprime(x):
    return a*cos(a*x**2)*2*x

def left(f, x, h):
    return (f(x) - f(x-h))/h

def central(f, x, h):
    return (f(x+h) - f(x-h))/h/2.0

def order5(f, x, h):
    return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))/h/12.0

# Daten generieren
# Testpunkte, an der die Ableitung angenaehert wird
xtest = arange(-pi, pi, 10000)

hs = logspace(-12,0,20)

left_prec    = []
central_prec = []
order5_prec  = []
for h in hs:
    left_prec.append(max( abs(left(f, x, h) - fprime(x)) for x in xtest))
    central_prec.append(max( abs(central(f, x, h) - fprime(x)) for x in xtest))
    order5_prec.append(max( abs(order5(f, x, h) - fprime(x)) for x in xtest))

figure = pyplot.figure(figsize=(4,4))
figure.subplots_adjust(bottom=0.14)

graph = figure.add_subplot(111)
graph.set_xlabel("h")
graph.set_xscale("log")
graph.set_yscale("log")
graph.set_xticks(logspace(-12, 0, 7))

graph.plot(hs, left_prec, "r-")
graph.plot(hs, central_prec, "g:")
graph.plot(hs, order5_prec, "b--")
# Zum Verifizieren der Ordnungen, mit geschaetzten Ableitungen
# graph.plot(hs, 100*(100*2*pi)**2*hs**2)
# graph.plot(hs, 100*(100*2*pi)**4*hs**4)

figure.savefig("num_diff.pdf")
