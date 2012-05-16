# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Fourierreihen
#
############################################
from numpy import *
import math
import matplotlib.pyplot as pyplot

x = linspace(0, 2*pi, 200)

figure = pyplot.figure(figsize=(8,4))

def f(x):
    return 2*(x < pi) - 1

def series(x,n):
    s = 0
    for k in range(1,n+1):
        s += sin((2*k-1)*x)/(2*k-1)
    return 4/pi*s

graph = figure.add_subplot(121)
graph.plot(x, f(x), "black", linewidth=2)
graph.plot(x, series(x,1), "b--", linewidth=1)
graph.plot(x, series(x,2), "r:", linewidth=1.5)
graph.plot(x, series(x,20), "g-", linewidth=1)
graph.axis([0,2*pi,-1.5,1.5])

def f(x):
    return (x < pi)*x + (x > pi)*(2*pi - x)

def series(x,n):
    s = pi/2
    for k in range(1,n+1):
        s += -4/pi/(2*k-1)**2*cos((2*k-1)*x)
    return s

graph = figure.add_subplot(122)
graph.plot(x, f(x), "black", linewidth=1)
graph.plot(x, series(x,1), "b--", linewidth=2)
graph.plot(x, series(x,2), "r:", linewidth=1.5)
graph.plot(x, series(x,20), "g-", linewidth=1)
graph.axis([0,2*pi,-0.5,3.5])

figure.savefig("fourier.pdf")
