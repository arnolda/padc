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
# baryzentrische Polynom-Interpolation der Rungefunktion
############################################
import numpy
import scipy.interpolate as ip
import matplotlib.pyplot as pyplot

# die zu interpolierende Funktion
def runge(x): return 1/(1+x**2)

steps = 7
xs = numpy.linspace(-5,5,steps)
ys = runge(xs)

omegas = numpy.ones_like(xs)
for i in range(steps):
    for k in range(steps):
        if i != k:
            omegas[i] /= (xs[i] - xs[k])
def mu(i, x):
    return omegas[i]/(x - xs[i])

def bary(x):
    enum  = 0
    denom = 0
    for i in range(steps):
        enum  += ys[i]*mu(i, x)
        denom += mu(i, x)

    return enum/denom

# Ausgabe
#################################
px = numpy.logspace(-10, -5, 200)

figure = pyplot.figure(figsize=(4,4))

graph = figure.add_subplot(111)
graph.set_xscale("log")
graph.set_yscale("log")
graph.plot(px, abs(runge(px) - bary(px)), "r")
graph.plot(px, abs(runge(px) - ip.lagrange(xs, ys)(px)), "k")

pyplot.show()
