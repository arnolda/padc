# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Taylorreihen des Sinus
#
############################################
import numpy
import math
import matplotlib.pyplot as pyplot

x = numpy.linspace(-numpy.pi, numpy.pi, 200)

figure = pyplot.figure(figsize=(4,4))

def symm(x):
    if x < -numpy.pi/2:
        return -numpy.pi-x
    elif x > numpy.pi/2:
        return numpy.pi-x
    else:
        return x

px = numpy.array([symm(v) for v in x ])

graph = figure.add_subplot(111)
graph.axis([-numpy.pi,numpy.pi,-1.5,1.5])
graph.plot(x, numpy.sin(px), "black", linewidth=2)
graph.plot(x, px, "b:", linewidth=1)
graph.plot(x, px - px**3/6, "g-", linewidth=1)
graph.plot(x, px  - px**3/6 + px**5/math.factorial(5), "r--", linewidth=1)

figure.savefig("sinus.pdf")
