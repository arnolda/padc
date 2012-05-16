# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
from scipy import *
import scipy.optimize as opt
import numpy.random as rand
import matplotlib.pyplot as pyplot

N=200

def err(c, x, y):
    return y - c[0]*sin(c[1]*x+c[2])

x = linspace(0, 2*pi, N+1)
rand.seed(1)
y = sin(x) + 0.1*sin(10*x) + rand.normal(0, 0.1, x.shape)
c, success = opt.leastsq(err, (0, 1, 0), args=(x[0:N], y[0:N]))

print c

# Ausgabe
#################################

figure = pyplot.figure(figsize=(4,4))

graph = figure.add_subplot(111)
graph.plot(x, y, "ko", markersize=1)
graph.plot(x, c[0]*sin(c[1]*x+c[2]), "r--",linewidth=1.5)
graph.axis([0, 2*pi,-1.1,1.1])

figure.savefig("leastsq.pdf")
