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
# Newtonverfahren
#
############################################
import math
import matplotlib.pyplot as plt
import numpy as np

def f(r, phi0):
    return np.exp(-r)/r - phi0

def fprime(r, phi0):
    return np.exp(-r)*(-r - 1)/r**2

def graphical_newton(x0, f, fprime, n):
    "Newtonverfahren"

    xn = x0

    res = [x0]

    for r in range(0, n):
        xn = xn - f(xn)/fprime(xn)
        res.append(xn)

    return res

#######################################################

phi0 = 2
n = 7
x0 = 0.05

rmin=0.01
rmax=0.5

path = graphical_newton(x0, lambda r: f(r, phi0),
                        lambda r: fprime(r, phi0), n)
        
figure = plt.figure(figsize=(8,4))

graph = figure.add_subplot(121)

x = np.linspace(rmin, rmax, 200)

graph.plot(x, f(x, phi0), "b-")
graph.plot(x, np.zeros(x.shape), "b:",linewidth=0.5)

for n in range(1, len(path)):
    xnn = path[n-1]
    xn = path[n]

    graph.plot((xnn, xnn), (0, f(xnn, phi0)), "k:", linewidth=1)
    graph.plot((xnn, xn), (f(xnn, phi0), 0), "r-", linewidth=2)
    if n < 7:
        graph.text(xnn, -1.3, "$x_%d$" % (n-1),
                   horizontalalignment='center')

print(abs(xn - xnn))

graph.axis([0,rmax,-1.5,20])

#######################################################

phi0 = 0.25
n = 7
x0 = 1.8

rmin=0.5
rmax=2

path = graphical_newton(x0, lambda r: f(r, phi0),
                        lambda r: fprime(r, phi0), n)
        
graph = figure.add_subplot(122)

x = np.linspace(0.01, rmax, 200)

graph.plot(x, f(x, phi0), "b-")
graph.plot(x, np.zeros(x.shape), "b:",linewidth=0.5)

for n in range(1, len(path)):
    xnn = path[n-1]
    xn = path[n]

    graph.plot((xnn, xnn), (0, f(xnn, phi0)), "k:", linewidth=1)
    graph.plot((xnn, xn), (f(xnn, phi0), 0), "r-", linewidth=2)
    if n < 5:
        graph.text(xnn, -0.06, "$x_%d$" % (n-1),
                   horizontalalignment='center')

graph.axis([rmin, rmax,-0.25,0.5])
graph.set_xticks(np.arange(rmin, rmax+0.01, 0.25))

figure.savefig("newton.pdf")
