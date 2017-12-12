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
# pi aus der Integration der Indikatorfunktion des Kreises ueber [-1,1]^2
#
############################################
from scipy import *
from numpy.random import *
import math
import matplotlib.pyplot as pyplot

def f2d(x, y):
    return 1.0*(x**2 + y**2 <= 1)

def f5d(x1, x2, x3, x4, x5):
    return 1.0*(x1**2 + x2**2 + x3**2 + x4**2 + x5**2 <= 1)

def trapez2d(f, Ntgt_list):
    res = []
    for Ntgt in Ntgt_list:
        N = int(ceil(Ntgt**0.5))
        h = 2.0/N
        pos = arange(0, N)*h + h/2 - 1.0
        S = 0
        for y in pos:
            S += sum(f(pos, y))

        res.append(S*(2.0/N)**2)

    return array(res)

def trapez5d(f, Ntgt_list):
    res = []
    for Ntgt in Ntgt_list:
        N = int(ceil(Ntgt**0.2))
        h = 2.0/N
        pos = arange(0, N)*h + h/2 - 1.0
        S = 0
        for x1 in pos:
            for x2 in pos:
                for x3 in pos:
                    for x4 in pos:
                        S += sum(f(x1, x2, x3, x4, pos))

        res.append(S*(2.0/N)**5)

    return array(res)

def mc2d(f, N_list):
    res = []
    for N in N_list:
        res.append(2.0**2*sum(f(uniform(-1,1,N), uniform(-1,1,N)))/N)
    return array(res)

def mc5d(f, N_list):
    res = []
    for N in N_list:
        res.append(2.0**5*sum(f(uniform(-1,1,N), uniform(-1,1,N),
                                uniform(-1,1,N), uniform(-1,1,N),
                                uniform(-1,1,N)))/N)
    return array(res)

def avg(f, N):
    res = f()
    for i in range(1,N):
        res += f()
    return res/N

##########################################

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.15)

##########################################

graph = figure.add_subplot(121)

x = linspace(-1,1,100)
graph.fill_between(x, sqrt(1-x**2), -sqrt(1-x**2), facecolor="#f0f0f5")

pts = zip(uniform(-1,1,200), uniform(-1,1,200))
inner = []
outer = []
for x, y in pts:
    if f2d(x,y) > 0:
        inner.append((x,y))
    else:
        outer.append((x,y))

x, y = zip(*outer)
graph.plot(x, y, "ro")
x, y = zip(*inner)
graph.plot(x, y, "g+")

##########################################

graph = figure.add_subplot(122)

N_range = logspace(0, 6, 10, dtype=int)

graph.set_xscale("log")
graph.set_yscale("log")
graph.plot(N_range, abs(trapez2d(f2d, N_range) - pi), "rx", markersize=3)
graph.plot(N_range, avg(lambda: abs(mc2d(f2d, N_range) - pi), 20), "bo", markersize=3)

res5d = 8*pi*pi/15
graph.plot(N_range, abs(trapez5d(f5d, N_range) - res5d), "g+", markersize=3)
graph.plot(N_range, avg(lambda: abs(mc5d(f5d, N_range) - res5d), 20), "yD", markersize=3)

graph.plot(N_range, N_range**-0.5, "k-",linewidth=0.5)
graph.set_xlabel("N")
graph.axis((10,1e6+1000, 1e-4, 10))
figure.savefig("pi.pdf")
