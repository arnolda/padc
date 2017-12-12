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
# mit Quasizufallszahlen
#
############################################
from scipy import *
from numpy.random import *
import math
import matplotlib.pyplot as pyplot

seed(123)

def vanderCorput(N, p):
    # zu wandelnde Zahlen
    numbers = arange(1,int(N)+1)
    # bitumgekehrtes Ergebnis
    result = zeros(N)
    # Wert der aktuellen, inversen Stelle
    frac = 1.0 / p

    # solange die groesste Zahl noch Stellen hat
    while numbers[-1] > 0:
        # unterste Stelle abschneiden
        digit = numbers % p
        numbers /= p
        # ... und zum Ergebnis hinzufuegen
        result += frac*digit
        frac /= p

    return result

def f(x, y):
    return 1.0*(x**2 + y**2 <= 1)

def trapez(f, N_list):
    res = []
    for Ntgt in N_list:
        N = int(ceil(Ntgt**0.5))
        h = 2.0/N
        pos = arange(0, N)*h + h/2 - 1.0
        S = 0
        for y in pos:
            S += sum(f(pos, y))

        res.append(S*(2.0/N)**2)

    return array(res)

def mc(f, N_list):
    res = []
    for N in N_list:
        res.append(2.0**2*sum(f(uniform(-1,1,N), uniform(-1,1,N)))/N)
    return array(res)

def qmc(f, N_list):
    res = []
    for N in N_list:
        res.append(2.0**2*sum(f(-1 + 2*array(vanderCorput(N,2)),
                                 -1 + 2*array(vanderCorput(N,3))))/N)
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

N = 200

x = -1 + 2*array(vanderCorput(N,2))
y = -1 + 2*array(vanderCorput(N,3))

# Punkte nach der x-Reihefolge markieren
color = []
for c in range(10):
    color.append("red")
for c in range(10,110):
    color.append("green")
for c in range(110,N):
    color.append("blue")

graph.scatter(x, y, c=color,linewidth=0)

graph.axis((-1,1,-1,1))

##########################################

graph = figure.add_subplot(122)

N_range = logspace(0, 6, 10, dtype=int)

graph.set_xscale("log")
graph.set_yscale("log")
graph.plot(N_range, abs(trapez(f, N_range) - pi), "gx", markersize=3)
graph.plot(N_range, avg(lambda: abs(mc(f, N_range) - pi), 20), "bo", markersize=3)
graph.plot(N_range, abs(qmc(f, N_range) - pi), "rD", markersize=3)

graph.set_xlabel("N")
graph.axis((10,1.1e6, 1e-5, 1))

figure.savefig("quasirandom.pdf")
