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
# Verwerfungsmethode fuer rho= x^2
#
############################################
import math
import matplotlib.pyplot as plt
import numpy as np

def rho(x):
    return x**2

def verwerfungsmethode(rho, n):
    """Ziehen eines rho-verteilten Zufallsvektors mit n Dimensionen.
    rho muss eine Pythonfunktion mit einem n-dimensionalen Vektor als
    einziges Argument.
    """
    global cnt
    while True:
        cnt += 1
        p = np.random.uniform(0, 1, n)
        u = np.random.uniform(0, 1)
        if u < rho(p):
            return p

##########################################

figure = plt.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.15,wspace=0.3, left=0.1,right=0.95)

##########################################

graph = figure.add_subplot(121)

x = np.linspace(0,1,100)

graph.plot(x, rho(x), "k-")

pts = zip(np.random.uniform(0,1,200), np.random.uniform(0,1,200))
inner = []
outer = []
for xx, yy in pts:
    if rho(xx) > yy:
        inner.append((xx,yy))
    else:
        outer.append((xx,yy))

xx, yy = zip(*outer)
graph.plot(xx, yy, "o", markeredgecolor="red", markerfacecolor="red")
xx, yy = zip(*inner)
graph.plot(xx, yy, "*", markeredgecolor="green", markerfacecolor="green")

##########################################

drawn = []
cnt = 0
for n in range(100000):
    drawn.append(verwerfungsmethode(rho, 1))

print("ZZ pro Zug", float(cnt)/100000)

histo, edges = np.histogram(np.array(drawn), bins=50, range=(0,1), density=True)

graph = figure.add_subplot(122)

graph.bar(edges[:-1], histo, width=(edges[1] - edges[0]),
          color="blue", edgecolor="blue")
graph.plot(x, 3*rho(x), "r--", linewidth=2)
graph.yaxis.set_label_text(r"$\rho(x)$")
graph.xaxis.set_label_text(r"$x$")

figure.savefig("rejection.pdf")
