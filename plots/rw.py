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
# Der Klassiker: Random Walk
#
############################################
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

np.random.seed(123)

# Anzahl Schritte
N = 100

##########################################

figure = plt.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.15,wspace=0.3, left=0.1,right=0.95)

##########################################

# random walk
def rw(N):
    def move(x):
        if x == 0: return -1
        elif x == 1: return 1
    moves = np.array([ move(x) for x in np.random.randint(0,2,N)])
    positions = np.concatenate(([0], np.cumsum(moves, axis=0)))
    return positions

graph = figure.add_subplot(121)

for lstyle, mstyle in (("r-", "ro"), ("b:", "bD"), ("g--", "gx")):
    pos = rw(N)
    graph.plot(range(N+1), pos, lstyle)
    graph.plot(range(N+1), pos, mstyle, markersize=3, linewidth=0)

graph.xaxis.set_label_text("$t$")
graph.yaxis.set_label_text("$x_t$")

##########################################

samples=1000
avgpos = np.zeros(N+1)
p5 = []
p80 = []
for cnt in range(samples):
    r = rw(N)
    avgpos += r**2
    p5.append(r[5])
    p80.append(r[80])

avgpos /= samples

graph = figure.add_subplot(222)

graph.plot(range(N+1), avgpos, "ro", markersize=3, linewidth=0)
graph.plot(range(N+1), range(N+1), "k--", linewidth=1)

graph.xaxis.set_label_text("$t$")
graph.yaxis.set_label_text("$<x_t^2>$")

##########################################

p5histo, edges = np.histogram(p5, bins=20, range=(-20,20), density=True)
p80histo, edges = np.histogram(p80, bins=20, range=(-20,20), density=True)

graph = figure.add_subplot(224)

width = (edges[1] - edges[0])
graph.bar(edges[:-1], p80histo, width=width,
          color="white", edgecolor="blue", linewidth=1)
graph.bar(edges[:-1], p5histo, width=0.9*width,
          color="#ff8080", edgecolor="red", linewidth=1)

x=np.linspace(-20,20,100)
for s2 in (5, 80):
    graph.plot(x, 1/np.sqrt(2*np.pi*s2)*np.exp(-0.5*x**2/s2), "k--",linewidth=1)

graph.xaxis.set_label_text("$x_t$")
graph.yaxis.set_label_text("$P(x_t)$")

figure.savefig("rw.pdf")
