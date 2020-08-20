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
# Chi^2-Test, statistischer Test von RNGs
#
############################################
import matplotlib.pyplot as plt
import numpy as np

from rng_tests import *

# Anzahl Samples
N = 100

# erwartete Varianz der Verteilung 
s = 1.0/np.sqrt(12*N)

# Darstellungsbereich
rang=(0.35, 0.6501)

def normal(x):
    return 1/np.sqrt(2*np.pi)/s*np.exp(-0.5*(x-0.5)**2/s**2)

##########################################

figure = plt.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.15,wspace=0.3, left=0.1,right=0.95)

##########################################

rng = Rand()
data = [ sum([rng.next() for k in range(N)])/N for x in range(10000)]

histo, edges = np.histogram(np.array(data), bins=50, range=rang, density=True)

graph = figure.add_subplot(121)

x = linspace(rang[0], rang[1],100)

graph.plot(x, normal(x), "r--",linewidth=2)
graph.bar(edges[:-1], histo, width=(edges[1] - edges[0]),
          color="white", edgecolor="blue", linewidth=2)
graph.axis(rang + (0, 14))

##########
# chi^2

histo, edges = np.histogram(np.array(data), bins=50, range=rang)

chi = 0
for i in range(len(histo)):
    # Mittelpunktregel
    E = len(data)*(edges[i+1] - edges[i])*normal(0.5*(edges[i] + edges[i+1]))
    chi += (histo[i] - E)**2 / E

print("chi^2 for rand:", chi)

##########################################

rng = Minstd()
data = [ sum([rng.next() for k in range(N)])/N for x in range(10000)]

histo, edges = np.histogram(np.array(data), bins=50, range=rang, density=True)

graph = figure.add_subplot(122)

x = np.linspace(rang[0], rang[1],100)

graph.plot(x, normal(x), "r--",linewidth=2)
graph.bar(edges[:-1], histo, width=(edges[1] - edges[0]),
          color="white", edgecolor="blue", linewidth=2)
graph.axis(rang + (0, 14))

##########
# chi^2

histo, edges = np.histogram(array(data), bins=50, range=rang, normed=False)

chi = 0
for i in range(len(histo)):
    # Mittelpunktregel
    E = len(data)*(edges[i+1] - edges[i])*normal(0.5*(edges[i] + edges[i+1]))
    chi += (histo[i] - E)**2 / E

print("chi^2 for minstd:", chi)

##########################################

figure.savefig("statistics_test.pdf")
