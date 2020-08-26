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
# LJ-Energiebeispiel
##############################################

import sys

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
import scipy.optimize

sys.path.append("..")

# Quadratkantenlaenge
L=10

# LJ-Potential, reduzierte Einheiten und quadratische Distanz als Eingabe
def lj(r2):
    return 4*(r2**-6 - r2**-3)

def minimag(p0, positions):
    "Min. Imag.-Trafo fuer einen Vektor gegen einen Satz Vektoren"
    d = positions.transpose().copy()
    for i in range(p0.shape[0]):
        d[i] -= p0[i]
        d[i] -= np.around(d[i]/L)*L
    return d

# Potential am Ort eines Teilchens
def energy(p0, particles):
    ds = minimag(p0, particles)
    dds = ds[0]*ds[0] + ds[1]*ds[1]
    return sum(lj(dds))

figure = plt.figure(figsize=(8,4))

# LJ-Potential
##################################

graph = figure.add_axes([0.1,0.1,0.375,0.75])

r = np.linspace(0.5, 2.5, 200)

m = 2**(1./6)
graph.plot(r, lj(r*r), "r-", linewidth=2)
graph.plot((m-0.1,m+0.1), (-1,-1), "k-", linewidth=1)
graph.plot((0,2.5), (0,0), "k-", linewidth=1)
graph.plot((m,m), (0,-1), "k--", linewidth=1)
graph.plot((1,1), (-1.5,2), "k--", linewidth=1)
graph.annotate("$\epsilon$", (m+0.02, -0.5), ha="left", size=16)
graph.annotate("$\sigma$", (1.02, 0.02), ha="left", va="bottom", size=16)

graph.axis((0,2.5,-1.5,2))
graph.xaxis.set_label_text("$r/\sigma$")
graph.yaxis.set_label_text("$\phi_{LJ}(r)/\epsilon$")

# Energielandschaft 
##################################

graph = figure.add_axes([0.52,0.1,0.475,0.75])

positions = np.loadtxt("lj_snapshot.coord")

# heatmap
n=100
rx = np.linspace(0, L, n)
ry = np.linspace(0, L, n)
x, y = np.meshgrid(rx, ry)
z = np.zeros_like(x)
cnt = 0
candidates = []
for k in range(n):
    for l in range(n):
        z[k, l] = energy(np.array((x[k, l], y[k, l])), positions)
        if z[k, l] < 1000:
            candidates.append(np.array((x[k, l], y[k, l])))
z = np.array(z)
print("Energiebereich", z.min(), z.max())

minima = []
for x0 in candidates:
    minimum = scipy.optimize.fmin(energy, x0, args=(positions,), disp=False)
    for mm in minima:
        if np.allclose(minimum, mm, 1e-2):
           break
    else:
        minima.append(minimum)
print(len(minima), "Minima gefunden")

im = graph.imshow(z, interpolation='bilinear', cmap=cm.gray,
                  origin='lower', norm = colors.Normalize(0,5000),
                  extent=[0,L,0,L])

# Position der anderen Teilchen
x = positions[:,0] % L
y = positions[:,1] % L
graph.plot(x, y, "r.", markersize=4)

# lokale Minima
x, y = zip(*minima)
graph.plot(x, y, "+", markeredgecolor="yellow",
           markerfacecolor="yellow", markersize=6)

figure.colorbar(im, shrink=0.8, extend="max")

figure.savefig("lj.pdf")
