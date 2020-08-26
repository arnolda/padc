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
# Illustration Box-Muller-Verfahren
#
############################################
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np

np.random.seed(123)

cm_x = cm.gray
cm_y = cm.hsv

def mix_color(c1, c2):
    return np.sqrt(np.array(colors.colorConverter.to_rgb(c1)) *
                np.array(colors.colorConverter.to_rgb(c2)))

# random points
pts = zip(np.random.uniform(0,1,1000), np.random.uniform(0,1,1000))

transforms = []
for u,v in pts:
    col = mix_color(cm_x(1-u), cm_y(v))
    r = np.sqrt(-2*np.log(u))
    phi = 2*np.pi*v
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    transforms.append((u, v, x, y, r, phi, col))

# markers
# avoid zero because of log
nspace = np.linspace(1e-4,1,100)
markers = [
    (0.1*np.ones_like(nspace), nspace, "m--"),
    (0.5*np.ones_like(nspace), nspace, "r--"),
    (0.8*np.ones_like(nspace), nspace, "k--"),
    (nspace, 0.5*np.ones_like(nspace),  "b-") ]

marker_transforms = []
for m in markers:
    transform = []

    uvec, vvec, col = m
    for u,v in zip(uvec, vvec):
        r = np.sqrt(-2*np.log(u))
        phi = 2*np.pi*v
        x = r*np.cos(phi)
        y = r*np.sin(phi)
        transform.append((u, v, x, y, r, phi))

    # back to uvec, vvec, xvec,...
    transform = zip(*transform)
    transform = list(transform)
    transform.append(col)
    marker_transforms.append(transform)

##########################################

figure = plt.figure(figsize=(8,8))
figure.subplots_adjust(bottom=0.15,wspace=0.3, left=0.1,right=0.95)

##########################################

graph = figure.add_subplot(221)

for u, v, x, y, r, phi, col in transforms:
    graph.plot(u, v, ".", markeredgecolor=col, markerfacecolor=col,
               markersize=3)

for m in marker_transforms:
    u, v, x, y, r, phi, col = m
    graph.plot(u, v, col, linewidth=2)

graph.xaxis.set_label_text(r"$u$")
graph.yaxis.set_label_text(r"$u'$")

##########################################

graph = figure.add_subplot(222)

for u, v, x, y, r, phi, col in transforms:
    graph.plot(x, y, ".", markeredgecolor=col, markerfacecolor=col,
               markersize=3)

for m in marker_transforms:
    u, v, x, y, r, phi, col = m
    graph.plot(x, y, col, linewidth=2)

graph.xaxis.set_label_text(r"$x$")
graph.yaxis.set_label_text(r"$y$")
graph.axis((-3,3,-3,3))

##########################################

graph = figure.add_subplot(223)

for u, v, x, y, r, phi, col in transforms:
    graph.plot(r, phi, ".", markeredgecolor=col, markerfacecolor=col,
               markersize=3)

for m in marker_transforms:
    u, v, x, y, r, phi, col = m
    graph.plot(r, phi, col, linewidth=2)

graph.xaxis.set_label_text(r"$r$")
graph.yaxis.set_label_text(r"$\phi$")
graph.axis((0,4, 0, 2*np.pi))

##########################################

z1, z2 = list(zip(*transforms))[2:4]
 
histox, edges = np.histogram(np.array(z1), bins=20, range=(-3,3), density=True)
histoy, edges = np.histogram(np.array(z2), bins=20, range=(-3,3), density=True)

graph = figure.add_subplot(224)

w = 0.5*(edges[1] - edges[0])
graph.bar(edges[:-1], histox, width=w,
          color="blue", edgecolor="blue")
graph.bar(edges[:-1] + w, histoy, width=w,
          color="green", edgecolor="green")
x = np.linspace(-3,3,100)
graph.plot(x, 1/np.sqrt(2*np.pi)*np.exp(-0.5*x**2), "r--", linewidth=2)

graph.xaxis.set_label_text(r"$z$")
graph.yaxis.set_label_text(r"$P(z)$")
graph.axis((-3,3,0,1))
figure.savefig("inversion.pdf")
