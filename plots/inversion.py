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
from scipy import *
from numpy.random import *
import math
import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
import matplotlib.cm as cm

seed(123)

cm_x = cm.gray
cm_y = cm.hsv

def mix_color(c1, c2):
    return sqrt(array(colors.colorConverter.to_rgb(c1)) *
                array(colors.colorConverter.to_rgb(c2)))

# random points
pts = zip(uniform(0,1,1000), uniform(0,1,1000))

transforms = []
for u,v in pts:
    col = mix_color(cm_x(1-u), cm_y(v))
    r = sqrt(-2*log(u))
    phi = 2*pi*v
    x = r*cos(phi)
    y = r*sin(phi)
    transforms.append((u, v, x, y, r, phi, col))

# markers
# avoid zero because of log
nspace = linspace(0.0001,1,100)
markers = [
    (0.1*ones_like(nspace), nspace, "m--"),
    (0.5*ones_like(nspace), nspace, "r--"),
    (0.8*ones_like(nspace), nspace, "k--"),
    (nspace, 0.5*ones_like(nspace),  "b-") ]

marker_transforms = []
for m in markers:
    transform = []

    uvec, vvec, col = m
    for u,v in zip(uvec, vvec):
        r = sqrt(-2*log(u))
        phi = 2*pi*v
        x = r*cos(phi)
        y = r*sin(phi)
        transform.append((u, v, x, y, r, phi))

    # back to uvec, vvec, xvec,...
    transform = zip(*transform)
    transform.append(col)
    marker_transforms.append(transform)

##########################################

figure = pyplot.figure(figsize=(8,8))
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
graph.axis((0,4, 0, 2*pi))

##########################################

z1, z2 = zip(*transforms)[2:4]
 
histox, edges = histogram(array(z1), bins=20, range=(-3,3), normed=True)
histoy, edges = histogram(array(z2), bins=20, range=(-3,3), normed=True)

graph = figure.add_subplot(224)

w = 0.5*(edges[1] - edges[0])
graph.bar(edges[:-1], histox, width=w,
          color="blue", edgecolor="blue")
graph.bar(edges[:-1] + w, histoy, width=w,
          color="green", edgecolor="green")
x = linspace(-3,3,100)
graph.plot(x, 1/sqrt(2*pi)*exp(-0.5*x**2), "r--", linewidth=2)

graph.xaxis.set_label_text(r"$z$")
graph.yaxis.set_label_text(r"$P(z)$")
graph.axis((-3,3,0,1))
figure.savefig("inversion.pdf")
