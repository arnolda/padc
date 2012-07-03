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
# CG-Verfahren
#
############################################
from scipy import *
from scipy.linalg import *
import matplotlib.pyplot as pyplot
import matplotlib.colors as colors

def conjugate_gradient(A, b, x0, tol=1e-10):
    x = x0
    path = [ x.copy() ]
    r = b - dot(A, x)
    r2 = dot(r,r)
    d = r.copy()
    while sqrt(r2) > tol:
        lmbda = dot(d, r) / dot(d, dot(A, d))
        x += lmbda*d
        path.append(x.copy())
        r -= lmbda*dot(A,d)
        oldr2 = r2
        r2 = dot(r,r)
        # orthogonalisieren
        d = r + r2/oldr2 * d


    return array(path).transpose()

############################################

figure = pyplot.figure(figsize=(4,4))
figure.subplots_adjust(left=0.05, right=0.95)
pyplot.gray()

A = array([[70,30],[30,50]])
b = array([0,0])

def f(x):
    return 0.5*dot(x, dot(A, x)) + dot(b, x)

graph = figure.add_subplot(111)

axis = (-1,1,-1,1)

path = conjugate_gradient(A, b, array((-0.5, -1)))

# contour
n=100
rx = linspace(axis[0], axis[1], n)
rx = linspace(axis[2], axis[3], n)
x, y = meshgrid(rx, rx)
z = zeros_like(x)
for k in range(n):
    for l in range(n):
        z[k, l] = max(1e-20, f([ x[k, l], y[k, l] ]))

levels = (1,10,30,50,100)
cont = graph.contour(x, y, z, levels=levels,
                     norm=colors.LogNorm(min(levels), max(levels)*2),linewidths=0.5)
graph.clabel(cont, inline=1, fontsize=10, fmt="%g")

# target cross
graph.plot(1, 1 , "+", markeredgecolor="black",markersize=20)

# path
graph.plot(path[0], path[1] , "-",color="#ffa0a0",linewidth=1)
graph.plot(path[0], path[1] , "o",
           markeredgecolor="red",markerfacecolor="red",markersize=3)

graph.axis(axis)

figure.savefig("cg.pdf")
