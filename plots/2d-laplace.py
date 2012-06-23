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
# Matrix zum 2d-Laplace-Operator
##############################################

from scipy import *
from scipy.linalg import *
import matplotlib.pyplot as pyplot

# Laplace 1. Ordnung auf NxN
N=5

def linindex(x, y): return (x % N) + N*(y % N)

h=1

# matrix
A=zeros((N*N, N*N))
# and corresponding points to plot
# 4/h^2
M = []
# -1/h^2
O = []
# 1 (just 1 point for gauge)
P = []

eqn=0
for y in range(N):
    for x in range(N):
        if eqn > 0:
            A[eqn, linindex(x,y)] = -4/h**2
            A[eqn, linindex(x+1,y)] = 1/h**2
            A[eqn, linindex(x-1,y)] = 1/h**2
            A[eqn, linindex(x,y+1)] = 1/h**2
            A[eqn, linindex(x,y-1)] = 1/h**2
            M.append((linindex(x,y)  , eqn))
            O.append((linindex(x+1,y), eqn))
            O.append((linindex(x-1,y), eqn))
            O.append((linindex(x,y+1), eqn))
            O.append((linindex(x,y-1), eqn))
            if eqn <= 2:
                # kommt in den Text als Beispiel der indizierten Speicherung
                print linindex(x+1,y), linindex(x-1,y),\
                    linindex(x,y+1), linindex(x,y-1)
                print A[eqn,:]
        else:
            # Normierung
            A[0, 0] = 1
            P.append((0, 0))

        eqn += 1

print det(A)

# Graphik
##########################################

figure = pyplot.figure(figsize=(4,4))
figure.subplots_adjust(left=0.1, bottom=0.05, right=1, top=0.95)

##########################################

graph = figure.add_subplot(111)

mx, my = zip(*M)
ox, oy = zip(*O)
px, py = zip(*P)

for y in range(N):
    graph.plot((y*N, y*N), (0,N*N), ":", color="#a0a0a0")

graph.plot(mx, my, "o", markersize=2,
           markerfacecolor="black", markeredgecolor="black")
graph.plot(ox, oy, "D", markersize=2,
           markerfacecolor="green", markeredgecolor="green")
graph.plot(px, py, "+", markersize=6,
           markerfacecolor="red", markeredgecolor="red")
graph.axis((-0.5, N*N-0.5, N*N-0.5, -0.5))

figure.savefig("2d-laplace.pdf")
