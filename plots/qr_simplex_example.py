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
# Polynom-Fit QR/Simplex
##############################################
from scipy import *
from scipy.linalg import *
from numpy.random import *
import matplotlib.pyplot as pyplot
import sys
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

from householder import householder
from simplex import simplex

# Fit eines Polynoms n-ten Grades an exp(x)
#############################################

n = 3
def f(x): return exp(x)
support = linspace(-1, 1, 10)

# inverse Reihenfolge fuer numpy's polyval
A = array([[x**(n-k-1) for k in range(n)] for x in support])
b = array([f(x) for x in support])

# Least-squares fit via QR
#############################################

Q, R = householder(A)

# oberes Quadrat
R1 = R[:R.shape[1]]
cqr = solve(R1.transpose(), dot(A.transpose(), b))
cqr = solve(R1, cqr)

# Chebyshev, Simplex
#############################################

e = ones(A.shape[0])
e = e.reshape((A.shape[0], 1))
z = zeros((A.shape[0], A.shape[0]))
I = identity(A.shape[0])
# anders als im Skript sind die Maximumsvariablen vorne
# macht es einfacher, sie loszuwerden
Asupper = concatenate((e, -e,  A, -A, -I,  z), axis=1)
Aslower = concatenate((e, -e, -A,  A,  z, -I), axis=1)
As = concatenate((Asupper, Aslower))
bs = concatenate((b, -b))
cs = zeros(As.shape[1])
cs[0] =  1
cs[1] = -1

csimplex = simplex(cs, As, bs)[2:2 + n]

# Ausgabe
#############################################

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, right=0.95)

# Fits
#############################################
graph = figure.add_subplot(121)

x = linspace(min(support), max(support), 200)
graph.plot(x, f(x), "k-", linewidth=2)
graph.plot(support, f(support), "ko", markersize=5)

graph.plot(x, polyval(cqr, x), "r--")
graph.plot(x, polyval(csimplex, x), "g:")

# Fehler
#############################################
graph = figure.add_subplot(122)

x = linspace(min(support), max(support), 200)
graph.plot(x, abs(polyval(cqr, x) - f(x)) , "r--", linewidth=2)
graph.plot(support, abs(polyval(cqr, support) - f(support)) , "ro",
           markersize=5, clip_on=False)
graph.plot(x, abs(polyval(csimplex, x) - f(x)), "g:", linewidth=2)
graph.plot(support, abs(polyval(csimplex, support) - f(support)),
           "gD", markersize=5, clip_on=False)

figure.savefig("qr_simplex_example.pdf")
pyplot.show()
