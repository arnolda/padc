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
# Bessel-Naeherung durch num. Diffgleichung
# und Integraldarstellung
#
############################################
from scipy import *
from scipy.linalg import *
from scipy.special import *
import math
import matplotlib.pyplot as pyplot

nu = 0

xmax = 15.0

# Numerisches Differential
############################################

# Anzahl Punkte Differential
N = 300

# Schrittweite Differential
h = xmax/N

# 3-Punkt-Stencil, Dirichlet-Randbedingung
############################################

x3   = arange(1,N+1)*h

# alles erstmal auf Null
b = numpy.zeros(N)
A = numpy.zeros((N, N))

# Startwerte an den Raendern
A[0, 0] = 1
b[0]    = jn(nu, h)
A[1, N -1] = 1
b[1]    = jn(nu, N*h)

# N-2 innere Punkte von 1 bis N-2 hinter den Startwerten, alle auf Null
for n in range(2, N):
    #             x^2 d^2f/dx^2 + x df/dx    + (x^2 - nu^2)f
    A[n, n - 2] = 0.5*n**2      - 0.5*n
    A[n, n - 1] =    -n**2                   + n**2*h**2 - nu**2
    A[n, n]     = 0.5*n**2      + 0.5*n

fx3 = solve(A, b)

# 3-Punkt-Stencil, linker Rand
############################################

# alles erstmal auf Null
b = numpy.zeros(N)
A = numpy.zeros((N, N))

# Startwerte an den Raendern
A[0, 0] = 1
b[0]    = jn(nu, h)
A[1, 1] = 1
b[1]    = jn(nu, 2*h)

# N-2 innere Punkte von 1 bis N-2 hinter den Startwerten, alle auf Null
for n in range(2, N):
    #             x^2 d^2f/dx^2 + x df/dx    + (x^2 - nu^2)f
    A[n, n - 2] = 0.5*n**2      - 0.5*n
    A[n, n - 1] =    -n**2                   + n**2*h**2 - nu**2
    A[n, n]     = 0.5*n**2      + 0.5*n

fx3l = solve(A, b)

# 5-Punkt-Stencil, Dirichlet
############################################

x5   = arange(0,N)*h

# alles wieder auf Null
b = numpy.zeros(N)
A = numpy.zeros((N, N))

# Startwerte an den Raendern, wie zuvor
A[0, 0]     = 1
b[0]        = jn(nu, 0)
A[1, N - 1] = 1
b[1]        = jn(nu, (N-1)*h)

# innere naechste 2 Punkte durch einfache Naeherung
for z,n in ((2,1),(3, N-2)):
    #             x^2 d^2f/dx^2 + x df/dx    + (x^2 - nu^2)f
    A[z, n - 1] = 0.5*n**2      - 0.5*n
    A[z, n    ] =    -n**2                   + n**2*h**2 - nu**2
    A[z, n + 1] = 0.5*n**2      + 0.5*n

# N-4 innere Punkte von 2 bis N-3 hinter den Startwerten, alle auf Null
for n in range(2, N-2):
    z = n + 2
    #              x^2 d^2f/dx^2 + x df/dx   + (x^2 - nu^2)f
    A[z, n - 2] = -1./12*n**2    + 1./12*n
    A[z, n - 1] =  4./3 *n**2    - 2./3 *n
    A[z, n    ] = -5./2 *n**2                + n**2*h**2 - nu**2
    A[z, n + 1] =  4./3 *n**2    + 2./3 *n
    A[z, n + 2] = -1./12*n**2    - 1./12*n

fx5 = solve(A, b)

############################################

figure = pyplot.figure(figsize=(4,4))

x = linspace(0, xmax, 200)

# diff equation
############################################

graph = figure.add_subplot(111)

graph.plot(x, jn(nu, x), "k-",linewidth=0.5)
graph.plot(x3, fx3, "g-.")
graph.plot(x3, fx3l, "r--")
graph.plot(x5, fx5, "b--", dashes=(1,8), dash_capstyle="round", linewidth=3)
graph.axis((0, xmax, -0.5, 1))

figure.savefig("bessel_diff.pdf")
