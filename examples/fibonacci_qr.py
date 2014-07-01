# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Demonstration QR-Algorithmus
######################################

from scipy import *
from scipy.linalg import *
from numpy.random import *
import sys

# Fibonacci-Matrix
A = array(((0,1),(1,1)))

# QR-Algorithmus
################################
print "QR-Algorithmus"

Ak = A.copy()
I = identity(Ak.shape[0])
for i in range(6):
    print "Iteration {}, Ak ist".format(i)
    print Ak

    shift = Ak[-1,-1]
    Q, R = qr(Ak - shift*I)
    Ak = dot(R, Q) + shift*I

print "Eigenwerte sollten {} und {} sein.".format(0.5*(1-sqrt(5)), 0.5*(1+sqrt(5)))

# Inverse Iteration
################################

# Iteration ueber die Eigenwerte
for l in diag(Ak):
    print "Inverse Iteration zum Eigenwert {}".format(l)
    Ashift = A - l*I
    # Startwert
    xk = ones(A.shape[0])
    for i in range(2):
        xk = solve(Ashift, xk)
        xk = xk / norm(xk)
        print "Iteration {}, xk ist".format(i)
        print xk

    print "Der Eigenvektor sollte sein:"
    v = array([1.0, l])
    print v/norm(v)
