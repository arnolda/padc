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

seed(123)

ex = "rand"

if ex == "fibo":
    # Fibonacci-Matrix
    A = array(((0,1),(1,1)))
    # was rauskommen soll
    l1 = 0.5*(1+sqrt(5))
    l2 = 0.5*(1-sqrt(5))
    v1 = array([1, l1])
    v1 = v1/norm(v1)
    v2 = array([1, l2])
    v2 = v2/norm(v2)

    # Anzahl Schritte QR
    n = 6
    # Anzahl Schritte inverse Iteration
    ni = 3
else:
    A = uniform(0,1,10*10)
    A = A.reshape((10, 10))
    # symmetrisieren fuer reelle Eigenwerte
    for i in range(10):
        for k in range(10):
            A[i, k] = A[k,i]
    n = 200
    ni = 5

    v1 = 0
    v2 = 0

I = identity(A.shape[0])


def digits(x, tgt):
    if tgt - x == 0:
        return 17
    return int(floor(-log(abs(tgt-x))/log(10)))

# QR-Algorithmus
################################

Ak = A.copy()
dim = Ak.shape[0]
for i in range(n):
    if ex == "fibo":
        sys.stdout.write("%d & $\\begin{pmatrix}\n" % i)
        sys.stdout.write("% 10.8f & % 10.8f\\\\\n" % (Ak[0,0], Ak[0,1]))
        sys.stdout.write("% 10.8f & % 10.8f\n" % (Ak[1,0], Ak[1,1]))
        sys.stdout.write("\\end{pmatrix}$\n")

        sys.stdout.write("& %d & %d\\\\\n" %
                         (digits(Ak[0,0], l2),
                          digits(Ak[1,1], l1)))
    else:
        print max([abs(Ak[i,k]) for i in range(A.shape[0]) for k in range(i)])
                
    shift = Ak[-1,-1]
    Q, R = qr(Ak - shift*I)
    Ak = dot(R, Q) + shift*I

# inverse Iteration
print "Inverse Iteration"

for l,v in (Ak[0,0], v2), (Ak[1,1], v1):
    x = ones(A.shape[0])
    Ashift = A - l*I
    for i in range(ni):
        x = solve(Ashift, x)
        x = x / norm(x)
        if ex == "fibo":
            if x[0] < 0:
                xp = -x
            else:
                xp = x
                sys.stdout.write("% 10.8f & % 10.8f & %d & %d\\\\\n" % \
                                     (xp[0],xp[1], digits(xp[0], v[0]),
                                      digits(xp[0], v[0])))
        else:
            print norm(l*x - dot(A, x))/norm(x)
