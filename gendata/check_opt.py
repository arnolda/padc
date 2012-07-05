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
# Test der Optimierungsverfahren
##############################################

from scipy import *
from scipy.linalg import *
from numpy.random import *

import sys
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

def rosenbrock(x):
    return (1.0-x[0])**2 + 100.0*(x[1]-x[0]**2)**2


def gradrosenbrock(x):
    return array((2*(200*x[0]**3 - 200*x[0]*x[1] + x[0] - 1),
                  200*(x[1]-x[0]**2)))

# Armijo
#########################################

from armijo import armijo_steepest_descent

if norm(armijo_steepest_descent(rosenbrock, gradrosenbrock, array((0,0)), tol=1e-6) - array((1,1))) > 1e-5:
    raise Exception("Steepest Descent hat das Minimum nicht gefunden")

# CG
#########################################

from cg import conjugate_gradient

A = uniform(0, 1, 100)
A = A.reshape((10,10))
A = 0.5*(A + A.transpose())

b = uniform(0, 1, 10)

res = solve(A, b)

if norm(conjugate_gradient(A, b, zeros_like(b)) - res) > 1e-5:
    raise Exception("Conjugate Gradient hat das Minimum nicht gefunden")
