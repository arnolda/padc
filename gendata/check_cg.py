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
# Test des CG-Verfahrens
##############################################

from scipy import *
from scipy.linalg import *
from numpy.random import *

import sys
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

from cg import conjugate_gradient

A = uniform(0, 1, 100)
A = A.reshape((10,10))
A = 0.5*(A + A.transpose())

b = uniform(0, 1, 10)

res = solve(A, b)

if norm(conjugate_gradient(A, b, zeros_like(b)) - res) > 1e-5:
    raise Exception("Conjugate Gradient hat das Minimum nicht gefunden")
