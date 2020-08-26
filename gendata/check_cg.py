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

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

from cg import conjugate_gradient

A = np.random.uniform(0, 1, 100)
A = A.reshape((10, 10))
A = 0.5 * (A + A.transpose())

b = np.random.uniform(0, 1, 10)

res = np.linalg.solve(A, b)

np.testing.assert_allclose(
    conjugate_gradient(
        A,
        b,
        np.zeros_like(b)),
    res,
    atol=1e-5,
    err_msg="Conjugate Gradient hat das Minimum nicht gefunden")
