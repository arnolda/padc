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
# Test des Hornerschemas
##############################################

import numpy as np
import sys
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

# check our horner against numpy's
#########################################

from horner import horner

p = [1, 3, 42, -3]
pr = p[:]
pr.reverse()

for x in np.arange(0, 1, 0.02):
    mine = horner(p, x)
    theirs = np.polyval(pr, x)
    np.testing.assert_almost_equal(
        mine,
        theirs,
        decimal=10,
        err_msg=f"horner: unexpected difference at {x}: {mine} != {theirs}")

# check horner newton step
#########################################

from horner_newton import poly_newton_step

p = [-2, 0, 1]

x = 1.0
for i in range(10):
    x = poly_newton_step(p, x)

np.testing.assert_almost_equal(
    x**2,
    2,
    decimal=10,
    err_msg="Newton step did not converge!")
