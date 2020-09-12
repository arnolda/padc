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
# Test des gedaempften Newtonverfahrens
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")


def f(x):
    return np.cos(x)


def df(x):
    return np.array((-np.sin(x),))


from damped_newton import damped_newton

# ueberschiesst ohne Daempfung
res = damped_newton(f, df, np.array([0.4]), 1e-6)

np.testing.assert_allclose(res, np.pi / 2, atol=1e-2)
