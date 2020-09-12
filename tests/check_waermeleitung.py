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
# Test des Romberg-Verfahrens
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import waermeleitung

N = waermeleitung.N = 21
L = waermeleitung.L

tnpns = waermeleitung.solve(waermeleitung.init_density(), waermeleitung.f, tmax=250)

np.testing.assert_allclose(tnpns[-1][1:], 0, atol=0.01)


def fhomogen(t, p):
    diff = waermeleitung.D*np.dot(waermeleitung.laplace_operator(), p)
    # Delta-Quelle bei L/2
    diff[N//2] += N / L
    return diff


tnpns = waermeleitung.solve(np.zeros(N), fhomogen, tmax=500)

x = 0.5 * L / N + np.linspace(0, L, N, endpoint=False)

expected = 10.5 - abs(x - 0.5 * L)

np.testing.assert_allclose(tnpns[-1][1:], expected, atol=0.1)
