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
# Test des VV-Loesers
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import vv


def F(_, y):
    a = y
    return -np.sin(a)

periode = 2*np.pi

result = vv.velocity_verlet(F, [0], [0.1], periode, 0.01)
_, final_a, final_da = result[-1]
np.testing.assert_allclose(final_a, 0, atol=0.01)
np.testing.assert_allclose(final_da, 0.1, atol=0.01)
