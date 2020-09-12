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
# Test der Verwerfungsmethode
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

def rho(x):
    return x**2

from verwerfungsmethode import verwerfungsmethode

res = np.array([verwerfungsmethode(rho, 1) for _ in range(100000)])

histo, edges = np.histogram(res, bins=50, range=(0, 1), density=True)

np.testing.assert_allclose(histo, 3/4*(edges[1:] + edges[:-1])**2, atol=1e-1)
