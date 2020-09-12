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

import romberg

pi =  4*romberg.romberg(lambda x: 1.0/(1.0 + x**2), 0.0, 1.0, 8)

np.testing.assert_allclose(pi, np.pi)
