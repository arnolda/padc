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
# Test Bessel-DGL
##############################################

import sys
from scipy import special
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import bessel

expected = special.jn(bessel.nu, bessel.x)

np.testing.assert_allclose(bessel.stencil3(), expected, atol=0.4)
np.testing.assert_allclose(bessel.stencil3l(), expected, atol=0.25)
np.testing.assert_allclose(bessel.stencil5(), expected, atol=0.1)
