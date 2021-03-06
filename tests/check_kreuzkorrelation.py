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
# Test der Kreuzkorrelation
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import kreuzkorrelation

x = np.arange(0, 2*np.pi, 0.1)
a = np.sin(x)
b = np.cos(x)

np.testing.assert_allclose(kreuzkorrelation.kreuzkorrelation(a, b), -0.5*np.sin(x), atol=0.01)
np.testing.assert_allclose(kreuzkorrelation.kreuzkorrelation(a, a), 0.5*np.cos(x), atol=0.01)
