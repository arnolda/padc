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
# Test der Hornerauswertung fuer die Newtonsche Darstellung
###########################################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import horner_newton_interpol
import neville

x = np.arange(0, 1, 0.1)

gamma = neville.neville(x, 2*x**3 + x)

xprobe = np.arange(0, 1, 0.01)
res = horner_newton_interpol.horner_newton(xprobe, x, gamma)

np.testing.assert_allclose(res, 2*xprobe**3 + xprobe, atol=0.01)
