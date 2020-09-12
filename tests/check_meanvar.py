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
# Test der Mittelwerte
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import meanvar

y = np.random.uniform(size=1000)

mean, var, err = meanvar.meanvar(y)

np.testing.assert_allclose(mean, 0.5, atol=0.01)
np.testing.assert_allclose(var, 1 / 12, atol=0.01)
np.testing.assert_allclose(err, 1 / 12 / np.sqrt(len(y)), atol=0.01)
