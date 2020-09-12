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
# Test der van der Corput Quasi-ZZ
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")


from vander_corput import vanderCorput

res = vanderCorput(8, 2)

assert np.all(res == [1/2, 1/4, 3/4, 1/8, 5/8, 3/8, 7/8, 1/16])

res = vanderCorput(8, 3)

np.testing.assert_allclose(res, [1/3, 2/3, 1/9, 4/9, 7/9, 2/9, 5/9, 8/9])
