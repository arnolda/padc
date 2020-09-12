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

from rw import rw

rws = [rw(10) for _ in range(100000)]
assert all(len(rw) == 11 for rw in rws)

end = np.array([rw[-1] for rw in rws])

np.testing.assert_allclose(np.mean(end), 0, atol=1e-1)
np.testing.assert_allclose(np.mean(end**2), 10, atol=1e-1)
