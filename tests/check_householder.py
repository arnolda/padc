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
# Test Householder
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import householder

n = 10
A = np.random.uniform(0, 1, n * n).reshape((n, n))

q, r = householder.householder(A)

np.testing.assert_allclose(q.conj().T.dot(q), np.identity(n), atol=1e-10)
np.testing.assert_allclose(q.dot(r), A)

rrest = r.copy()

for i in range(n):
    for j in range(i, n):
        rrest[i, j] = 0

np.testing.assert_allclose(rrest, np.zeros((n, n)), atol=1e-10)
