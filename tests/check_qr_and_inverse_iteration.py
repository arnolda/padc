# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Test QR-Code und inverse Iteration
######################################

import sys
import numpy as np
sys.path.append("..")

from qr import qr_eigenvalues
from inverse_iteration import inverse_iteration

np.random.seed(123)

A = np.random.uniform(0, 1, 10 * 10)
A = A.reshape((10, 10))
# symmetrisieren fuer reelle Eigenwerte
for i in range(10):
    for k in range(10):
        A[i, k] = A[k, i]

ews = qr_eigenvalues(A, 1e-5)

for ew in ews:
    x = inverse_iteration(A, ew, 1e-5)
    np.testing.assert_allclose(
        ew * x,
        np.dot(
            A,
            x),
        atol=1e-10,
        err_msg="qr oder inverse_iteration haben ein Problem")

print("Alle Eigenvektoren und Eigenwerte passen")
