# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
import numpy as np


def qr_eigenvalues(A, tolerance):
    n = A.shape[0]
    I = np.identity(n)
    Ak = A.copy()
    while max([abs(Ak[i, k])
               for i in range(n) for k in range(i)]) >= tolerance:
        shift = Ak[-1, -1]
        Q, R = np.linalg.qr(Ak - shift * I)
        Ak = np.dot(R, Q) + shift * I

    return np.array([Ak[i, i] for i in range(n)])
