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


def conjugate_gradient(A, b, x0, tol=1e-10):
    x = x0.copy()
    r = b - np.dot(A, x)
    rr = np.dot(r, r)
    d = r.copy()
    while np.sqrt(rr) > tol:
        lmbda = np.dot(d, r) / np.dot(d, np.dot(A, d))
        x += lmbda * d
        r -= lmbda * np.dot(A, d)
        oldrr = rr
        rr = np.dot(r, r)
        # orthogonalisieren
        d = r + rr / oldrr * d
    return x
