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


def inverse_iteration(A, l, tolerance):
    n = A.shape[0]
    x = np.ones(n)
    Ashift = A - l * np.identity(n)
    while np.linalg.norm(l * x - np.dot(A, x)) >= tolerance:
        x = np.linalg.solve(Ashift, x)
        x = x / np.linalg.norm(x)
    return x
