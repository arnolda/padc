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


def sor_step(A, b, omega, x):
    n = A.shape[0]
    for j in range(n):
        s = b[j]
        # diese Werte sind schon neu
        s -= np.dot(A[j, :j], x[ :j])
        # und diese Werte sind noch alt
        s -= np.dot(A[j, j + 1:], x[j + 1:])
        # Ueberschreiben von x[j]!
        x[j] = omega * s / A[j, j] + (1 - omega) * x[j]
