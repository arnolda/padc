# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
from scipy import *
from scipy.linalg import *

def qr_eigenwerte(A, tolerance):
    n = A.shape[0]
    I = identity(n)
    Ak = A.copy()
    converged = False
    while not converged:
        shift = Ak[-1,-1]
        Q, R = qr(Ak - shift*I)
        Ak = dot(R, Q) + shift*I

        if max([abs(Ak[i,k]) for i in range(n) \
                             for k in range(i)]) < tolerance:
            converged = True
    return array([Ak[i,i] for i in range(n)])
