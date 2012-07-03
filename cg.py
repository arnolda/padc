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

def conjugate_gradient(A, b, x0, tol=1e-10):
    x = x0.copy()
    r = b - dot(A, x)
    rr = dot(r,r)
    d = r.copy()
    while sqrt(rr) > tol:
        lmbda = dot(d, r) / dot(d, dot(A, d))
        x += lmbda*d
        r -= lmbda*dot(A,d)
        oldrr = rr
        rr = dot(r,r)
        # orthogonalisieren
        d = r + rr/oldrr * d
    return x
