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

def inverse_iteration(A, l, tolerance):
    n = A.shape[0]
    x = ones(n)
    Ashift = A - l*identity(n)
    converged = False
    while not converged:
        x = solve(Ashift, x)
        x = x / norm(x)
        if norm(l*x - dot(A, x)) < tolerance:
            converged = True
    return x
