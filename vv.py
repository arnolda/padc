# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
from numpy import *

def velocity_verlet(acc, x0, v0, tmax, h):
    def step(acc, tn, xn, vn, h):
        vn += 0.5*h*acc(tn, xn)
        xn += h*vn
        vn += 0.5*h*acc(tn, xn)

    # Startwert und -zeit
    tn = 0.0
    xn = x0.copy()
    vn = v0.copy()
    # Ergebnisvektor
    result = [ concatenate(((tn,), xn.copy(), vn.copy())) ]
    while tn < tmax:
        step(acc, tn, xn, vn, h)
        tn += h
        result.append(concatenate(((tn,), xn.copy(), vn.copy())))
    return array(result)
