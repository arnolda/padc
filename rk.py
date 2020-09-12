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


def rk_step(hc, hA, hb, F, yn, tn):
    """Ein einzelner Runge-Kutta-Schritt."""
    k = []
    for i in range(len(hc)):
        k.append(F(tn + hc[i], yn + np.dot(hA[i, :i], k[:i])))
    return yn + np.dot(hb, k)


def rk_explicit(verfahren, F, y0, tmax, h):
    """Das ganze explizite Runge-Kutta-Verfahren."""
    hA = h * verfahren['A']
    hc = h * verfahren['c']
    hb = h * verfahren['b']

    # Startwert und -zeit
    tn = 0.0
    yn = y0.copy()
    # Ergebnisvektor mit Zeit und Punkten
    result = [np.concatenate(((tn,), yn.copy()))]
    while tn < tmax:
        yn = rk_step(hc, hA, hb, F, yn, tn)
        tn += h
        result.append(np.concatenate(((tn,), yn.copy())))
    return np.array(result)
