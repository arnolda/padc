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


def multiplyleft(a, i, j, c, s):
    """Matrix a mit Givens-Rotation von links multiplizieren."""
    # eigentlich nur zwei Zeilen verdrehen
    ai_gedreht = +c * a[i, :] + s * a[j, :]
    aj_gedreht = -s * a[i, :] + c * a[j, :]
    a[i, :] = ai_gedreht
    a[j, :] = aj_gedreht


def multiplyright(a, i, j, c, s):
    """Matrix a mit Givens-Rotation von rechts multiplizieren."""
    # eigentlich nur zwei Spalten verdrehen
    ai_gedreht = +c * a[:, i] + s * a[:, j]
    aj_gedreht = -s * a[:, i] + c * a[:, j]
    a[:, i] = ai_gedreht
    a[:, j] = aj_gedreht


def calccs(x, y):
    """Stabilisierte Bestimmung von Sinus und Cosinus fuer die Rotation."""
    if abs(x) >= abs(y):
        r = np.sqrt(1 + (y / x)**2)
        c = np.sign(x) / r
        s = y / (abs(x) * r)
    else:
        r = np.sqrt(1 + (x / y)**2)
        c = x / (abs(y) * r)
        s = np.sign(y) / r
    return c, s


def givens(a):
    """Givens-QR-Zerlegung fuer eine n x m-Matrix a."""
    r = a.copy()
    q = np.identity(r.shape[0])
    # Doppelschleife ueber alle Subdiagonalelemente von a bzw. r
    for k in range(min(r.shape[0] - 1, r.shape[1])):
        for i in range(r.shape[0] - 1, k, -1):
            c, s = calccs(r[i - 1, k], r[i, k])
            # Matrix r updaten
            multiplyleft(r, i - 1, i, c, s)
            # und unitaeren Teil q, hermitesch durch Vz-wechsel im Sinus
            multiplyright(q, i - 1, i, c, s)
    return q, r
