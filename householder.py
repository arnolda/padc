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

def householder(a):
    "Householder-QR-Zerlegung fuer eine n x m-Matrix a"

    def multiplyleft(a, v):
        "a mit Householder-Spiegelung S_v von links multiplizieren"
        v_H_a = dot(v.transpose().conj(), a)
        for k in range(a.shape[1]):
            a[:,k] -= 2*v_H_a[k] * v

    def multiplyright(a, v):
        "a mit Householder-Spiegelung S_v von rechts multiplizieren"
        a_v = dot(a, v)
        for k in range(a.shape[1]):
            a[:,k] -= 2*v[k].conj() * a_v

    r = a.copy()
    q = identity(r.shape[0])
    # Schleife ueber alle Spalten von a bzw. r,
    # die wenigstens ein Subdiagonalelement haben
    for k in range(min(r.shape) - 1):
        # Kopie von a_k, um Spiegelvektor zu bauen
        v = r[:,k].copy()
        # schon bearbeitete Komponenten auf Null
        v[:k] = 0
        # Spiegelvektor @$a + \text{sgn}(a_1)\norm{a} e_1$@ berechnen
        v[k] += sign(v[k])*norm(v)
        # normieren
        v = v / norm(v)

        # Matrix r updaten
        multiplyleft(r, v)
        # und unitaeren Teil q
        multiplyright(q, v)

    return q, r
