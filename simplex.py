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

def rank1update(Abinv, t, y):
    """
    Update von Ab^{-1}, wenn in Ab Spalte t durch Ab*y ersetzt wird
    """
    v = y.copy()
    pivot = v[t]
    v[t] -= 1
    v /= pivot
    u = Abinv[t]
    for i in range(Abinv.shape[0]):
        Abinv[:,i] -= v*u[i]

def minelement(x):
    """
    sucht die kleinste Komponente von x und deren Index
    """
    s = 0
    minimum = x[0]
    for i in range(1,len(x)):
        if x[i] < minimum:
            s, minimum = i, x[i]
    return minimum, s

def maxelement(x):
    """
    sucht die groesste Komponente von x und deren Index
    """
    s = 0
    maximum = x[0]
    for i in range(1,len(x)):
        if x[i] > maximum:
            s, maximum = i, x[i]
    return maximum, s

def minposelement(x, y, eps):
    """
    sucht den Index der kleinsten Komponente von x_i/y_i,
    wobei y_i > 0 sein soll.
    """
    s = -1
    for i in range(0,len(x)):
        if y[i] > eps:
            v = x[i] / y[i]
            if s == -1 or v < minimum:
                # erstes passendes Element
                # oder spaeter ein kleineres
                s, minimum = i, v
    if s >= 0: return s
    else:      return None
# @\newpage@

def phase2(c, A, b, basis, Abinv, eps=1e-10):
    """
    Phase 2 des Simplexverfahrens. c ist der Kostenvektor der
    Zielfunktion c^Tx, die Nebenbedingungen Ax=b, A habe maximalen
    Zeilenrang. basis ist die Menge der Basiskoordinaten der aktuellen
    Ecke, und Abinv die Inverse von A eingeschraenkt auf die
    Basiskoordinaten. Liefert das Optimum x zurueck und passt die
    Basis sowie Abinv an.  eps ist die Toleranz des Algorithmus.
    """

    # Anzahl Variablen
    n = A.shape[1]
    # und Gleichungen
    m = A.shape[0]

    # Nichtbasis berechnen
    nichtbasis = []
    for i in range(n):
        if not i in basis:
            nichtbasis.append(i)

    while True:
        An = A[:, nichtbasis]
        cn = c[nichtbasis]
        cb = c[basis]
        xb = dot(Abinv, b)

        # reduzierte Kosten
        r = cn - dot(An.transpose(), dot(Abinv.transpose(), cb))
        # beste Abstiegskoordinate s suchen
        minc, spos = minelement(r)
        if minc > -eps:
            # keine Abstiegsrichtung, Minimum gefunden!
            x = zeros(n)
            x[basis] = xb
            return x
        s = nichtbasis[spos]

        # rauszuwerfende Variable suchen
        Abinvas = dot(Abinv, A[:,s])
        t = minposelement(xb, Abinvas, eps)

        if t == None:
            # zulaessige Menge unbegrenzt!
            raise Exception("zulaessige Menge unbeschraenkt!")

        # Austausch von j_t und s
        nichtbasis[spos] = basis[t]
        basis[t] = s

        # Update von Abinv
        rank1update(Abinv, t, Abinvas)
# @\newpage@

def phase1(c, A, b, eps=1e-10):
    """
    Phase 1 des Simplexverfahrens. c ist der Kostenvektor der
    Zielfunktion c^Tx, die Nebenbedingungen Ax=b. Liefert eine
    zulaessige Basis, die zugehoerige Inverse sowie A und b zurueck,
    wobei A und b nun maximalen Zeilenrange haben. eps ist
    die Toleranz des Algorithmus.
    """

    # Anzahl Variablen und Gleichungen
    m, n = A.shape
    
    # Problem erweitern, damit wir eine Loesung kennen
    b = b.copy()
    A = concatenate((A, identity(m)),axis=1)
    c = concatenate((zeros(n), ones(m)))
    # Ax = b positiv machen
    for i in range(m):
        if b[i] < 0:
            b[i] = -b[i]
            A[i,:n] = -A[i,:n]
    # sichere Ecke
    basis = range(n, n+m)
    # Inverse
    Abinv = identity(m)

    # Loesung mit Hilfe von Phase 2 suchen
    x = phase2(c, A, b, basis, Abinv, eps)
    if dot(c, x) > eps:
        raise Exception("zulaessige Menge ist leer!")

    while True:
        # pruefen, ob Basis noch Schattenvariablen enthaelt
        maxb, t = maxelement(basis)
        if maxb <= n: break

        # echte Ersatzvariable suchen, die nicht in der Basis ist
        for s in range(n):
            if s in basis: continue
            Abinvas = dot(Abinv, A[:,s])
            if Abinvas[t] < 0:
                # Ein Tauschpartner!
                basis[t] = s
                rank1update(Abinv, t, Abinvas)
                break
        else:
            # Schleife durchgelaufen, ohne Ersatzvariable zu finden
            # -> Matrix linear abhaengig, q-te Zeile streichen
            del basis[t]
            q = maxb - n
            A = delete(A, q, 0)
            b = delete(b, q, 0)
            Abinv = delete(Abinv, q, 1)
            Abinv = delete(Abinv, t, 0)

    return basis, Abinv, A[:,:n], b
# @\newpage@

def simplex(c, A, b, eps=1e-10):
    """
    Simplexverfahren zur Minimierung von c^Tx unter der
    Nebenbedingungen Ax=b und x >= 0. Liefert das Optimum x
    zurueck. eps ist die Toleranz des Algorithmus.
    """

    basis, Abinv, A, b = phase1(c, A, b, eps)
    return phase2(c, A, b, basis, Abinv, eps)
