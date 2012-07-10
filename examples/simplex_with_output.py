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
import sys

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
        cn = c[:, nichtbasis]
        cb = c[:, basis]
        xb = dot(Abinv, b)

        x = zeros(n)
        x[basis] = xb

        print
        print "**** Neuer Simplex-Schritt"
        print "Basis ist", basis
        print "und x = ", x

        # reduzierte Kosten
        r = cn - dot(An.transpose(), dot(Abinv.transpose(), cb))

        print "reduzierte Kosten sind", r

        # beste Abstiegskoordinate s suchen
        minc, spos = minelement(r)
        if minc > -eps:
            # keine Abstiegsrichtung, Minimum gefunden!
            print "x war minimal, das wars!"
            print
            return x
        s = nichtbasis[spos]

        print "Koordinate %d ist beste Abstiegsrichtung" % s

        # rauszuwerfende Variable suchen
        Abinvas = dot(Abinv, A[:,s])
        t = minposelement(xb, Abinvas, eps)

        for i in range(0,len(basis)):
            if Abinvas[i] > eps:
                print "%d-te Basiskoordinate %d erlaubt maximal %f" % \
                    (i, basis[i], x[i] / Abinvas[i])

        if t == None:
            # zulaessige Menge unbegrenzt!
            print "keine Basisgrenze gefunden, zulaessige Menge unbeschraenkt!"
            return

        print "... und %d-te Basiskoordinate %d muss raus" % (t, basis[t])

        # Austausch von j_t und s
        nichtbasis[spos] = basis[t]
        basis[t] = s

        # Update von Abinv
        rank1update(Abinv, t, Abinvas)

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
    
    print "Problem:"
    print "A=", A
    print "b=", b
    print "c=", c
    print

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

    print "Erweitertes Problem:"
    print "A=", A
    print "b=", b
    print "c=", c
    print

    print "****** Phase II in Phase I"
    # Loesung mit Hilfe von Phase 2 suchen
    x = phase2(c, A, b, basis, Abinv, eps)
    print "****** fertig"
    print

    if dot(c, x) > eps:
        print "zulaessige Menge ist leer! Minimum ist", dot(c, x)
        print "Zugehoeriger Punkt ist ", x
        return

    while True:
        # pruefen, ob Basis noch Schattenvariablen enthaelt
        maxb, t = maxelement(basis)
        if maxb <= n: break

        print "Kuenstliche Variable %d muss raus" % maxb

        # echte Ersatzvariable suchen, die nicht in der Basis ist
        for s in range(n):
            if s in basis: continue
            Abinvas = dot(Abinv, A[:,s])
            if Abinvas[t] < 0:
                print "Austausch mit %d-ter Koordinate %d" % (t, basis[t])

                # Ein Tauschpartner!
                basis[t] = s

                print "Neue Basis ist", basis

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

            print "Kein Partner, Zeile %d-%d = %d gestrichen" % (maxb, n, t)
        print

    return basis, Abinv, A[:,:n], b

def simplex(c, A, b, eps=1e-10):
    """
    Simplexverfahren zur Minimierung von c^Tx unter der
    Nebenbedingungen Ax=b und x >= 0. Liefert das Optimum x
    zurueck. eps ist die Toleranz des Algorithmus.
    """

    print "***************"
    print "Phase I"
    print "***************"
    basis, Abinv, A, b = phase1(c, A, b, eps)
    print "***************"
    print "Phase II"
    print "***************"
    return phase2(c, A, b, basis, Abinv, eps)

##############################################
# Verschiedene Testprobleme
##############################################

if len(sys.argv) > 1:
    task = sys.argv[1]
else:
    task = "1"

if task == "1":
    # Abbildung 9.3.1
    # die zweite Bedingung ist offenbar sinnlos
    # und wird in der ersten Phase des Simplex eliminiert
    A = array(((1, 1, 1),
               (2, 2, 2)))
    b = array((1, 2))
    c = (-0.3, 0.4, 0)

elif task == "2":
    # etwas komplizierter, die klassische Kostenoptimierung
    #
    # Bauer hat 40ha die er mit Kartoffeln, Rueben oder Weizen bebauen
    # kann.  Er kann 312 Tage im Jahr arbeiten und 2400 Euro
    # investieren.  Ausserdem moechte der Bauer genauso viele
    # Kartoffeln wie Weizen anzubauen, weil dann der Weizenanbau kaum
    # zusaetzlich kostet.
    #
    # Kartoffeln kosten pro ha  7 Tage Arbeit und  40 Euro,
    #            und bringen 100 Euro pro ha Erloes
    # Rueben     kosten pro ha 12 Tage Arbeit und 120 Euro
    #            und bringen 250 Euro pro ha Erloes
    # Weizen     kostet pro ha  7 Tage Arbeit und  10 Euro
    #            und bringt 120 Euro pro ha Erloes
    #
    # In ha-Gleichungen:
    #    K + R + W <= 40
    #    40 K + 120 R + 10 W <= 2400
    #    7 K + 12 R + 7 W <= 312
    #    K = W
    #
    #    Gewinn ist 100 K + 250 R + 120 W
    #
    A = array(((1,   1,  1, 1,0,0),
               (40,120, 10, 0,1,0),
               ( 7, 12,  7, 0,0,1),
               ( 1,  0, -1, 0,0,0),
               ), dtype=float)
    b = array((40, 2400, 312, 0), dtype=float)
    c = array((-100, -250, -120, 0, 0, 0), dtype=float)

elif task == "3":
    # unbeschraenkte zulaessige Menge und Zielfunktion
    #
    # y  >=   5 + x  <=>  x - y <= -5
    # y  >= -10 + 2x <=> 2x - y <= 10
    # 
    A = array(((1, -1, -1,  0),
               (2, -1,  0, -1)), dtype=float)
    b = array((-5, 10), dtype=float)
    c = array((-1, -1, 0, 0), dtype=float)

elif task == "4":
    # leere zulaessige Menge
    #
    # y >= 1-x       <=>   x + y >= 1
    # y >= 2+x       <=> -x  + y >= 2
    # y <= 3/2 + x/2 <=> -2x + 2y <= 3
    #
    A = array((( 1, 1, -1, 0, 0),
               (-1, 1,  0,-1, 0),
               (-2, 2,  0, 0, 1)), dtype=float)
    b = array((1, 2, 3), dtype=float)
    c = array((1, 1, 0, 0, 0), dtype=float)

print simplex(array(c), A, b)
