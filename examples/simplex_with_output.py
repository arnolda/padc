# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
"""
Beispiele für das Simplexverfahren
"""
import sys
import numpy as np


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
        Abinv[:, i] -= v * u[i]


def minelement(x):
    """
    Sucht die kleinste Komponente von x und deren Index
    """
    s = np.argmin(x)
    return x[s], s


def maxelement(x):
    """
    sucht die groesste Komponente von x und deren Index
    """
    s = np.argmax(x)
    return x[s], s


def minposelement(x, y, eps):
    """
    sucht den Index der kleinsten Komponente von x_i/y_i,
    wobei y_i > 0 sein soll.
    """
    s, minimum = None, None
    for i, x_cur in enumerate(x):
        if y[i] > eps:
            v = x_cur / y[i]
            if s is None or v < minimum:
                # erstes passendes Element
                # oder spaeter ein kleineres
                s, minimum = i, v
    if s >= 0:
        return s
    return None


def phase2(c, A, b, basis, Abinv, eps=1e-10):
    """
    Phase 2 des Simplexverfahrens. c ist der Kostenvektor der
    Zielfunktion c^Tx, die Nebenbedingungen Ax=b, A habe maximalen
    Zeilenrang. basis ist die Menge der Basiskoordinaten der aktuellen
    Ecke, und Abinv die Inverse von A eingeschraenkt auf die
    Basiskoordinaten. Liefert das Optimum x zurueck und passt die
    Basis sowie Abinv an.  eps ist die Toleranz des Algorithmus.
    """
    nichtbasis = [i for i in range(A.shape[1]) if i not in basis]

    while True:
        An = A[:, nichtbasis]
        cn = c[nichtbasis]
        cb = c[basis]
        xb = np.dot(Abinv, b)

        x = np.zeros(A.shape[1])
        x[basis] = xb

        print("**** Neuer Simplex-Schritt")
        print(f"Basis ist {basis}")
        print(f"und x = {x}")

        # reduzierte Kosten
        r = cn - np.dot(An.transpose(), np.dot(Abinv.transpose(), cb))

        print("reduzierte Kosten sind", r)

        # beste Abstiegskoordinate s suchen
        minc, spos = minelement(r)
        if minc > -eps:
            # keine Abstiegsrichtung, Minimum gefunden!
            print("x war minimal, das wars!")
            return x
        s = nichtbasis[spos]

        print(f"Koordinate {s} ist beste Abstiegsrichtung")

        # rauszuwerfende Variable suchen
        Abinvas = np.dot(Abinv, A[:, s])
        t = minposelement(xb, Abinvas, eps)

        for i, val in enumerate(Abinvas):
            if val > eps:
                print(f"{i}-te Basiskoordinate {basis[i]} "
                      f"erlaubt maximal {x[i] / val}")

        if t is None:
            # zulaessige Menge unbegrenzt!
            print("Keine Basisgrenze gefunden, "
                  "zulaessige Menge unbeschraenkt!")
            return None

        print("... und %d-te Basiskoordinate %d muss raus" % (t, basis[t]))

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

    print("Problem:")
    print("A=", A)
    print("b=", b)
    print("c=", c)

    # Problem erweitern, damit wir eine Loesung kennen
    b = b.copy()
    A = np.concatenate((A, np.identity(m)), axis=1)
    c = np.concatenate((np.zeros(n), np.ones(m)))
    # Ax = b positiv machen
    for i in range(m):
        if b[i] < 0:
            b[i] = -b[i]
            A[i, :n] = -A[i, :n]
    # sichere Ecke
    basis = list(range(n, n + m))
    # Inverse
    Abinv = np.identity(m)

    print("Erweitertes Problem:")
    print("A=", A)
    print("b=", b)
    print("c=", c)

    print("****** Phase II in Phase I")
    # Loesung mit Hilfe von Phase 2 suchen
    x = phase2(c, A, b, basis, Abinv, eps)
    print("****** fertig")

    if np.dot(c, x) > eps:
        print("zulaessige Menge ist leer! Minimum ist", np.dot(c, x))
        print("Zugehoeriger Punkt ist ", x)
        return None, None, None, None

    while True:
        # pruefen, ob Basis noch Schattenvariablen enthaelt
        maxb, t = maxelement(basis)
        if maxb <= n:
            break

        print("Kuenstliche Variable %d muss raus" % maxb)

        # echte Ersatzvariable suchen, die nicht in der Basis ist
        for s in range(n):
            if s in basis:
                continue
            Abinvas = np.dot(Abinv, A[:, s])
            if Abinvas[t] < 0:
                print("Austausch mit %d-ter Koordinate %d" % (t, basis[t]))

                # Ein Tauschpartner!
                basis[t] = s

                print("Neue Basis ist", basis)

                rank1update(Abinv, t, Abinvas)
                break
        else:
            # Schleife durchgelaufen, ohne Ersatzvariable zu finden
            # -> Matrix linear abhaengig, q-te Zeile streichen
            del basis[t]
            q = maxb - n
            A = np.delete(A, q, 0)
            b = np.delete(b, q, 0)
            Abinv = np.delete(Abinv, q, 1)
            Abinv = np.delete(Abinv, t, 0)

            print("Kein Partner, Zeile %d-%d = %d gestrichen" % (maxb, n, t))

        print()

    return basis, Abinv, A[:, :n], b


def simplex(c, A, b, eps=1e-10):
    """
    Simplexverfahren zur Minimierung von c^Tx unter der
    Nebenbedingungen Ax=b und x >= 0. Liefert das Optimum x
    zurueck. eps ist die Toleranz des Algorithmus.
    """

    print("***************")
    print("Phase I")
    print("***************")
    basis, Abinv, A, b = phase1(c, A, b, eps)
    if basis is None:
        # keine zulässige Lösung
        return None
    print("***************")
    print("Phase II")
    print("***************")
    return phase2(c, A, b, basis, Abinv, eps)


def main():
    """Verschiedene Testprobleme"""

    if len(sys.argv) > 1:
        task = sys.argv[1]
    else:
        task = "1"

    if task == "1":
        # Abbildung 9.3.1
        # die zweite Bedingung ist offenbar sinnlos
        # und wird in der ersten Phase des Simplex eliminiert
        A = np.array(((1, 1, 1),
                      (2, 2, 2)))
        b = np.array((1, 2))
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
        A = np.array((( 1,   1,  1, 1, 0, 0),
                      (40, 120, 10, 0, 1, 0),
                      ( 7,  12,  7, 0, 0, 1),
                      ( 1,   0, -1, 0, 0, 0)), dtype=float)
        b = np.array((40, 2400, 312, 0), dtype=float)
        c = np.array((-100, -250, -120, 0, 0, 0), dtype=float)

    elif task == "3":
        # unbeschraenkte zulaessige Menge und Zielfunktion
        #
        # y  >=   5 + x  <=>  x - y <= -5
        # y  >= -10 + 2x <=> 2x - y <= 10
        #
        A = np.array(((1, -1, -1,  0),
                      (2, -1,  0, -1)), dtype=float)
        b = np.array((-5, 10), dtype=float)
        c = np.array((-1, -1, 0, 0), dtype=float)

    elif task == "4":
        # leere zulaessige Menge
        #
        # y >= 1-x       <=>   x + y >= 1
        # y >= 2+x       <=> -x  + y >= 2
        # y <= 3/2 + x/2 <=> -2x + 2y <= 3
        #
        A = np.array((( 1, 1, -1,  0, 0),
                      (-1, 1,  0, -1, 0),
                      (-2, 2,  0,  0, 1)), dtype=float)
        b = np.array((1, 2, 3), dtype=float)
        c = np.array((1, 1, 0, 0, 0), dtype=float)

    print(simplex(np.array(c), A, b))


if __name__ == "__main__":
    main()
