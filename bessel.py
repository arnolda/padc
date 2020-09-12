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
Bessel-Naeherung durch Loesen der diskretisierten DGL.
"""
import numpy as np
from scipy import special
import matplotlib.pyplot as plt

# Ordnung der Besselfunktion
nu = 0
# Rechter Rand
xmax = 15.0
# Anzahl Punkte Differential
N = 31
# Schrittweite Differential
h = xmax / (N - 1)
# Stuetzpunkte
x = np.arange(0, N) * h


def stencil3():
    """3-Punkt-Stencil, Dirichlet-Randbedingung."""
    b = np.zeros(N)
    A = np.zeros((N, N))

    # Startwerte an den Raendern
    A[0, 0] = 1
    b[0] = special.jn(nu, 0)
    A[1, N - 1] = 1
    b[1] = special.jn(nu, (N - 1) * h)

    # N-2 innere Punkte von 1 bis N-2 hinter den Startwerten
    for n in range(1, N - 1):
        #                 x^2 d^2f/dx^2   + x df/dx    + (x^2 - nu^2)f
        A[n + 1, n - 1] =      n**2       - 0.5 * n
        A[n + 1, n]     = -2 * n**2                    + n**2 * h**2 - nu**2
        A[n + 1, n + 1] =      n**2       + 0.5 * n

    return np.linalg.solve(A, b)


def stencil3l():
    """3-Punkt-Stencil, nur linker Rand."""
    b = np.zeros(N)
    A = np.zeros((N, N))

    # Startwert am Rand
    A[0, 0] = 1
    b[0] = special.jn(nu, 0)

    # Ableitung vorgeben
    A[1, 0] = -1.0 / h
    A[1, 2] =  1.0 / h
    # Formel fuer die Ableitung der Besselfunktionen
    if nu == 0:
        b[1] = -special.jn(1, 0)
    else:
        b[1] = 0.5 * (special.jn(nu - 1, 0) - special.jn(nu + 1, 0))

    # N-2 innere Punkte von 1 bis N-2
    for n in range(1, N - 1):
        #                 x^2 d^2f/dx^2 + x df/dx    + (x^2 - nu^2)f
        A[n + 1, n - 1] =      n**2     - 0.5 * n
        A[n + 1, n]     = -2 * n**2                  + n**2 * h**2 - nu**2
        A[n + 1, n + 1] =      n**2     + 0.5 * n

    return np.linalg.solve(A, b)


def stencil5():
    """5-Punkt-Stencil, Dirichlet-Randbedingung."""
    b = np.zeros(N)
    A = np.zeros((N, N))

    # Startwerte an den Raendern, wie zuvor
    A[0, 0] = 1
    b[0] = special.jn(nu, 0)
    A[1, N - 1] = 1
    b[1] = special.jn(nu, (N - 1) * h)

    # innere naechste 2 Punkte durch einfache Naeherung
    for z, n in ((2, 1), (3, N - 2)):
        #             x^2 d^2f/dx^2  + x df/dx    + (x^2 - nu^2)f
        A[z, n - 1] =      n**2      - 0.5 * n
        A[z, n]     = -2 * n**2                   + n**2 * h**2 - nu**2
        A[z, n + 1] =      n**2      + 0.5 * n

    # N-4 innere Punkte von 2 bis N-3 hinter den Startwerten, alle auf Null
    for n in range(2, N - 2):
        z = n + 2
        #              x^2 d^2f/dx^2     + x df/dx     + (x^2 - nu^2)f
        A[z, n - 2] = -1. / 12 * n**2    + 1. / 12 * n
        A[z, n - 1] =  4. / 3  * n**2    - 2. / 3  * n
        A[z, n]     = -5. / 2  * n**2                  + n**2 * h**2 - nu**2
        A[z, n + 1] =  4. / 3  * n**2    + 2. / 3  * n
        A[z, n + 2] = -1. / 12 * n**2    - 1. / 12 * n

    return np.linalg.solve(A, b)


def plot():
    """Ausgabe als Graphen."""

    xfine = np.linspace(0, xmax, 200)

    plt.plot(xfine, special.jn(nu, xfine), "k-", linewidth=0.5)
    plt.plot(x, stencil3(), "gD")
    plt.plot(x, stencil3l(), "r^")
    plt.plot(x, stencil5(), "bo")
    plt.show()


if __name__ == "__main__":
    plot()
