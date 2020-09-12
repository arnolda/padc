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
Einfacher iterativer Poisson-Boltzmann-Loeser.
"""
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

L    = 5.0    # Kantenlaenge des Quadrats
N    = 50     # Punkte der Diskretisierung
lb   = 0.7    # Bjerrumlaenge
cinf = 1e-3   # Salzkonzentration am Rand, 1-mmolar
tol  = 1e-2   # maximales Residuum
h    = L / N  # Schrittweite


def linindex(x, y):
    """Fortran/NumPy-artige linearisierte Indizierung."""
    return y + N * x


def fixed_charge_density():
    """
    Fixes rho und zugaenglichen Bereich fuer eine Ladung im Zentrum
    aufsetzen. Die Ionen koennen nicht in die fixe Ladung eindringen.
    """
    rho_fix, chi = np.zeros(N * N), np.zeros(N * N)
    r = 0.1  # Radius der fixen Ladung im Zentrum
    for i in range(N):
        for k in range(N):
            x, y = k * h, i * h
            d = np.sqrt((x - 0.5 * L)**2 + (y - 0.5 * L)**2)
            if d <= r:
                rho_fix[linindex(i, k)] += 1.0 / np.pi / r**2
            chi[linindex(i, k)] = (rho_fix[linindex(i, k)] == 0.0)
    return rho_fix, chi


def laplace_operator():
    """Diskreten 2d-Laplace-Operator mit neutralem 0-Rand berechnen."""
    operator = np.zeros((N * N, N * N))
    for y in range(N):
        for x in range(N):
            eqn = linindex(x, y)
            operator[eqn, linindex(x, y)] = -4 / h**2
            if x < N - 1:
                operator[eqn, linindex(x + 1, y)] = 1 / h**2
            if x > 0:
                operator[eqn, linindex(x - 1, y)] = 1 / h**2
            if y < N - 1:
                operator[eqn, linindex(x, y + 1)] = 1 / h**2
            if y > 0:
                operator[eqn, linindex(x, y - 1)] = 1 / h**2
    return operator


def solve(rho_fix, chi):
    """Iterativer Loeser."""
    laplace = laplace_operator()

    psi = np.zeros(N * N)  # Potential
    lu = linalg.lu_factor(laplace)

    while True:
        # aktuelle vollstaendige Ladungsdichte
        rho = rho_fix + cinf * 2 * np.sinh(-psi) * chi
        residual = np.dot(laplace, psi) / (4 * np.pi * lb) + rho
        print(f"Residuum ist {np.linalg.norm(residual)}")
        if np.linalg.norm(residual) < tol:
            break
        psi = linalg.lu_solve(lu, -4 * np.pi * lb * rho)

    return psi


def plot(psi, chi):
    """Ausgabe der resultierende positiven Ionendichte."""
    n = (cinf * np.exp(-psi) * chi).reshape((N, N))
    im = plt.imshow(n, origin="lower", extent=(0, L, 0, L))
    plt.colorbar(im)
    plt.show()


if __name__ == "__main__":
    rho_fix, chi = fixed_charge_density()
    psi = solve(rho_fix, chi)
    plot(psi, chi)
