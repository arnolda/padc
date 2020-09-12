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
2d-Poisson mittels Fourier.
"""
import numpy as np
import matplotlib.pyplot as plt


def charge_density(sigma2, L, N):
    """
    Ladungsdichte einer gaussschen Ladung im Zentrum.
    Das Gitter hat N Punkte pro Kante der Laenge L.
    """
    rho = np.zeros((N, N))
    prefactor = 0.5 * L**2 / np.pi / sigma2
    for i in range(N):
        for k in range(N):
            x, y = k * L / N, i * L / N
            d2 = (x - 0.5 * L)**2 + (y - 0.5 * L)**2
            rho[i, k] += prefactor * np.exp(-0.5 * d2 / sigma2)
    # neutralisieren
    rho -= sum(sum(rho)) / N**2
    return rho


def solve(rho, L, N):
    """Loesung per Fouriertransformation."""
    rho_fft = np.fft.fft2(rho)
    # Laplace-Operator (2 pi n / L)^2
    for nx in range(N // 2 + 1):
        for ny in range(N // 2 + 1):
            if nx == 0 and ny == 0:
                rho_fft[nx, ny] = 0
            else:
                n2 = (nx**2 + ny**2) * (2 * np.pi / N)**2
                rho_fft[nx, ny] /= n2
                if 0 < nx < N // 2:
                    rho_fft[N - nx, ny] /= n2
                if 0 < ny < N // 2:
                    rho_fft[nx, N - ny] /= n2
                if 0 < nx < N // 2 and 0 < ny < N // 2:
                    rho_fft[N - nx, N - ny] /= n2
    return np.real(np.fft.ifft2(rho_fft))


def plot():
    """Ausgabe als Heatmap."""
    L, N = 1.0, 50
    rho = charge_density(0.01, L, N)
    psi = solve(rho, L, N)

    im = plt.imshow(psi, interpolation="bilinear", origin="lower",
                    extent=(0, L, 0, L))
    plt.colorbar(im)
    plt.show()


if __name__ == "__main__":
    plot()
