# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Waermeleitungsgleichung mittels finiter Differenzen und RK4
#############################################################
import numpy as np
import matplotlib.pyplot as plt
import rk
import rk_tableaus

D = 0.5    # Diffusionskonstante
L = 20.0   # Kantenlaenge Simulationsbox
N = 100    # Punkte der Raumdiskretisierung
tmax = 80  # Zeitraum
dt = 0.01  # Zeitschritt


def laplace_operator():
    """Laplace mit 0-Randbedingung."""
    h = L / N
    laplace = np.zeros((N, N))
    for i in range(N):
        if i > 0:
            laplace[i, i - 1] = 1.0 / h**2
        laplace[i, i] = -2.0 / h**2
        if i < N - 1:
            laplace[i, i + 1] = 1.0 / h**2
    return laplace


def init_density():
    """Startdichte: ein Teilchen in der Mitte."""
    p0 = np.zeros(N)
    p0[N // 2] = N / L
    return p0


def f(t, p):
    """p' = f(t, p) = D*Laplace p."""
    return D * np.dot(laplace_operator(), p)


def solve(p0, f, tmax):
    """Integration der DGL bei Startdichte p0"""
    return rk.rk_explicit(rk_tableaus.rk_klassisch, f, p0, tmax, dt)


def plot():
    """Ausgabe."""
    tnpns = solve(init_density(), f, tmax)
    # Umpacken in getrennte Arrays fuer pyplot
    ts = [t for t, *_ in tnpns]
    ps = [p for _, *p in tnpns]
    mass = [sum(p) for p in ps]

    figure = plt.figure(figsize=(8, 6))
    x = 0.5 * L / N + np.linspace(-L / 2.0, L / 2.0, N, endpoint=False)

    # links: Verlauf
    graph = figure.add_subplot(121)
    for step in (int(5.0 / dt), int(20.0 / dt), int(tmax / dt)):
        graph.plot(x, ps[step], label=f"t={ts[step]}")
    graph.legend()

    # rechts: Masse
    graph = figure.add_subplot(122)
    graph.plot(ts, L / N * np.array(mass))

    plt.show()


if __name__ == "__main__":
    plot()
