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
Simulation der Bahn eines Fadenpendels.
"""
import numpy as np
import matplotlib.pyplot as plt

# Laenge des Pendelarms
l = 1
# Erdbeschleunigung
g = 9.81


def F(a):
    """Kraft, die auf die Kugel wirkt"""
    return -g / l * np.sin(a)


def integrate(methode, a, da, dt, T):
    """
    Pendel von t=0 bis t=T in Zeitschritten der Laenge dt integrieren.
    Die Startposition und - winkelgeschwindigkeit sind a und da.
    """
    t = 0  # Zeit
    while t < T:
        if methode == "simple":
            da += F(a) * dt
            a += da * dt
        elif methode == "velocity-verlet":
            da += 0.5 * F(a) * dt
            a += da * dt
            da += 0.5 * F(a) * dt
        t += dt
        yield (t, a, 0.5 * (l * da)**2 + g * (l - l * np.cos(a)))


def plot():
    """Ausgabe von Winkel und Energie in Graphen."""
    tn, an, En = zip(
        *integrate("velocity-verlet", a=0.1, da=0.0, dt=0.01, T=2))

    ausgabe = plt.figure(figsize=(8, 4))

    loesung = ausgabe.add_subplot(121)
    loesung.set_xlabel("T")
    loesung.set_ylabel("Winkel")
    loesung.plot(tn, an)

    energie = ausgabe.add_subplot(122)
    energie.set_xlabel("Zeit")
    energie.set_ylabel("Energie")
    energie.plot(tn, En)

    plt.show()


if __name__ == "__main__":
    plot()
