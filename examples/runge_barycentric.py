# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
#
"""
Baryzentrische Polynom-Interpolation der Rungefunktion
"""
import argparse
import itertools as it
import numpy as np
import scipy.interpolate as ip
import matplotlib.pyplot as plt


def runge(x):
    "Die zu interpolierende Funktion"
    return 1 / (1 + x**2)


class Barycentric:
    """Baryzentrische Polynom-Interpolation"""
    def __init__(self, steps, function):
        self.steps = 7
        self.x = np.linspace(-5, 5, steps)
        self.y = function(self.x)

        self.omegas = np.ones_like(self.x)
        for i in range(steps):
            for k in range(steps):
                if i != k:
                    self.omegas[i] /= (self.x[i] - self.x[k])

    def mu(self, i, x):
        "i-tes Stützpolynom"
        return self.omegas[i] / (x - self.x[i])

    def get(self, x):
        "Interpoliertes Polynom an Stelle x auswerten"
        enum = 0
        denom = 0
        for i in range(self.steps):
            enum += self.y[i] * self.mu(i, x)
            denom += self.mu(i, x)

        return enum / denom


def main():
    """Ausgabe"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", dest="limit",
                        help=u"Zeitlimit in 0.1s", metavar=u"Limit",
                        type=int)
    args = parser.parse_args()

    bary = Barycentric(7, runge)

    pxs = np.logspace(-10, -5, 200)

    figure = plt.figure(figsize=(4, 4))

    graph = figure.add_subplot(111)
    graph.set_xscale("log")
    graph.set_yscale("log")
    graph.plot(pxs, abs(runge(pxs) - bary.get(pxs)), "r")
    graph.plot(pxs, abs(runge(pxs) - ip.lagrange(bary.x, bary.y)(pxs)), "k")

    plt.show(block=False)
    for _ in range(args.limit) if args.limit else it.count():
        plt.pause(.1)


if __name__ == "__main__":
    main()
