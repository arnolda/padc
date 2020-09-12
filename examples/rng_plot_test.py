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
Pseudozufallszahlentests
"""
import argparse
import itertools as it
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Randmine:
    """RAND-Zufallszahlengenerator mit selbstbestimmten Koeffizienten"""
    def __init__(self):
        self.state = 123

    def next(self):
        "Nächste Zahl berechnen und zurückgeben"
        m = 1 << 32
        a = 49
        b = 1975
        self.state = (a * self.state + b) % m
        return float(self.state & ((1 << 31) - 1)) / ((1 << 31) - 1)


class Rand:
    """Standard-RAND-Zufallszahlengenerator"""
    def __init__(self):
        self.state = 123

    def next(self):
        "Nächste Zahl berechnen und zurückgeben"
        m = 1 << 32
        a = 1103515245
        b = 12345
        self.state = (a * self.state + b) % m
        return float(self.state & ((1 << 31) - 1)) / ((1 << 31) - 1)


class Randu:
    """RANDU-Zufallszahlengenerator"""
    def __init__(self):
        self.state = 123

    def next(self):
        "Nächste Zahl berechnen und zurückgeben"
        m = 1 << 31
        a = 65539
        self.state = (self.state * a) % m
        return float(self.state - 1) / float(m - 1)


class Minstd:
    """Minstd-Zufallszahlengenerator"""
    def __init__(self):
        self.state = 123

    def next(self):
        "Nächste Zahl berechnen und zurückgeben"
        m = (1 << 31) - 1
        a = 16807
        self.state = (self.state * a) % m
        return float(self.state - 1) / float(m - 1)


def main():
    "Hauptroutine"
    parser = argparse.ArgumentParser()
    parser.add_argument("--class", dest="klass",
                        help="Zufallszahlengenerator", metavar="Name",
                        type=str,
                        default="randu")
    parser.add_argument("--limit", dest="limit",
                        help=u"Zeitlimit in 0.1s", metavar=u"Limit",
                        type=int)
    args = parser.parse_args()

    rng_name = args.klass.lower()
    if rng_name == "randmine":
        RngType = Randmine
    elif rng_name == "rand":
        RngType = Randmine
    elif rng_name == "randu":
        RngType = Randu
    elif rng_name == "minstd":
        RngType = Minstd
    else:
        raise ValueError("invalid RNG")

    figure = plt.figure(figsize=(4, 4))

    rng = RngType()
    data = [rng.next() for x in range(1000)]

    graph = Axes3D(figure)

    graph.scatter(data[:-2], data[1:-1], data[2:], s=1, marker="o",
                  edgecolors="blue", facecolors="blue")

    plt.show(block=False)
    for _ in range(args.limit) if args.limit else it.count():
        plt.pause(.1)


if __name__ == "__main__":
    main()
