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
Diffusionslimitierte Aggregation (DLA) in 2D
"""
import argparse
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

np.random.seed(123)


class DLA:
    """Die eigentliche DLA:

    Arguments:
        M: Größe der Welt
        R: Startradius, waechst mit der Zeit
    """
    def __init__(self, M=100, R=5):
        self.M = M
        self.R = R
        # Zentrum, wo der Keim sitzt
        self.ctr = (M + 1) // 2
        # 2D-Simulationsgitter, 1 bedeutet besetzt
        self.world = np.zeros((M, M))
        # Der Keim ist schon da
        self.world[self.ctr, self.ctr] = 1

    def center_dist(self, x, y):
        "2-Distanz zum Zentrum"
        return np.sqrt((x - self.ctr)**2 + (y - self.ctr)**2)

    def step(self):
        """Ein Teilchen anlagern.

        Returns:
           Ein Tupel aus Pfad und ob angelagert wurde oder nicht.
        """

        path = []

        # Startpunkt auf Kreis um das Zentrum
        phi = np.random.uniform(0, 2 * np.pi)
        x = int(round(self.ctr + self.R * np.cos(phi)))
        y = int(round(self.ctr + self.R * np.sin(phi) + 0.5))
        path.append((x, y))

        attached = False
        running = True
        while running:
            xalt, yalt = x, y
            move = np.random.randint(0, 4)
            if move == 0:
                x += -1
            elif move == 1:
                x += 1
            elif move == 2:
                y += -1
            elif move == 3:
                y += 1

            # Position merken
            path.append((x, y))

            if self.center_dist(x, y) >= 2 * self.R or \
               not 0 <= x < self.M or not 0 <= y < self.M:
                # zu weit draussen, abbrechen und neu starten
                running = False
            elif self.world[x, y] == 1:
                # Angelagert!
                self.world[xalt, yalt] = 1
                # Radius anpassen
                self.R = max(self.R, self.center_dist(xalt, yalt) + 5)
                running = False
                attached = True

        return path, attached


def animate(_frame, dla, axis, cur_path, circle):
    """Plotroutine"""

    # neuen Punkt generieren
    path_data, attached = dla.step()

    # Pfad und Kreis updaten
    x, y = zip(*path_data)
    cur_path.set_data(x, y)
    phi = np.linspace(0, 2 * np.pi, 100)
    circle.set_data(
        [dla.ctr + dla.R * np.cos(phi)],
        [dla.ctr + dla.R * np.sin(phi)])

    # und letzten Punkt permanent ausgeben, falls angelagert
    if attached:
        axis.scatter(x[-1], y[-1], s=4, linewidth=0)

    # Groesse anpasse
    axis.axis((0, dla.M - 1, 0, dla.M - 1))

    if dla.R >= 0.4 * dla.M:
        print(f"Radius {dla.R} ist zu gross geworden")
        return


def main():
    """Hauptroutine"""
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", dest="limit",
                        help=u"Zeitlimit in 0.1s", metavar=u"Limit",
                        type=int)
    args = parser.parse_args()

    fig, axis = plt.subplots(1, 1, figsize=(4, 4))

    dla = DLA(100, 5)

    cur_path, = axis.plot([], [], "r-")
    circle, = axis.plot([], [], "k:")

    anim = FuncAnimation(
        fig, animate, fargs=[dla, axis, cur_path, circle],
        frames=None, interval=20, repeat=False)
    plt.show(block=False)
    for _ in range(args.limit) if args.limit else it.count():
        plt.pause(.1)
    return anim


if __name__ == "__main__":
    main()
