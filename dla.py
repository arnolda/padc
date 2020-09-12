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
Diffusionslimitierte Aggregation (DLA) in 2D.
"""
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(123)

# Laenge einer Gitterseite
M = 201
# Anzahl der anzulagernden Teilchen
N = 1500


def dla():
    """Anlagerung der Punkte."""
    # Zentrum, wo der Keim sitzt
    ctr = (M + 1) // 2
    # 2D-Simulationsgitter, 1 bedeutet besetzt
    world = np.zeros((M, M))
    # Der Keim ist schon da
    world[ctr, ctr] = 1

    # Anlagerung
    ##########################################

    # angelagerte Punkte fuer die Ausgabe
    pt_history = []
    # Startradius, waechst mit der Zeit
    R = 5

    def center_dist(x, y):
        """Abstand zum Zentrum."""
        return np.sqrt((x - ctr)**2 + (y - ctr)**2)

    for n in range(N):
        print(f"Anfuegen von Teilchen {n}")
        attached = False
        while not attached:
            # Startpunkt auf Kreis um das Zentrum
            phi = np.random.uniform(0, 2 * np.pi)
            x = int(round(ctr + R * np.cos(phi)))
            y = int(round(ctr + R * np.sin(phi) + 0.5))

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

                if center_dist(x, y) >= 2 * R or \
                        x < 0 or x >= M or y < 0 or y >= M:
                    # zu weit draussen, abbrechen und neu starten
                    running = False
                elif world[x, y] == 1:
                    # Angelagert!
                    world[xalt, yalt] = 1
                    # Angelargertes Teilchen merken
                    pt_history.append((xalt, yalt))
                    # Radius anpassen
                    R = max(R, center_dist(xalt, yalt) + 5)
                    if R >= 0.4 * M:
                        raise Exception("Radius bzw. N zu gross!")
                    running = False
                    attached = True
    return ctr, R, pt_history


def plot():
    """Ausgabe des Aggregats."""

    # Liste von Paaren in Paar von Listen verwandeln
    ctr, R, points = dla()
    xx, yy = zip(*points)

    plt.hsv()
    plt.scatter(xx, yy, c=range(len(xx)))
    plt.axis((ctr - R, ctr + R, ctr - R, ctr + R))
    plt.show()


if __name__ == "__main__":
    plot()
