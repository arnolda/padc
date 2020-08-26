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
# Diffusionslimitierte Aggregation (DLA) in 2D
#
##############################################
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation

np.random.seed(123)

M = 100

# Startradius, waechst mit der Zeit
R = 5

# Zentrum, wo der Keim sitzt
ctr = (M + 1) // 2

# 2D-Simulationsgitter, 1 bedeutet besetzt
world = np.zeros((M, M))
# Der Keim ist schon da
world[ctr, ctr] = 1


def dla():
    global M, R, ctr

    def center_dist(x, y):
        return np.sqrt((x - ctr)**2 + (y - ctr)**2)

    path = []

    # Startpunkt auf Kreis um das Zentrum
    phi = np.random.uniform(0, 2 * np.pi)
    x = int(round(ctr + R * np.cos(phi)))
    y = int(round(ctr + R * np.sin(phi) + 0.5))
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

        if center_dist(x, y) >= 2 * R or \
                x < 0 or x >= M or y < 0 or y >= M:
            # zu weit draussen, abbrechen und neu starten
            running = False
        elif world[x, y] == 1:
            # Angelagert!
            world[xalt, yalt] = 1
            # Radius anpassen
            R = max(R, center_dist(xalt, yalt) + 5)
            running = False
            attached = True

    return path, R, attached

# Ausgabe
##########################################
fig, axis = plt.subplots(1, 1, figsize=(4, 4))

pfad, = axis.plot([], [], "r-")
kreis, = axis.plot([], [], "k:")

# Plot-Routine
##########################################

def animate(*args):
    global axis, pfad, kreis, ctr

    # neuen Punkt generieren
    pts, R, attached = dla()

    # Pfad und Kreis updaten
    x, y = zip(*pts)
    pfad.set_data(x, y)
    kreis.set_data([ctr + R * np.cos(phi) for phi in np.linspace(0, 2 * np.pi, 100)],
                   [ctr + R * np.sin(phi) for phi in np.linspace(0, 2 * np.pi, 100)])

    # und letzten Punkt permanent ausgeben, falls angelagert
    if attached:
        axis.scatter(x[-1], y[-1], s=4, linewidth=0)

    # Groesse anpasse
    axis.axis((0, M - 1, 0, M - 1))

    if R >= 0.4 * M:
        print(f"Radius {R} ist zu gross geworden")
        return

    return pfad, kreis


anim = FuncAnimation(fig, animate, frames=100, interval=20, blit=True, repeat=False)
plt.show(block=False)
plt.pause(3)
plt.close()
