# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Diffusionslimitierte Aggregation (DLA) in 2D
##############################################
from scipy import *
from numpy.random import *
import matplotlib.pyplot as pyplot

seed(123)

# Laenge einer Gitterseite
M = 201
# Anzahl der anzulagernden Teilchen
N = 1500

# Zentrum, wo der Keim sitzt
ctr = (M+1)/2
# 2D-Simulationsgitter, 1 bedeutet besetzt
world = zeros((M,M))
# Der Keim ist schon da
world[ctr,ctr] = 1

# Anlagerung
##########################################

# angelagerte Punkte fuer die Ausgabe
pt_history = []
# Startradius, waechst mit der Zeit
R = 5

def center_dist(x, y):
    return sqrt((x - ctr)**2 + (y - ctr)**2)

for n in range(N):
    print "Anfuegen von Teilchen", n
    attached = False
    while not attached:
        # Startpunkt auf Kreis um das Zentrum
        phi = uniform(0, 2*pi)
        x = int(round(ctr + R*cos(phi)))
        y = int(round(ctr + R*sin(phi) + 0.5))

        running = True
        while running:
            xalt, yalt = x, y
            move = randint(0,4)
            if   move == 0: x += -1
            elif move == 1: x +=  1
            elif move == 2: y += -1
            elif move == 3: y +=  1

            if center_dist(x, y) >= 2*R or \
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
                if R >= 0.4*M:
                    raise Exception("Radius bzw. N zu gross!")
                running = False
                attached = True

# Ausgabe
##########################################

# Liste von Paaren in Paar von Listen verwandeln
xx, yy = zip(*pt_history)

pyplot.hsv()
pyplot.scatter(xx, yy, c=range(len(xx)))
pyplot.axis((ctr-R,ctr+R,ctr-R,ctr+R))
pyplot.show()
