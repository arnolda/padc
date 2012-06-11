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
from scipy import *
from numpy.random import *
import matplotlib.pyplot as pyplot

seed(123)

def dla(M, N):
    """DLA-Historie fuerr ein Gitter von M**2 Plaetzen und N Teilchen
    """
    # Zentrum, wo der Keim sitzt
    ctr = (M+1)/2
    # 2D-Simulationsgitter, 1 bedeutet besetzt
    world = zeros((M,M))
    # Der Keim ist schon da
    world[ctr,ctr] = 1
    # angelagerte Punkte fuer die Ausgabe
    pt_history = []

    # Startradius, waechst mit der Zeit
    R = 5

    # Anlagerung
    ##########################################

    def center_dist(x, y):
        return sqrt((x - ctr)**2 + (y - ctr)**2)

    for n in range(N):
        if n % 100 == 0:
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
                        raise Exception("Radius ist zu gross geworden, weniger Teilchen anlagern!")
                    running = False
                    attached = True
    return pt_history, ctr, R

# Ausgabe
##########################################

pyplot.hsv()
figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)

##########################################

graph = figure.add_subplot(121)

pts, ctr, R = dla(101,500)
# Graph auf etwa die Groesse des Clusters bringen
R=R-4
xx, yy = zip(*pts)

graph.scatter(xx, yy, s=2, c=range(len(xx)), linewidth=0)
graph.axis((ctr-R,ctr+R,ctr-R,ctr+R))

##########################################

graph = figure.add_subplot(122)

pts, ctr, R = dla(201,2000)
# Graph auf etwa die Groesse des Clusters bringen
R=R-4
xx, yy = zip(*pts)

graph.scatter(xx, yy, s=1, c=range(len(xx)), linewidth=0)
graph.axis((ctr-R,ctr+R,ctr-R,ctr+R))

figure.savefig("dla.pdf")
