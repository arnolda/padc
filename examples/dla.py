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
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
import matplotlib.cm as cm

seed(123)

M = 100

# Startradius, waechst mit der Zeit
R = 5

# Zentrum, wo der Keim sitzt
ctr = (M + 1)/2

# 2D-Simulationsgitter, 1 bedeutet besetzt
world = zeros((M, M))
# Der Keim ist schon da
world[ctr,ctr] = 1

def dla():
    global M, R, ctr

    def center_dist(x, y):
        return sqrt((x - ctr)**2 + (y - ctr)**2)

    path = []

    # Startpunkt auf Kreis um das Zentrum
    phi = uniform(0, 2*pi)
    x = int(round(ctr + R*cos(phi)))
    y = int(round(ctr + R*sin(phi) + 0.5))
    path.append((x,y))

    attached = False
    running = True
    while running:
        xalt, yalt = x, y
        move = randint(0,4)
        if   move == 0: x += -1
        elif move == 1: x +=  1
        elif move == 2: y += -1
        elif move == 3: y +=  1
            
        # Position merken
        path.append((x,y))

        if center_dist(x, y) >= 2*R or \
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

pyplot.hsv()
figure = pyplot.figure(figsize=(4,4))
figure.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)

graph = figure.add_subplot(111)

pfad, = graph.plot([],[], "r-")
kreis, = graph.plot([],[], "k:")

# Plot-Routine
##########################################

# Fuer die Farbe
cnt = 0.0

def animate():
    global graph, pfad, kreis, ctr, cnt

    # neuen Punkt generieren
    pts, R, attached = dla()

    # Pfad und Kreis updaten
    x,y = zip(*pts)
    pfad.set_data(x, y)
    kreis.set_data([ctr + R*cos(phi) for phi in linspace(0,2*pi,100)],
                   [ctr + R*sin(phi) for phi in linspace(0,2*pi,100)])

    # und letzten Punkt permanent ausgeben, falls angelagert
    if attached:
        graph.scatter(x[-1], y[-1], s=4, c=cm.hsv(cnt), linewidth=0)
        cnt += 0.1

    # Groesse anpasse
    graph.axis((0,M-1,0,M-1))

    # periodisch die Funktion wieder aufrufen, wenn noch ok
    if R >= 0.4*M:
        print "Radius %d ist zu gross geworden" % R
    else:
        figure.canvas.manager.window.after(200, animate)

    figure.canvas.draw()

figure.canvas.manager.window.after_idle(animate)
pyplot.show()
