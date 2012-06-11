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
# Perkolation
##############################################
from scipy import *
from numpy.random import *
import matplotlib.pyplot as pyplot

seed(123)

def try_percolation(M, p):
    # 2D-Simulationsgitter, 1 bedeutet besetzt
    world = zeros((M,M))

    # Besetzung
    ##########################################
    # Liste der besetzen Punkte, fuer die Ausgabe    
    structure = []

    for x in range(M):
        for y in range(M):
            u = uniform(0, 1)
            world[x,y] = (u < p)
            if world[x,y]:
                structure.append((x,y))

    # Flutung von oben nach unten
    ##########################################
    # Liste der unbearbeiteten Punkte
    connected_undone = []
    # ... und der bearbeiteten wie unbearbeiteten
    connected = []

    # Initialisierung von oben
    for x in range(M):
        if world[x, M-1]:
            connected_undone.append((x, M-1))
            connected.append((x, M-1))

    # Durchfluten
    while connected_undone:
        # zu bearbeitendes Element holen und als getan markieren
        cx, cy = connected_undone.pop()
        for (nx, ny) in (cx - 1, cy), (cx + 1, cy), \
                (cx, cy - 1), (cx, cy + 1):
            # 1. Index gueltig (wichtig gerade am Anfang!)
            # 2. besetzt
            # 3. noch nicht gesehen
            if nx >= 0 and nx < M and ny >= 0 and ny < M and \
                    world[nx, ny] and \
                    (nx, ny) not in connected:
                connected.append((nx, ny))
                connected_undone.append((nx, ny))

    # Ausgang unten suchen
    percolating = False
    for x in range(M):
        if (x, 0) in connected:
            percolating = True

    return percolating, structure, connected

# Simulationen
##########################################

ps = linspace(0.5, 0.7, 7)

threshold = 0.6

def sims(M, ps, samples):
    probs = []
    for p in ps:
        successful = 0
        for cnt in range(samples):
            if cnt % 10 == 0: print "Test ", cnt
            percolating, structure, connected = try_percolation(M, p)
            if percolating:
                successful += 1

        probs.append(float(successful)/samples)

    return probs

# Graphik
##########################################

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, bottom=0.15, wspace=0.3, right=0.95, top=0.95)

##########################################

graph = figure.add_subplot(121)

M=50
percolating, structure, connected = try_percolation(M, 0.6)

sx, sy = zip(*structure)
cx, cy = zip(*connected)

graph.plot(sx, sy, "k+", markersize=4)
graph.plot(cx, cy, "bo", markersize=4)
graph.axis((-0.5,M-0.5, -0.5,M-0.5))

##########################################

graph = figure.add_subplot(122)

probs100 = sims(100, ps, 10)
probs20 = sims(20, ps, 100)
probs5 = sims(5, ps, 100)

graph.plot(ps, probs5, "g-.")
graph.plot(ps, probs5, "g*")
graph.plot(ps, probs20, "b--")
graph.plot(ps, probs20, "bo")
graph.plot(ps, probs100, "r:")
graph.plot(ps, probs100, "rD")
graph.xaxis.set_label_text("$p$")
graph.xaxis.set_ticks(linspace(0.4,0.8,5))
graph.yaxis.set_label_text("Perkolationsw-keit")
graph.axis((0.495,0.705,0,1))

figure.savefig("percolation.pdf")
