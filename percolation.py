# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Perkolation
##############################################
from scipy import *
from numpy.random import *
import matplotlib.pyplot as pyplot

seed(123)

# Laenge einer Gitterseite
M = 50
# Wahrscheinlichkeit einer Gitterzelle, besetzt zu sein
p = 0.55
# wie oft wir messen, um die W-keit zu bestimmen
samples = 100

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

# Ausgabe
##########################################

successful = 0
for cnt in range(samples):
    print "Test ", cnt
    percolating, structure, connected = try_percolation(M, p)
    if percolating:
        successful += 1
        # Erfolg merken fuer die Ausgabe
        succ_structure = structure
        succ_connected = connected

print "Perkolations-W-keit", float(successful)/samples

if successful > 0:
    # Liste von Paaren in Paar von Listen verwandeln
    xs, ys = zip(*succ_structure)
    xc, yc = zip(*succ_connected)

    pyplot.plot(xs, ys, "k+")
    pyplot.plot(xc, yc, "bo")

    pyplot.show()
