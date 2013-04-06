# -*- coding: utf-8 -*-
# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugänglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
#
# Simulation der Bahn eines Pendels
# mit Hilfe des Velocity-Verlet-Integrators.
#
############################################
import numpy
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyplot
from optparse import OptionParser

parser = OptionParser()

parser.add_option("--l", dest = "l",
                  help = "Länge des Pendels", metavar = "Länge",
                  type = float,
                  default = 1)
parser.add_option("--dt", dest = "dt",
                  help = "Integrationszeitschritt", metavar = "Zeitschritt",
                  type = float,
                  default = 0.1)
parser.add_option("--start_a", dest = "start_a",
                  help = "Anfangsauslenkung", metavar = "Winkel",
                  type = float,
                  default = 0.1)
parser.add_option("--start_da", dest = "start_da",
                  help = "Anfangsgeschwindigkeit", metavar = "Geschwindigkeit",
                  type = float,
                  default = 0)
parser.add_option("--verlet", dest = "method",
                  help = "Einfachen Integrator statt Velocity-Verlet benutzen",
                  action = "store_const", const = "verlet",
                  default = "simple")

options, args = parser.parse_args()

# Laenge des Pendelarms
l=options.l
# Erdbeschleunigung
g = 9.81
# Zeitschritt
dt = options.dt
# Anzahl Schritte zwischen zwei replots
steps = 10
# Dauer (Breite) der zyklischen Graphen
# in Zeiteinheiten
speicher = 5
# Startposition
start_a = options.start_a
# Startwinkelgeschwindigkeit
start_da = options.start_da
# Integrationsmethode, verlet oder simple
method=options.method

# Graphikfenster mit 4 Subgraphen
#################################
abb = pyplot.figure()
# Pendelskizze
graph_pos = abb.add_subplot(221, title="Pendel")
kurve_arm, kurve_gewicht = graph_pos.plot([],[], "b", [],[], "ro")
graph_pos.axis((-l-0.1,l+0.1,-0.1,2*l+0.1))

graph_v   = abb.add_subplot(222, title="Winkelgeschwindigkeit")
kurve_v, = graph_v.plot([],[])
graph_a   = abb.add_subplot(223, title="Winkel")
kurve_a, = graph_a.plot([],[])
graph_E   = abb.add_subplot(224, title="Energie")
kurve_E, = graph_E.plot([],[])

# Initialisierung
#################

# Speicher fuer plots
ring_a = []
ring_v = []
ring_E = []
ring_t = []
# Position
a = start_a
# Winkelgeschwindigkeit
da = start_da
# Zeit
t = 0

# Hilfsroutinen
##############################

def pendelkraft(a):
    "berechnet die Kraft auf das Pendel"
    return -g/l*numpy.sin(a)

def vv_propagation(kraft):
    "propagiert das System einen Zeitschritt mit Hilfe des VV"

    global a, da, t

    da += 0.5*kraft(a)*dt
    a += da*dt
    da += 0.5*kraft(a)*dt

    t += dt

def simple_propagation(kraft):
    "propagiert das System einen Zeitschritt erster Ordnung"

    global a, da, t

    # simpler Integrator
    da += kraft(a)*dt
    a += da*dt

    t += dt

if method == "simple":
    propagation = simple_propagation
else:
    propagation = vv_propagation

# die Hauptroutine, die staendig die Pendelposition
# in der Zeit propagiert
###################################################

def animate():

    # periodisch die Funktion animate aufrufen
    abb.canvas.manager.window.after(100, animate)

    # einige Schritte propagieren
    for i in range(steps):
        global propagation
        propagation(pendelkraft)

        # Aufnahme der Observablen in Ringpuffern
        #########################################
        global ring_a, ring_da, ring_E, ring_t
        # Position
        ring_a.append(a)
        # Geschwindigkeit
        ring_v.append(da)
        # Gesamtenergie
        ring_E.append(0.5*(l*da)**2 + g*(l - l*numpy.cos(a)))
        # Zeit
        ring_t.append(t)

        if ring_t[-1] - ring_t[0] > speicher:
            del ring_t[0], ring_E[0], ring_a[0], ring_v[0]

    # Ausgabe etwas seltener, um Zeit zu sparen
    ###########################################

    # Skizze der Pendelbewegung
    kurve_gewicht.set_data([l*numpy.sin(a)], [l - l*numpy.cos(a)])
    kurve_arm.set_data([0, l*numpy.sin(a)], [l, l - l*numpy.cos(a)])

    # Observablen
    t_min, t_max = ring_t[0], ring_t[-1]

    kurve_a.set_data(ring_t, ring_a)
    graph_a.axis((t_min, t_max, numpy.min(ring_a), numpy.max(ring_a)))

    kurve_v.set_data(ring_t, ring_v)
    graph_v.axis((t_min, t_max, numpy.min(ring_v), numpy.max(ring_v)))
    # die Energie ist praktisch konstant, daher
    # den Bereich etwas vergroessern
    E_min, E_max = numpy.min(ring_E), numpy.max(ring_E)
    E_min -= 0.01
    E_max += 0.01
    kurve_E.set_data(ring_t, ring_E)
    graph_E.axis((t_min, t_max, E_min, E_max))

    abb.canvas.draw()

abb.canvas.manager.window.after_idle(animate)
pyplot.show()
