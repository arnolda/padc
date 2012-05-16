# -*- coding: utf-8 -*-
#
# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Simulation der Bahn eines Pendels
# mit Hilfe des Velocity-Verlet-Integrators.
#
############################################
import numpy
import matplotlib.pyplot as pyplot

# Laenge des Pendelarms
l=1
# Erdbeschleunigung
g = 9.81
# Zeitschritt
dt1 = 0.1
dt2 = 0.01
# Zeitraum
T = 2
# Startpositionen
start_a1 = 0.1
start_a2 = 1.5
# Startwinkelgeschwindigkeit
start_da = 0

# Daten erzeugen
#################################

def generate(start_a, start_da, dt, integrator):
    global g, l, T
    def F(a):
        return -g/l*numpy.sin(a)
    def E(a, da):
        return 0.5*(l*da)**2 + g*(l - l*numpy.cos(a))

    # Position
    a = start_a
    # Winkelgeschwindigkeit
    da = start_da
    # Zeit
    t = 0

    # Ausgabe
    tn = []
    an = []
    En = []

    while t < T:
        if integrator == "simple":
            # simpler Integrator
            da += F(a)*dt
            a += da*dt
        else:
            # ein Velocity-Verlet-Schritt
            da += 0.5*F(a)*dt
            a += da*dt
            da += 0.5*F(a)*dt

        t += dt

        tn.append(t)
        an.append(a)
        En.append(E(a, da))

    return tn, an, En

# 4 Graphiken
#################################
omega = numpy.sqrt(g/l)

loesung = pyplot.figure(figsize=(8,4))
energie = pyplot.figure(figsize=(8,4))

loesung_1 = loesung.add_subplot(121)
loesung_1.set_xlabel("T")
loesung_1.set_ylabel(u"α")

t=numpy.arange(0,2,0.01)
loesung_1.plot(t, start_a1*numpy.cos(omega*t), "g", lw=1)

loesung_2 = loesung.add_subplot(122)
loesung_2.set_xlabel("T")

loesung_2.plot(t, start_a2*numpy.cos(omega*t), "g", lw=1)

energie_1 = energie.add_subplot(121)
energie_1.set_xlabel("T")
energie_1.set_ylabel("E")
energie_2 = energie.add_subplot(122)
energie_2.set_xlabel("T")

for start_a, lo, en in (start_a1,loesung_1,energie_1), \
        (start_a2,loesung_2,energie_2):
    for dt, method, sym in (dt1,"simple", "b+"), (dt1,"verlet", "bx"), \
            (dt2,"simple", "--r"):
        tn, an, En = generate(start_a, start_da, dt, method)
        lo.plot(tn, an, sym, ms=6, lw=1.5)
        en.plot(tn, En, sym, ms=6, lw=1.5)

for start_a, lo, en in (start_a1,loesung_1,energie_1), \
        (start_a2,loesung_2,energie_2):
    dt, method, sym = dt2, "verlet", "g--"
    tn, an, En = generate(start_a, start_da, dt, method)
    en.plot(tn, En, sym, ms=5, lw=2)

loesung_1.axis((0,2,-0.1,0.1))
loesung_2.axis((0,2,-1.5,1.5))
loesung.savefig("pendel_loesung.pdf")
energie_1.axis((0,2,0.04,0.06))
energie_2.axis((0,2,7,11))
energie.savefig("pendel_energie.pdf")
