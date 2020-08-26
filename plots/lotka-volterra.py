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
#
# Lotka-Volterra-DGL
##############################################
import matplotlib.pyplot as plt
import numpy as np
import sys
# da liegt der Integrator da er Teil des Skripts ist
sys.path.append("..")

# Ein Beutetier vermehrt sich einmal pro Jahr
A =  2.0
# und ein Raeuber frisst es etwa einmal alle 10 Jahre
B =  0.1
# 
# die Raeuber sterben einmal pro Monat, ohne etwas zu fressen
C = 12.0
# aber vermehren sich etwas, wenn sie entsprechend fressen
D = 0.1
#
# diese einfache Funktion ist nicht zeitabhaengig
def f(t, y):
    return np.array((+A*y[0] - B*y[0]*y[1],
                  -C*y[1] + D*y[0]*y[1]))

print("Stationaer:", C/D, A/B)
# Loeser
#############################################

from rk import rk_explicit, euler, rk_klassisch

y0 = np.array((100,1))

# Ausgabe
#############################################

figure = plt.figure(figsize=(8,4))
figure.subplots_adjust(left=0.15, right=0.95,wspace=0.3)

# links oben, Expliziter Euler
#############################################
graph = figure.add_subplot(221)

tnyns_ee = rk_explicit(euler, f, y0, 50, 1./365)

t, beute, raeuber = zip(*tnyns_ee)

graph.plot(t, beute, "b:", linewidth=2)
graph.plot(t, raeuber, "r--", linewidth=2)
graph.axis((0,10,0,300))
graph.yaxis.set_label_text("Populationen")

# links unten, RK
#############################################
graph = figure.add_subplot(223)

tnyns_rk = rk_explicit(rk_klassisch, f, y0, 50, 1./365)

t, beute, raeuber = zip(*tnyns_rk)

graph.plot(t, beute, "b:", linewidth=2)
graph.plot(t, raeuber, "r--", linewidth=2)
graph.axis((0,10,0,300))
graph.xaxis.set_label_text("Zeit in Jahren")
graph.yaxis.set_label_text("Populationen")

# rechts, Phasenraum
#############################################
graph = figure.add_subplot(122)

t, ee_beute, ee_raeuber = zip(*tnyns_ee)
t, rk_beute, rk_raeuber = zip(*tnyns_rk)

graph.plot(ee_beute, ee_raeuber, "b--", linewidth=2)
graph.plot(rk_beute, rk_raeuber, "r-", linewidth=2)
graph.xaxis.set_label_text("Beute")
graph.yaxis.set_label_text(u"Räuber")
graph.axis((0,400,0,500))

figure.savefig("lotka-volterra.pdf")
