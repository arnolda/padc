# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Geschwindigkeitsautokorrelationsfunktion
#
############################################
import math
import matplotlib.pyplot as pyplot
import numpy as np

vac1  = np.loadtxt("v_ac_1.data").transpose()
vac05 = np.loadtxt("v_ac_0.5.data").transpose()

# Abschaetzen von tau durch mutige Integration
def estimate_tau(d):
    t, v = d
    dt = t[1] - t[0]
    I = 0
    for f in v:
        I += f*dt
    return I/v[0]

print(f"vac 1 -> tau = {estimate_tau(vac1)}")
print(f"vac 0.5 -> tau = {estimate_tau(vac05)}")

figure = pyplot.figure(figsize=(4,4))
figure.subplots_adjust(left=0.16)

graph = figure.add_subplot(111)
graph.set_xlabel("t")
graph.set_ylabel("C(v,v)")
graph.plot(vac1[0], vac1[1],   "b-")
graph.plot(vac05[0], vac05[1], "r:")

figure.savefig("v_ac.pdf")
