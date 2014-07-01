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
# Waermeleitungsgleichung
##############################################
from scipy import *
from scipy.linalg import *
import matplotlib.pyplot as pyplot
import sys
# da liegt der Integrator, da er Teil des Skripts ist
sys.path.append("..")

# Diffusionskonstante
D = 0.5

# Kantenlaenge Simulationsbox
L = 20.0
# innere Punkte der Raumdiskretisierung
N = 199

# Zeitraum
tmax = 80
# Zeitschritt
dt = 0.01

# Initialisierung
#############################################

h = L/(N+1)

# Raumdiskretisierung, Laplace mit 0-Randbedingung
Laplace = zeros((N,N))

for i in range(N):
    if i > 0:
        Laplace[i, i-1] = 1.0/h**2
    Laplace[i, i] = -2.0/h**2
    if i < N-1:
        Laplace[i, i+1] = 1.0/h**2

# p' = f(t, p)
def fhomogen(t, p):
    return D*dot(Laplace, p)

from rk import rk_explicit, rk_klassisch

def unpack(tnpns):
    # Nettomasse
    mass = []
    ps = []
    ts = []
    for pt in tnpns:
        mass.append(sum(pt[1:]))
        ps.append(pt[1:])
        ts.append(pt[0])

    return ts, ps, mass

# Ausgabe
#############################################

figure = pyplot.figure(figsize=(8,6))
figure.subplots_adjust(left=0.15, right=0.95, top=0.95, wspace=0.3)

x = linspace(-L/2.0, L/2.0, N+2)

# links: homogen, ein Teilchen
#############################################

# Am Anfang ein Teilchen in der Mitte
p0 = zeros(N)
p0[N/2] = 1.0/h

tnpns = rk_explicit(rk_klassisch, fhomogen, p0, tmax, dt)
ts, ps, mass = unpack(tnpns)

# oben: Verlauf
#############################################
graph = figure.add_subplot(221)

for step, style in ((int(5.0/dt), "r-"), (int(20.0/dt), "g--"), (int(tmax/dt), "b:")):
    graph.plot(x, [0] + list(ps[step]) + [0], style, linewidth=2, label=("t=%0.1f" % ts[step]))

graph.axis((-L/2.0,L/2.0,0,0.2))
graph.xaxis.set_label_text("$x$")
graph.yaxis.set_label_text("$p(x,t)$")


# unten, Masse
#############################################
graph = figure.add_subplot(223)

graph.plot(ts, h*array(mass),  "r:", linewidth=2, label="dt=%.3f" % dt)

for ddt, style in ((1.3*dt, "y--"), (1.4*dt, "k-")):
    tnpns2 = rk_explicit(rk_klassisch, fhomogen, p0, tmax, ddt)
    ts2, ps2, mass2 = unpack(tnpns2)
    graph.plot(ts2, h*array(mass2), style, linewidth=2, label="dt=%.3f" % ddt)

graph.axis((0,tmax,0.0,2.0))
graph.xaxis.set_label_text("$t$")
graph.yaxis.set_label_text("$m(t)$")
graph.legend(prop={"size": 10})

# rechts: Quelle in der Mitte, Senke am Rand
#############################################

# p' = f(t, p)
def fhomogen(t, p):
    diff = D*dot(Laplace, p)
    # Delta-Quellen bei L/2 und L/4
    diff[N/2] += 1.0/h
    diff[N/5] += 0.5/h
    return diff

# Am Anfang nix
p0 = zeros(N)

tmax = 500

tnpns = rk_explicit(rk_klassisch, fhomogen, p0, tmax, dt)
ts, ps, mass = unpack(tnpns)

# oben: Verlauf
#############################################
graph = figure.add_subplot(222)

for step, style in ((int(5.0/dt), "r-"), (int(20.0/dt), "g--"), (int(80.0/dt), "b:"), (int(tmax/dt), "k-")):
    graph.plot(x, [0] + list(ps[step]) + [0], style, linewidth=2, label=("t=%0.0f" % ts[step]))

graph.axis((-L/2.0,L/2.0,0,20))
graph.xaxis.set_label_text("$x$")
graph.yaxis.set_label_text("$p(x,t)$")
graph.legend(prop={"size": 10})

# unten, Masse
#############################################
graph = figure.add_subplot(224)

graph.plot(ts, h*array(mass),  "r:", linewidth=2)
graph.xaxis.set_label_text("$t$")
graph.yaxis.set_label_text("$m(t)$")
graph.axis((0,tmax,0,150))

figure.savefig("waermeleitung.pdf")
