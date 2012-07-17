# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Waermeleitungsgleichung mittels finiter Differenzen und RK4
#############################################################
from scipy import *
import matplotlib.pyplot as pyplot
from rk import rk_explicit, rk_klassisch

D = 0.5   # Diffusionskonstante
L = 20.0  # Kantenlaenge Simulationsbox
N = 200   # Punkte der Raumdiskretisierung
tmax = 80 # Zeitraum
dt = 0.01 # Zeitschritt

# Raumdiskretisierung, Laplace mit 0-Randbedingung
h = L/N
Laplace = zeros((N,N))
for i in range(N):
    if i > 0:   Laplace[i, i-1] =  1.0/h**2
    Laplace[i, i]               = -2.0/h**2
    if i < N-1: Laplace[i, i+1] =  1.0/h**2

# Startdichte: ein Teilchen in der Mitte
p0 = zeros(N)
p0[N/2] = 1.0/h

# p' = f(t, p) = D*Laplace p
def f(t, p): return D*dot(Laplace, p)
tnpns = rk_explicit(rk_klassisch, f, p0, tmax, dt)

# Umpacken in getrennte Arrays fuer pyplot
# Nettomasse gleich mit berechnen
ts, ps, mass = [], [], []
for pt in tnpns:
    ts.append(pt[0])
    ps.append(pt[1:])
    mass.append(sum(pt[1:]))

# Ausgabe
#############################################
figure = pyplot.figure(figsize=(8,6))
x = linspace(-L/2.0, L/2.0, N, endpoint=False)

# links: Verlauf
graph = figure.add_subplot(121)
for step in (int(5.0/dt), int(20.0/dt), int(tmax/dt)):
    graph.plot(x, ps[step], label=("t=%0.1f" % ts[step]))
graph.legend()

# rechts: Masse
graph = figure.add_subplot(122)
graph.plot(ts, h*array(mass))

pyplot.show()
