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
# Waermeleitungsgleichung in 1d mittels FEM und FDM
###################################################
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
N = 9
# Schrittweite
h = L/(N+1)

# rechte Seite: Delta-Quellen bei L/2 und L/5
f_delta = zeros(N)
f_delta[N/2] = 1.0/h
f_delta[N/5] = 0.5/h

# FEM
#############################################

def integral_v(p, q):
    if p == q:          return 2./3*h
    elif abs(p-q) == 1: return 1./6*h
    else:               return 0

def integral_dv(p, q):
    if p == q:          return -2./h
    elif abs(p-q) == 1: return  1./h
    else:               return 0

lhs = array([[integral_dv(p, q) for q in xrange(N)] for p in xrange(N)])
# das h korrigiert, das wir hier echte delta-Quellen rechnen
rhs = h*f_delta

pfem_delta = [0] + list(solve(-D*lhs, rhs)) + [0]

rhs = array([sum([f_delta[p]*integral_v(p, q) for p in xrange(N)]) for q in xrange(N)])

print "folded rhs"
print rhs

pfem_hut = [0] + list(solve(-D*lhs, rhs)) + [0]

# Ausgabe
#############################################

figure = pyplot.figure(figsize=(8,3))
figure.subplots_adjust(left=0.1, right=0.95, top=0.95, hspace=0.3, wspace=0.3)

x = linspace(-L/2.0, L/2.0, N+2)

# links: Loesungen
#############################################

graph = figure.add_subplot(121)

graph.plot(x, pfem_delta, "b-")
graph.plot(x, pfem_hut, "r--")
graph.plot(x, pfem_delta, "bD")
graph.plot(x, pfem_hut, "ro")

graph.axis((-L/2.0,L/2.0,0,20))
graph.xaxis.set_label_text("$x$")
graph.yaxis.set_label_text("$p(x)$")

# rechts: finite Elemente
#############################################

graph = figure.add_subplot(222)

form = (0, 1, 0)
pos = array((-h, 0, h))
graph.plot(-h + pos, form, "b--", linewidth=2)
graph.plot( h + pos, form, "r--", linewidth=2)
graph.plot(     pos, form, "k-", linewidth=2)

graph.axis((-5,5,0,1))
graph.xaxis.set_label_text("$x$")
graph.yaxis.set_label_text("$v_p(x)$")

graph = figure.add_subplot(224)

form = (0, 0,-1./h, -1./h, 1./h, 1./h, 0, 0)
pos = array([-10,-h, -h, 0, 0, h, h, 10])
graph.plot(-h + pos, form, "b--", linewidth=2)
graph.plot( h + pos, form, "r--", linewidth=2)
graph.plot(     pos, form, "k-", linewidth=2)

graph.axis((-5,5,-1,1))
graph.xaxis.set_label_text("$x$")
graph.yaxis.set_label_text("$v'_p(x)$")

figure.savefig("fem.pdf")
