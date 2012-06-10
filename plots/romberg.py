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
# Besselfunktion ueber Besselintegral, Richardson-Extrapolation
#
############################################
from scipy import *
from scipy.special import *
from scipy.interpolate import *
import math
import matplotlib.pyplot as pyplot

a = 0
b = 1
# kleinste und groesste zu beruecksichtigende Ordnung
# bei der Extrapolation
kmin = 0
kmax = 5
k = 12
 
def f(tau):
    return exp(-tau**2)

tgt = erf(b)*sqrt(pi)/2

# Romberg-Stuetzpunkte erzeugen
############################################

h = b - a
tf = h*(0.5*f(a) + 0.5*f(b))

tf_list = [ tf ]
h_list  = [ h ]

for i in range(0, k-1):
    print "tf=", tf
    tf = 0.5*tf + 0.5*h*sum(f(a + (arange(0,2**i) + 0.5)*h))
    h *= 0.5
    tf_list.append(tf)
    h_list.append(h)

tf = array(tf_list)
h = array(h_list)

############################################


h_cont = logspace(log(min(h)), log(max(h)), 100, base=exp(1))

print "tgt=", tgt
print "error=", lagrange(h**2, tf)(0) - tgt
print "best input=", tf[kmax] - tgt
est = lagrange(h[kmin:kmax]**2, tf[kmin:kmax])

figure = pyplot.figure(figsize=(6,4))
figure.subplots_adjust(left=0.2, bottom=0.15)

graph = figure.add_subplot(111)
graph.set_xlabel("h")
graph.set_ylabel("$\Delta T_f(h)$")
graph.ticklabel_format(scilimits=(3,3))

graph.plot(h, tgt - tf, "k+", markersize=4)
graph.plot(h[kmin:kmax], tgt - tf[kmin:kmax], "bo", markersize=6)
graph.plot(h_cont, tgt - est(h_cont**2), "r:")
graph.axis((min(h_cont), 0.251, 1e-8, 4e-3))

inset = pyplot.axes([0.3, 0.47, 0.3, 0.4])
inset.set_xscale("log")
#inset.set_xlabel("h")
inset.set_yscale("log")
#inset.set_ylabel("$\Delta T_f(h)$")
inset.set_yticks([1e-8, 1e-6, 1e-4, 1e-2])
inset.tick_params(labelsize='small')
inset.plot(h, tgt - tf, "k+", markersize=4)
inset.plot(h[kmin:kmax], tgt - tf[kmin:kmax], "bo", markersize=6)
inset.plot(h_cont, tgt - est(h_cont**2), "r:")
inset.axis((min(h_cont), max(h), 1e-8, 1e-1))

# graph = figure.add_subplot(122)
# graph.set_xscale("log")
# graph.set_xlabel("h")
# graph.set_yscale("log")
# graph.set_ylabel("$\Delta T_f(h)$")

# graph.plot(h, tgt - tf, "k+", markersize=4)
# graph.plot(h[kmin:kmax], tgt - tf[kmin:kmax], "bo", markersize=6)
# graph.plot(h_cont, tgt - est(h_cont**2), "r:")


figure.savefig("romberg.pdf")
