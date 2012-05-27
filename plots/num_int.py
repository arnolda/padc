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
# Bessel-Naeherung durch num. Diffgleichung
# und Integraldarstellung
#
############################################
from scipy import *
from scipy.interpolate import *
import math
import matplotlib.pyplot as pyplot

a = 0.1
b = 0.9

x = linspace(0, 1, 200)
xab = linspace(a, b, 200)
m = 0.5*(a+b)

def f(x):
    return 0.1*sin(pi*x)**2 + 1.5*(x-0.3)**2

figure = pyplot.figure(figsize=(8,2))

# Trapezregel
############################################

def prepare_graph(pos, title, support):
    est = lagrange(support, f(support))
    graph = figure.add_subplot(pos)
    graph.text(0.05, 0.85, title, ha="left")
    graph.set_yticks(())
    if m in support:
        l  = (a,m,b)
        lt = ("a", "m", "b")
    else:
        l  = (a,b)
        lt = ("a", "b")
    graph.set_xticks(l)
    graph.set_xticklabels(lt)
    graph.xaxis.set_ticks_position("bottom")
    graph.plot(x, f(x), "k-")
    graph.fill_between(xab, est(xab), color="lightblue", facecolor="lightblue")
    graph.plot(xab, est(xab), "b-")
    graph.plot(support, est(support), "k.")
    graph.axis((0, 1, 0, 1))

prepare_graph(141,
              "Trapezregel",
              array((a, b)))
prepare_graph(142,
              "Simpsonregel",
              array((a, m, b)))
prepare_graph(143,
              "Rechteckregel",
              array((a,)))
prepare_graph(144,
              "Mittelpunktsregel",
              array((m,)))
figure.subplots_adjust(left=0.02, right=0.98)

figure.savefig("num_int.pdf")
