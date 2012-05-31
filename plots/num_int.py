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

# Trapezregel
############################################

def prepare_graph(file, title, support, est, xticks):
    figure = pyplot.figure(figsize=(2,2))
    graph = figure.add_subplot(111)

    graph.text(0.05, 0.95, title, size=10, ha="left", va="top")

    graph.set_yticks(())

    l, lt = xticks
    graph.set_xticks(l)
    graph.set_xticklabels(lt)
    graph.xaxis.set_ticks_position("bottom")
    for tick in graph.xaxis.get_major_ticks():
        tick.label1.set_va('baseline')
        tick.label1.set_y(-0.05)

    graph.plot(x, f(x), "k-")

    graph.fill_between(xab, est(xab), color="lightblue", facecolor="lightblue")
    graph.plot(xab, est(xab), "b-")
    graph.plot(support, est(support), "k.")

    graph.axis((0, 1, 0, 1))

    figure.savefig(file)

def prepare_simple_graph(title, support):
    if m in support:
        ticks = ((a,m,b), ("$a$", "$m$", "$b$"))
    else:
        ticks = ((a,b), ("$a$", "$b$"))

    est = lagrange(support, f(support))

    prepare_graph(title.lower() + ".pdf", title, support, est, ticks)


prepare_simple_graph("Trapezregel",
                     array((a, b)))
prepare_simple_graph("Simpsonregel",
                     array((a, m, b)))
prepare_simple_graph("Rechteckregel",
                     array((a,)))
prepare_simple_graph("Mittelpunktsregel",
                     array((m,)))

N=5
support = linspace(a, b, N)
labels = [ "$a$" ] + [ "$x_%d$" % i for i in range(1,N-1) ] + [ "$b$" ]

est = interp1d(support, f(support), "linear")

prepare_graph("trapezzusammen.pdf", "Zusammengesetzte\nTrapezregel",
              support, est, (support, labels))
