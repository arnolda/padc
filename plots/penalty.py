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
# Straffunktionsverfahren
##############################################

from scipy import *
from scipy.linalg import *
from numpy.random import *
import sys
# da liegen die Armijo-Methode, da sie Teil des Skripts ist
sys.path.append("..")
import matplotlib.pyplot as pyplot

from armijo import armijo_steepest_descent

def penalty(f, gradf, g, gradg, x0, sigma_exp=2, tol=0.1):
    x = x0.copy()
    path = [ x.copy()]
    converged = False
    step = 1.0
    while not converged:
        sigma = .1*step**2.0

        # Funktion + aktuelle Straffunktionen
        def q(x):
            res = f(x)
            for lg in g(x):
                res += minimum(0, sigma*lg)**2
            return res
        def gradq(x):
            res = gradf(x)
            for lg, lgg in zip(g(x), gradg(x)):
                res += 2*minimum(0, sigma*lg)*sigma*lgg
            return res

        # steilster Abstieg mit Armijo-Schrittweite
        x = armijo_steepest_descent(q, gradq, x, tol)

        path.append(x.copy())

        if min(g(x)) > -tol:
            converged = True
        step += 1.0

    return array(path).transpose()

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.06, bottom=0.1, right=0.95, top=0.95)

# Beispiel: kuerzester Abstand zu einem Halbkreis
################################################

graph = figure.add_subplot(121)

axis = (-2.5,2.5,-2.5, 1.5)

# 1. NB: im Kreis vom Radius 2 um 0
# g(x) = 2**2 - x^Tx
# 2. NB: nicht oberhalb von 0
# g(x) = - x.y 
def g(x):
    return array((2.0**2 - dot(x, x), -10*x[1]))
def gradg(x):
    return array(((-2.0*x[0], -2.0*x[1]), (0, -10.0)))

# Funktion: qudratischer Abstand zu Punkt p
def f(x, p):
    xx = x - p
    return dot(xx, xx)
def gradf(x, p):
    xx = x - p
    return 2.0*xx

def showpath(p, color, sym):
    path = penalty(lambda x: f(x, p), lambda x: gradf(x, p),
                    g, gradg, array((0.0,-1.0)))
    # Zielpunkt, freies Minimum
    graph.plot(p[0], p[1], sym, markeredgecolor=color,
               markerfacecolor="white", markersize=8,clip_on=False)
    # Punkte und Linie
    graph.plot(path[0], path[1], "-", color=color, linewidth=1)
    graph.plot(path[0], path[1], sym,
               markeredgecolor=color,markerfacecolor=color,
               markersize=4)

# Ausgewaehlte Levelsets, 0-Kante
levelsets = ((1, "k:"), (2, "k--"), (5, "k-"))

for step, style in levelsets:
    sigma = .1*step**2.0

    def contour(niveau, style):
        r = sqrt(4.0 - niveau)
        x, y = [], []
        for phi in linspace(0, 2*pi, 400):
            xx = r*cos(phi)
            yy = -r*sin(phi)
            if yy <= -niveau/10.0:
                x.append(xx)
                y.append(yy)

        graph.plot(x, y, style)
        graph.plot((x[0], x[-1]), (-niveau/10.0, -niveau/10.0), style)

    contour(-0.5/sigma, style)

# im Inneren
showpath(array((-1,-1)), "red", "o")
# oben
showpath(array((-1,1)), "green", "^")
# links oben
showpath(array((2.5,1)), "blue", "v")

graph.axis(axis)

# Penaltyfunktion entlang y-Achse
##################################

graph = figure.add_subplot(122)

def q(x, ff, sigma):
    res = ff*f(x, array((1,1)))
    for lg in g(x):
        res += minimum(0, sigma*lg)**2
    return res

for step, style in levelsets:
    sigma = .1*step**2.0

    x = linspace(-3, 2,400)
    y = [ q(array((0, xx)), 1, sigma) for xx in x]
    graph.plot(x, y, style, color="green",linewidth=1)

    y = [ q(array((0, xx)), 0, sigma) for xx in x]
    graph.plot(x, y, style,linewidth=2)

graph.axis((-3,2,0,10))

figure.savefig("penalty.pdf")
