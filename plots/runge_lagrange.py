# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Polynom-Interpolation der Rungefunktion
#
############################################
import numpy
import matplotlib.pyplot as pyplot

# die zu interpolierende Funktion
def runge(x): return 1/(1+x**2)

# Neville-Aitken fuers Polynom
def neville(x, y):
    n = len(x)
    gamma = y.copy()
    # speichert gamma_{i-k+1,k} gemaess Schema
    # die ersten k Eintraege sind die alten, also
    # gamma[0] = gamma_{0,0} = gamma_0
    for k in range(1, n):
        for i in range(n-k-1, -1, -1):
            gamma[i+k] = (gamma[i+k] - gamma[i+k-1])/(x[i+k] - x[i])
    return gamma

# Hornerschema zum Auswerten
def horner(x0, x, coeff):
    r = 0
    for k in range(len(x)-1, -1, -1):
        r = r*(x0-x[k]) + coeff[k];
    return r

# interpolierende Polynome erzeugen
coeff_a = []
coeff_c = []
support_a = []
support_c = []

for steps in [3, 5, 7, 11]:
    # aequidistant
    x = numpy.linspace(-5,5,steps)
    y = runge(x)

    support_a.append(x)
    coeff_a.append(neville(x, y))

    # Chebyshev
    x = numpy.zeros(steps)
    for i in range(steps):
        x[i] = 5*numpy.cos((2.*i + 1)/(2.*steps)*numpy.pi)
    y = runge(x)

    support_c.append(x)
    coeff_c.append(neville(x, y))

allsupp_a = []
for l in support_a: allsupp_a += list(l)
allsupp_a = numpy.array(list(set(allsupp_a)))
allsupp_c = []
for l in support_c: allsupp_c += list(l)
allsupp_c = numpy.array(list(set(allsupp_c)))

# Ausgabe
#################################
px = numpy.linspace(-8, 8, 200)

figure = pyplot.figure(figsize=(8,4))

# lagrange
graph = figure.add_subplot(121)
graph.plot(px, runge(px), "black")
graph.plot(px, horner(px, support_a[0], coeff_a[0]), "red")
graph.plot(px, horner(px, support_a[1], coeff_a[1]), "purple")
graph.plot(px, horner(px, support_a[2], coeff_a[2]), "blue")
graph.plot(px, horner(px, support_a[3], coeff_a[3]), "green")
graph.plot(allsupp_a, runge(allsupp_a), "ko", markersize=3)
graph.axis([-5.1,5.1,-0.5,1.5])

# Chebyshev
graph = figure.add_subplot(122)
graph.plot(px, runge(px), "black")
graph.plot(px, horner(px, support_c[0], coeff_c[0]), "red")
graph.plot(px, horner(px, support_c[1], coeff_c[1]), "purple")
graph.plot(px, horner(px, support_c[2], coeff_c[2]), "blue")
graph.plot(px, horner(px, support_c[3], coeff_c[3]), "green")
graph.plot(allsupp_c, runge(allsupp_c), "ko", markersize=3)
graph.axis([-5.1,5.1,-0.5,1.5])

figure.savefig("runge_lagrange.pdf")
