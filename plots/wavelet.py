# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
from scipy import *
import numpy.random as rand
import matplotlib.pyplot as pyplot

N=200

x = linspace(0, 2*pi, N+1)

rand.seed(1)
y = sin(x) + 0.1*sin(10*x) + rand.normal(0, 0.1, x.shape)

# da wir nur Zweierpotenzen so einfach hinbekommen, ergaenzen wir mit
# Nullen auf die naechste Zweierpotenz. Da alles lokal ist, schadet das
# nix
steps = int(ceil(log(N)/log(2)))
ypadded = zeros(2**steps)
ypadded[:N+1] = y

# Haar-Wavelet-Transformation
#############################

def haar_trafo(data):
    # Ausgabe ist ein Vektor wie folgt:
    # (f,phi_0), (f, psi_00), (f, psi_10), (f, psi_11),
    # (f, psi_20), (f, psi_21), (f, psi_22), ...

    # Daten mit kleinstem Integrationsschritt multiplizieren
    c = data.copy() / len(data)
    # Temporaerer Puffer, um benoetigte Werte nicht zu ueberschreiben
    ctmp = zeros(c.shape)
    width = len(c)/2
    while width >= 1:
        for n in range(width):
            tmp1 = c[2*n]
            tmp2 = c[2*n+1]
            # Detail
            ctmp[width + n] = tmp1 - tmp2
            # Downsampling
            ctmp[n]         = tmp1 + tmp2
        # Puffer zurueckschreiben
        c[:2*width] = ctmp[:2*width]
        width = width / 2
    return c

def inverse_haar_trafo(coeff, lstart, lend):
    # berechnet einen Ausschnitt an Detailstufen.

    # linke Punkte der Intervalle des feinsten Levels
    x = arange(0,2**(lend))/2.0**(lend)
    # Ausgabewerte
    y = zeros(2**lend)
    # phi mitnehmen auf der niedrigsten Stufe
    if lstart == 0: y[0] = coeff[0]
    width = 1
    coeffstart = 1
    for level in range(lend):
        for n in range(width-1, -1, -1):
            tmp = y[n]
            if level >= lstart:
                y[2*n]     = tmp + width*coeff[coeffstart + n]
                y[2*n + 1] = tmp - width*coeff[coeffstart + n]
            else:
                y[2*n]     = tmp
                y[2*n + 1] = tmp
        coeffstart += width
        width = width * 2
    return x, y

coeff = haar_trafo(ypadded)

print "coeff=", coeff

# Ausgabe
#################################

def pretty_inverse(coeff, lstart, lend):
    x, y = inverse_haar_trafo(coeff, lstart, lend)

    rangex = 2*pi/N*len(coeff)

    # Anfang und Ende der Stufen der Treppen
    xx = zeros(2*len(x))
    yy = zeros(2*len(x))
    for n in range(len(x)):
        xx[2*n]     = rangex*x[n]
        yy[2*n]     = y[n]
        if n + 1 < len(x):
            xx[2*n + 1] = rangex*x[n + 1]
        else:
            xx[2*n + 1] = rangex*(2*x[n] - x[n-1])
        yy[2*n + 1] = y[n]

    return xx, yy

figure = pyplot.figure(figsize=(8,4))

graph = figure.add_subplot(221)
graph.plot(x, y, "ko", markersize=1)
graph.axis([0, 2*pi,-1.1,1.1])

graph = figure.add_subplot(222)
xx, yy = pretty_inverse(coeff, 0, 4)
graph.plot(xx, yy, "k")
graph.axis([0, 2*pi,-1.1,1.1])

graph = figure.add_subplot(223)
xx, yy = pretty_inverse(coeff, 4, 6)
graph.plot(xx, yy, "k")
graph.axis([0, 2*pi,-1.1,1.1])

graph = figure.add_subplot(224)
xx, yy = pretty_inverse(coeff, 6, 8)
graph.plot(xx, yy, "k")
graph.axis([0, 2*pi,-1.1,1.1])

figure.savefig("wavelet.pdf")
