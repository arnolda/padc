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
# Besselfunktion ueber Besselintegral
#
############################################
import math
import matplotlib.pyplot as pyplot
import numpy as np
import scipy.interpolate
import scipy.special

nu = 0

xmax = 15.0

# Integral
# j(nu,x) = int_0^pi cos(nu*tau - x sin(tau)) dtau
############################################

def trapez(x, N):
    """Zusammengesetzte Trapezregel fuer das Besselintegral
    mit N Stuetzstellen auf [0,pi]"""

    h = np.pi/N

    res = np.zeros(x.shape)

    for n in range(len(x)):
        xn = x[n]

        # Trapezregel
        # linker Rand tau = 0
        # und rechter Rand tau = pi
        I = 0.5 + 0.5*np.cos(nu*np.pi)
        # Rest
        tau = np.arange(1, N)*h
        I += sum(np.cos(nu*tau - xn*np.sin(tau)))

        res[n] = h*I/np.pi

    return res

def simpson(x, N):
    """Zusammengesetzte Simpsonregel fuer das Besselintegral
    mit N Stuetzstellen auf [0,pi]"""

    h = np.pi/N

    res = np.zeros(x.shape)

    for n in range(len(x)):
        xn = x[n]

        # Trapezregel
        # linker Rand tau = 0, tau = h
        # und rechter Rand tau = pi
        I = 1 + 4*np.cos(nu*h - xn*np.sin(h)) + np.cos(nu*np.pi)
        # Rest, doppelte Schrittweite
        tau = np.arange(1, N/2)*2*h
        I += sum(  2*np.cos(nu*tau       - xn*np.sin(tau))
                 + 4*np.cos(nu*(tau + h) - xn*np.sin(tau + h)))

        res[n] = h/3*I/np.pi

    return res

############################################

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, right=0.95, wspace=.2)

# Links, Integrale
############################################

x = np.linspace(0, xmax, 200)

graph = figure.add_subplot(121)

graph.plot(x, scipy.special.jn(nu, x), "k-",linewidth=0.5)
xint = np.arange(0,xmax, 1)
graph.plot(xint, trapez(xint, 6), "r.")
xint = np.arange(0.25,xmax, 1)
graph.plot(xint, trapez(xint, 20), "b+", markersize=3)
xint = np.arange(0.5,xmax, 1)
graph.plot(xint, simpson(xint, 6), "gx", markersize=3)
xint = np.arange(0.75,xmax, 1)
graph.plot(xint, simpson(xint, 20), "yD", markeredgewidth=0, markersize=3)
graph.axis((0, xmax, -0.5, 1))

# Rechts, Interpolierende
############################################

N = 6
xn = 12

tau = np.linspace(0, np.pi, 200)

graph = figure.add_subplot(122)

def f(tau):
    return np.cos(nu*tau - xn*np.sin(tau))

support = np.linspace(0,np.pi, N+1)

trapez = scipy.interpolate.interp1d(support, f(support), "linear")
    
graph.plot(tau, f(tau), "k-", linewidth=0.5)
graph.plot(tau, trapez(tau), "r:")
# Simpson-Parabeln
h = np.pi/N
for n in range(N//2):
    lsupp = np.array((2*n*h, (2*n+1)*h, (2*n+2)*h))
    ltau = np.linspace(2*n*h, (2*n+2)*h, 400//N)
    graph.plot(ltau, scipy.interpolate.lagrange(lsupp, f(lsupp))(ltau), "g--")
graph.plot(support, f(support), "k.")
graph.axis((0, np.pi, -1, 1.5))

figure.savefig("bessel_int.pdf")
