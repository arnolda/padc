# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
from numpy import *
import numpy.fft as fft
import matplotlib.pyplot as pyplot

# Gauss vs geshifteter Gauss
############################

T = 16.0
N = 32
ncont = 1000

def gauss(x):
    return 1.0/sqrt(2.0*pi) * exp(-0.5 * x**2)

# diskretes Gitter
t_grid = linspace(-T/2, T/2, N, endpoint=False)
t_grid_shifted = linspace(0, T, N, endpoint=False)

# "kontinuierliche" Naeherung zum Malen
t_cont = linspace(-T/2, T/2, ncont)
t_cont_shifted = linspace(0, T, ncont)

# Ausgabe
#################################

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, bottom=0.12, right=0.98, top=0.95,wspace=0.15,hspace=.4)

# einfacher Gauss und Diskretisierung
graph = figure.add_subplot(121)

g = gauss(t_cont)
graph.plot(t_cont, g, 'k-')
g = gauss(t_grid)
graph.plot(t_grid, g, 'bo', clip_on=False)
graph.axis([-T/2, T/2, 0, 0.5])
graph.set_xlabel("t")

# und die Transformierte dazu, mit und ohne Twiddlefaktor
graph = figure.add_subplot(122)

# erst mit
ixs = arange(-N/2, N/2)
ghat = T/N * (-1)**ixs * 1./sqrt(2*pi) * fft.fft(g)
omega=2*pi/T
graph.plot(omega*ixs, real(ghat[ixs]), 'bo-')

# dann ohne
ghat = T/N * 1./sqrt(2*pi) * fft.fft(g)
graph.plot(omega*ixs, real(ghat[ixs]), 'rD-')

graph.axis([-omega*N/2, omega*N/2, -0.5, 0.5])
graph.set_xlabel("$\omega$")

figure.savefig("dft.pdf")
