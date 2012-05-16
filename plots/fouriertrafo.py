# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Fourierreihen
#
############################################
from scipy import *
import numpy.fft as fft
import math
import matplotlib.pyplot as pyplot

N=10000
T=2*pi*5

x = linspace(-T/2, T/2, N)
o = linspace(0,N*2*pi/T,N)

def contfft(data):
    return ((-1)**arange(N))*T/sqrt(2*pi)/N*fft.fft(data)

def icontfft(data):
    return fft.ifft(sqrt(2*pi)*N/T*((-1)**arange(N))*data)

# Beispiel 1: Glaettung einer hochfrequenten Stoerung
# mit einer Gaussglocke
yf = cos(x) + 0.25*sin(7.5*x)
f  = contfft(yf)

s=0.5
yg = 1/sqrt(2*pi)/s*exp(-x**2/2/s**2)
g  = contfft(yg)

h=sqrt(2*pi)*f*g
yh=real(icontfft(h))

figure = pyplot.figure(figsize=(8,8))

graph = figure.add_subplot(221)
graph.plot(x, yf, "k", linewidth=0.5)
graph.plot(x, yg, "b--", linewidth=1.5)
graph.plot(x, yh, "r:", linewidth=2)
graph.axis([-pi,pi,-1.3,1.3])

graph = figure.add_subplot(222)
graph.plot(o, abs(f), "k", linewidth=0.5)
graph.plot(o, abs(g), "b--", linewidth=1)
graph.plot(o, abs(h), "r:", linewidth=2)
graph.axis([0,10,0,2])

# Beispiel 2: Mittelpassfilter
# konstruiert im k-Raum
o0 = 7.5
s = 1
# k-Raum-Konstruktion
g  = 1/sqrt(2*pi)/s*exp(-(o-o0)**2/2/s**2)
# g = (o < 9) * (o > 6)

# Symmetrisierung fuer ein reelles Signal
g = g + concatenate(((0,), g[N:0:-1]))

# Ruecktrafo, das real gleicht numerische Fehler im Imaginaerteil aus. 
yg = icontfft(g)
print "max imag=", max(imag(yg))
yg = real(yg)

h=sqrt(2*pi)*f*g
yh=real(icontfft(h))

graph = figure.add_subplot(223)
graph.plot(x, yf, "k", linewidth=0.5)
graph.plot(x, yg, "b--", linewidth=1.5)
graph.plot(x, yh, "r:", linewidth=2)
graph.axis([-pi,pi,-1.3,1.3])

graph = figure.add_subplot(224)
graph.plot(o, abs(f), "k", linewidth=0.5)
graph.plot(o, abs(g), "b--", linewidth=1.5)
graph.plot(o, abs(h), "r:", linewidth=2)
graph.axis([0,10,-1,1])

figure.savefig("fouriertrafo.pdf")
