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
import numpy.fft as fft
import matplotlib.pyplot as pyplot

# zu hochfrequenter Sinus
####################

N=14

x = linspace(0, 2*pi, N+1)

fnyquist = N/2

print "fnyquist=", fnyquist

f = 13

def ref(x):
    return sin(f*x)

y = ref(x)

dft = fft.fft(y[0:N])/N

ifft = fft.ifft(dft*N)
print max(imag(ifft))
# nur formal, ist bis auf numerische Fehler reell (siehe print)
ifft = real(ifft)
# erstes Element verdoppeln, periodischer Ringschluss
ifft = concatenate((ifft, (ifft[0],)))

# Ausgabe
#################################

px = linspace(0, 2*pi, 100*N)

figure = pyplot.figure(figsize=(8,4))

graph = figure.add_subplot(121)
graph.plot(px, ref(px), "r-")
graph.plot(x, ifft, "ro", markersize=3)
graph.plot(x, ifft, "k:")
graph.axis([0, 2*pi,-1.1,1.1])

graph = figure.add_subplot(122)
graph.bar(arange(N/2) - 0.5, abs(dft[0:N/2])**2, color="k", linewidth=0)
graph.axis([0, N/2,0,0.5])

figure.savefig("fftalias.pdf")
