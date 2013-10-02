# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Sukzessive Substitution
#
############################################
from numpy import *
import matplotlib.pyplot as pyplot

# cross correlation
def crosscorr(A, B):
    ftA = fft.rfft(A - A.mean())
    ftB = fft.rfft(B - B.mean())
    return 1./len(A)*fft.irfft(ftA.conj()*ftB)

N = 100000
nf = 1000
na = 333
x = arange(0, N) * 6*pi / nf

f0 = sin(x)
f1 = sin(1.1*x)
f2 = cos(x) + 0.5

acf00 = crosscorr(f0, f0)
acf01 = crosscorr(f0, f1)
acf02 = crosscorr(f0, f2)

# Ausgabe
#########

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.92,wspace=0.2)
graph = figure.add_subplot(121)
graph.plot(x[:nf], f0[:nf], 'r-', label="$f_1(t)=\sin(t)$")
graph.plot(x[:nf], f1[:nf], 'g--',label="$f_2(t)=\sin(1,1 t)$")
graph.plot(x[:nf], f2[:nf], 'b:', label="$f_3(t)=\cos(t) + 0,5$")
graph.legend(bbox_to_anchor=(1.1, 1.1))
graph.axis([0, x[nf], -1.0, 2.5])

graph = figure.add_subplot(122)
graph.plot(x[:na], acf00[:na], 'r-', label="$C(f_1, f_1)(t)$")
graph.plot(x[:na], acf01[:na], 'g--',label="$C(f_1, f_2)(t)$")
graph.plot(x[:na], acf02[:na], 'b:', label="$C(f_1, f_3)(t)$")
graph.legend(bbox_to_anchor=(1.1, 1.1))
graph.axis([0, x[na], -1.0, 1.0])

figure.savefig('kreuzkorrelation.pdf')
