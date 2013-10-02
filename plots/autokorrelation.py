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
    ftA = fft.rfft(A)
    ftB = fft.rfft(B)
    return 1./len(A)*fft.irfft(ftA.conj()*ftB)

def acf(A):
    N = len(A)
    Ap = empty(N*2)
    Ap[:N] = (A-A.mean())/(A.std()/sqrt(2))
    Ap[N:] = 0.0
    return crosscorr(Ap, Ap)[:N]

def correlated_noise(N, tau):
    rho = exp(-1./tau)
    noise_factor = sqrt(1.0-rho**2)
    A = noise_factor*random.standard_normal(N)
    A[0] = 0.0
    for i in range(1, N):
        A[i] += rho*A[i-1]
    return A

N = 10000
nf = 1000
na = 100
x = linspace(0.0, N, N)
random.seed(17)

f0 = correlated_noise(N, 20)
acf0 = acf(f0)

f1 = correlated_noise(N, 1)
f1 += 10.0
acf1 = acf(f1)

f2 = correlated_noise(N, 1)
f2 -= 50.*x/N
f2 += 8.0
acf2 = acf(f2)

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95,wspace=0.2)
graph = figure.add_subplot(121)
graph.plot(x[:nf], f0[:nf], 'r-')
graph.plot(x[:nf], f1[:nf], 'g.')
graph.plot(x[:nf], f2[:nf], 'b:')

graph = figure.add_subplot(122)
graph.plot(x[:na], acf0[:na], 'r-')
graph.plot(x[:na], acf1[:na], 'g.')
graph.plot(x[:na], acf2[:na], 'b:')
graph.axis([0, na, -0.1, 1.1])

figure.savefig('autokorrelation.pdf')
