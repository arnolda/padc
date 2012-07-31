# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# 2d-Poisson mittels Fourier
##############################################
from scipy import *
from numpy.fft import *
import matplotlib.pyplot as pyplot

L=1.0 # Kantenlaenge des Quadrats
N=50  # Punkte der Diskretisierung, gerade fuer FFT

h = L/N
# rho aufsetzen
rho = zeros((N, N))
s2 = 0.01 # Quadrat des Radius der Ladung im Zentrum
for i in range(N):
    for k in range(N):
        x, y = k*h, i*h
        d2 = (x-0.5*L)**2 + (y-0.5*L)**2
        rho[i, k] += 1/sqrt(2*pi*s2)*exp(-0.5*d2/s2)
# neutralisieren
rho -= sum(rho)/N**2

# Loesung per Fouriertransformation
##############################################

rho_fft = fft2(rho)
# Laplace-Operator (2 pi/n L)^2
for nx in range(N/2 + 1):
    for ny in range(N/2+1):
        if nx == 0 and ny == 0:
            rho_fft[nx,ny] = 0
        else:
            n2 = (nx**2 + ny**2)*(2*pi/L)**2
            rho_fft[      nx,   ny] /= n2
            if nx > 0 and nx < N/2:
                rho_fft[N-nx,   ny] /= n2
            if ny > 0 and ny < N/2:
                rho_fft[  nx, N-ny] /= n2
            if nx > 0 and ny > 0 and nx < N/2 and ny < N/2:
                rho_fft[N-nx, N-ny] /= n2
psi_fft = ifft2(rho_fft)
psi_fft = real(psi_fft)

# Ausgabe
#############################################

im = pyplot.imshow(psi_fft, interpolation="bilinear", origin="lower",
                   extent=(0,L,0,L))
pyplot.colorbar(im)
pyplot.show()
