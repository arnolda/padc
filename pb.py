# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Einfacher iterativer Poisson-Boltzmann-Loeser
##############################################
from scipy import *
from scipy.linalg import *
import matplotlib.pyplot as pyplot

L = 1.0     # Kantenlaenge des Quadrats
N = 50      # Punkte der Diskretisierung
eps = 1     # Dielektrische Konstante
cinf = 1    # Salzkonzentration am Rand
tol = 1e-4  # maximales Residuum
h = L/N     # Schrittweite

# Fortran/NumPy-artige Indizierung
def linindex(x, y): return y + N*x
# Fixes rho und zugaenglichen Bereich aufsetzen.
# Ueberall dort, wo feste Ladung sitzt, geht es nicht hin
rho_fix, chi = zeros(N*N), zeros(N*N)
r = 0.1 # Radius der fixen Ladung im Zentrum
for i in range(N):
    for k in range(N):
        x, y = k*h, i*h
        d = sqrt((x-0.5*L)**2 + (y-0.5*L)**2)
        if d <= r:     rho_fix[linindex(i, k)] += 2/pi/r**2
        chi[linindex(i,k)] = (rho_fix[linindex(i,k)] == 0.0)
# 2d-Laplace, 0-Rand (Neutral am Rand)
Laplace=zeros((N*N, N*N))
for y in range(N):
    for x in range(N):
        eqn = linindex(x,y)
        Laplace[eqn, linindex(x,y)] = -4/h**2
        if x < N-1: Laplace[eqn, linindex(x+1,y)] = 1/h**2
        if x > 0:   Laplace[eqn, linindex(x-1,y)] = 1/h**2
        if y < N-1: Laplace[eqn, linindex(x,y+1)] = 1/h**2
        if y > 0:   Laplace[eqn, linindex(x,y-1)] = 1/h**2

# Iterativer Loeser
##############################################
psi = zeros(N*N) # Potential
while True:
    # aktuelle vollstaendige Ladungsdichte
    rho = rho_fix + cinf*2*sinh(-psi)*chi
    residual = eps*dot(Laplace, psi) + rho
    print "Residuum ist", norm(residual)
    if norm(residual) < tol: break
    psi = solve(Laplace, -rho)/eps

# Ausgabe der resultierende positiven Ionendichte
##############################################
n = (cinf*exp(-psi)*chi).reshape((N,N))
im = pyplot.imshow(n, origin="lower", extent=(0,L,0,L))
pyplot.colorbar(im)
pyplot.show()
