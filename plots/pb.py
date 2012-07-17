# 
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
# Einfacher Poisson-Boltzmann-Loeser
##############################################
from scipy import *
from scipy.linalg import *
import matplotlib.pyplot as pyplot
import matplotlib.colors as colors

L = 1.0     # Kantenlaenge des Quadrats
N = 50      # Punkte der Diskretisierung
eps = 1     # Dielektrische Konstante
cinf = 1    # Salzkonzentration am Rand
tol = 1e-4  # maximales Residuum

# Schrittweite
h = L/N
# Indizierung, Fortran/NumPy-artig
def linindex(x, y): return y + N*x

# Fixes rho und zugaenglichen Bereich aufsetzen. Ueberall dort, wo
# feste Ladung sitzt, geht es nicht hin
rho_fix, chi = zeros(N*N), zeros(N*N)
# Radius der beiden Ladungen
r = 0.1
for i in range(N):
    for k in range(N):
        x, y = k*h, i*h

        d = sqrt((x-0.5*L)**2 + (y-0.6*L)**2)
        if d <= r:     rho_fix[linindex(i, k)] += -2/pi/r**2
        d = sqrt((x-0.37*L)**2 + (y-0.4*L)**2)
        if d <= r:     rho_fix[linindex(i, k)] += 1/pi/r**2
        d = sqrt((x-0.63*L)**2 + (y-0.4*L)**2)
        if d <= r:     rho_fix[linindex(i, k)] += 1/pi/r**2

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

# Potential
psi = zeros(N*N)

converged = False
while not converged:
    rho = rho_fix + cinf*2*sinh(-psi)*chi
    residual = eps*dot(Laplace, psi) + rho

    print "Residual is", norm(residual)

    if norm(residual) < tol:
        converged = True
        break

    psi = solve(Laplace, -rho)/eps

# Ausgabe
#############################################

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.05,top=0.95,
                       left=0.05,right=0.9,wspace=0.3)

# links, positive Ionen
#############################################

# Codierung
# 0 verboten
# relevanter Bereich fuer den Colorbar:
rng = (0.8, 1.5)

# Colormap mit Schwarz fuer verbotene Dichten und weiss-blau sonst
white = 0.8/1.5
cmdict = {'red':   ((0.0,   0.0, 0.0),
                    (white, 0.0, 1.0),
                    (1.0,   0.0, 0.0),
                    ),
          'green': ((0.0,   0.0, 0.0),
                    (white, 0.0, 1.0),
                    (1.0,   0.0, 0.0),
                    ),
          'blue':  ((0.0,   0.0, 0.0),
                    (white, 0.0, 1.0),
                    (1.0,   1.0, 1.0),
                    )}
cm_densities = colors.LinearSegmentedColormap('densmap', cmdict, 256)

graph = figure.add_subplot(121)

n = cinf*exp(psi)
print "Reichweite in der Dichte", amin(n), amax(n)

n = (n*chi).reshape((N,N))

im = graph.imshow(n, interpolation="bilinear", origin="lower",
                  cmap = cm_densities,
                  extent=(0,L,0,L), norm = colors.Normalize(0, rng[1]))
figure.colorbar(im,shrink=0.66, boundaries=linspace(rng[0], rng[1], 33))

graph.text(0.5*L, 0.6*L, "2-", fontsize=14, family="bold",
           va="center", ha="center", color="#a0a0ff")
graph.text(0.37*L, 0.4*L, "+", fontsize=14, family="bold",
           va="center", ha="center", color="#ffa0a0")
graph.text(0.63*L, 0.4*L, "+", fontsize=14, family="bold",
           va="center", ha="center", color="#ffa0a0")

#        d = sqrt((x-0.37*L)**2 + (y-0.4*L)**2)
#        if d <= r:     rho_fix[linindex(i, k)] += -1/pi/r**2
#        d = sqrt((x-0.63*L)**2 + (y-0.4*L)**2)
#        if d <= r:     rho_fix[linindex(i, k)] += -1/pi/r**2

# rechts, Nettoladung
#############################################

# Codierung
# -1 im verbotenen Bereich
# relevanter Bereich fuer den Colorbar:
rng = (-0.41, +0.2)

# Entsprechende Colormap mit Rueckskalierung
# weiss bei neutral, also 0 urspruenglich
white = 1.0/(1 + rng[1])
# rot und blau symmetrisch
blue = (1.0 + rng[0])/(1 + rng[1])

cmdict = {'red':   ((0.0,   0.0, 0.0),
                    (blue,  0.0, 0.0),
                    (white, 1.0, 1.0),
                    (1.0,   1.0, 1.0),
                    ),
          'green': ((0.0,   0.0, 0.0),
                    (blue,   0.0, 0.0),
                    (white, 1.0, 1.0),
                    (1.0,   0.0, 0.0),
                    ),
          'blue':  ((0.0,   0.0, 0.0),
                    (blue,  0.0, 1.0),
                    (white, 1.0, 1.0),
                    (1.0,   0.0, 0.0),
                    )}
cm_densities = colors.LinearSegmentedColormap('densmap', cmdict, 256)

graph = figure.add_subplot(122)

n = psi.reshape((N,N))

print "Reichweite im Potential", amin(n), amax(n)

im = graph.imshow(n, interpolation="bilinear", origin="lower",
                  cmap = cm_densities,
                  extent=(0,L,0,L), norm = colors.Normalize(-1, rng[1]))
figure.colorbar(im,shrink=0.66, boundaries=linspace(rng[0], rng[1], 101))

figure.savefig("pb.pdf")
