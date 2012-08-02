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

L    = 2.0     # Kantenlaenge des Quadrats
N    = 50      # Punkte der Diskretisierung
lb   = 0.7     # Bjerrumlaenge
cinf = 0.2*0.6 # Molaritaet der Salzkonzentration am Rand
tol  = 1e-2    # maximales Residuum
h    = L/N     # Schrittweite

# Schrittweite
h = L/N
# Indizierung, Fortran/NumPy-artig
def linindex(x, y): return y + N*x

# Fixes rho und zugaenglichen Bereich aufsetzen. Ueberall dort, wo
# feste Ladung sitzt, geht es nicht hin
rho_fix, chi = zeros(N*N), zeros(N*N)
# Radius der drei Ladungen
r = 0.15
# Position und Ladung
chrg=(( 0,   0.7, -2, "2-", "#8080ff"),
      ( 0.5,-0.7,  1, "+",  "#ff8080"),
      (-0.5,-0.7,  1, "+",  "#ff8080"))

for i in range(N):
    for k in range(N):
        x, y = k*h -0.5*L, i*h-0.5*L
        
        for px, py, q, str, col in chrg:
            d = sqrt((x-px)**2 + (y-py)**2)
            if d <= r:
                rho_fix[linindex(i, k)] += q/pi/r**2
                
        chi[linindex(i,k)] = (rho_fix[linindex(i,k)] == 0.0)

print "Nettoladung", sum(rho_fix*h**2)

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

lu = lu_factor(Laplace)

while True:
    rho = rho_fix + cinf*2*sinh(-psi)*chi
    residual = dot(Laplace, psi)/(4*pi*lb) + rho

    print "Residual is", norm(residual)

    if norm(residual) < tol:
        break

    psi = lu_solve(lu, -4*pi*lb*rho)

# Direkter Loeser, nicht PB
##############################################

psi_fix = solve(Laplace, -4*pi*lb*rho_fix)

# Ausgabe
#############################################

figure = pyplot.figure(figsize=(8,8))
figure.subplots_adjust(bottom=0.05,top=0.95, left=0.05,right=0.9)

print "Reichweite in der festen Ladungsdichte", amin(rho_fix), amax(rho_fix)
n = 2*cinf*sinh(psi)*chi
print "Reichweite in der Ionenladungsdichte ausserhalb der Ladungen", amin(n), amax(n)
print "Reichweite im PB-Potential ausserhalb der Ladungen", amin(psi*chi), amax(psi*chi)
print "Reichweite im reinen Potential ausserhalb der Ladungen", amin(psi_fix*chi), amax(psi_fix*chi)

# links oben, negative Ionen
#############################################

# Codierung
# relevanter Bereich fuer den Colorbar:
rng = (0, 1.0)

# Colormap weiss-blau
cmdict = {'red':   ((0.0,   1.0, 1.0),
                    (1.0,   0.0, 0.0),
                    ),
          'green': ((0.0,   1.0, 1.0),
                    (1.0,   0.0, 0.0),
                    ),
          'blue':  ((0.0,   1.0, 1.0),
                    (1.0,   1.0, 1.0),
                    )}
cm_densities = colors.LinearSegmentedColormap('densmap', cmdict, 512)

graph = figure.add_subplot(221)

n = cinf*exp(psi)*chi
print "Reichweite in der Dichte", amin(n), amax(n)

n = n.reshape((N,N))

im = graph.imshow(n, interpolation="bilinear", origin="lower",
                  cmap = cm_densities,
                  extent=(0,L,0,L), norm = colors.Normalize(0, rng[1]))
figure.colorbar(im,shrink=0.66)

for px, py, q, chstring, color in chrg:
    graph.text(0.5*L+px, 0.5*L+py, chstring, fontsize=14, weight="bold",
               va="center", ha="center", color=color)

# rechts oben, positive Ionen
#############################################

# Colormap rot statt blau
cmdict = { 'red':   cmdict['blue'],
           'green': cmdict['green'],
           'blue':  cmdict['red'] }

cm_densities = colors.LinearSegmentedColormap('densmap', cmdict, 256)

graph = figure.add_subplot(222)

n = cinf*exp(-psi)*chi
print "Reichweite in der Dichte", amin(n), amax(n)

n = n.reshape((N,N))

im = graph.imshow(n, interpolation="bilinear", origin="lower",
                  cmap = cm_densities,
                  extent=(0,L,0,L), norm = colors.Normalize(0, rng[1]))
figure.colorbar(im,shrink=0.66, extend="max")

for px, py, q, chstring, color in chrg:
    graph.text(0.5*L+px, 0.5*L+py, chstring, fontsize=14, weight="bold",
               va="center", ha="center", color=color)

# links unten, Potential
#############################################

# relevanter Bereich fuer den Colorbar:
rng = (-1.5, +2.5)

# Entsprechende Colormap mit Rueckskalierung
# weiss bei neutral, also 0 urspruenglich
white = -rng[0]/(rng[1] - rng[0])

cmdict = {'red':   ((0.0,   0.0, 0.0),
                    (white, 1.0, 1.0),
                    (1.0,   1.0, 1.0),
                    ),
          'green': ((0.0,   0.0, 0.0),
                    (white, 1.0, 1.0),
                    (1.0,   0.0, 0.0),
                    ),
          'blue':  ((0.0,   1.0, 1.0),
                    (white, 1.0, 1.0),
                    (1.0,   0.0, 0.0),
                    )}
cm_densities = colors.LinearSegmentedColormap('densmap', cmdict, 256)

graph = figure.add_subplot(223)

n = psi.reshape((N,N))

im = graph.imshow(n, interpolation="bilinear", origin="lower",
                  cmap = cm_densities,
                  extent=(0,L,0,L), norm = colors.Normalize(rng[0], rng[1]))
figure.colorbar(im,shrink=0.66)

# rechts unten, Potentiallinien
#############################################

graph = figure.add_subplot(4,2,6)

psi = psi.reshape((N,N))
psi_fix = psi_fix.reshape((N,N))

graph.plot(linspace(0, L, N+1), concatenate((psi[:, N/2], (0,))), "r--")
graph.plot(linspace(0, L, N+1), concatenate((psi_fix[:, N/2], (0,))), "k-")
#graph.axis(())

graph = figure.add_subplot(4,2,8)

graph.plot(linspace(0, L, N+1), concatenate((psi[N/4, :], (0,))), "r--")
graph.plot(linspace(0, L, N+1), concatenate((psi_fix[N/4, :], (0,))), "k-")
#graph.axis(())

figure.savefig("pb.pdf")
