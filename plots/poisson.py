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
# 2d-Poisson mittels Fourier
##############################################
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import numpy as np

L=1.0    # Kantenlaenge des Quadrats
N=50     # Punkte der Diskretisierung, gerade fuer FFT
eps=0.02 # dielektrische Konstante

h = L/N
# rho aufsetzen
rho = np.zeros((N, N))
s2 = 0.01
for i in range(N):
    y = i*h
    for k in range(N):
        x = k*h
        d2 = (x-0.7*L)**2 + (y-0.5*L)**2
        rho[i, k] += 1/np.sqrt(2*np.pi*s2)*np.exp(-0.5*d2/s2)
        d2 = (x-0.3*L)**2 + (y-0.3*L)**2
        rho[i, k] -= 1/np.sqrt(2*np.pi*s2)*np.exp(-0.5*d2/s2)
# neutralisieren, wegen Diskretisierungsfehlern
rho -= sum(rho)/N**2

# 2d-Laplace, periodisch
##############################################
Laplace=np.zeros((N*N, N*N))

def linindex(x, y): return (x % N) + N*(y % N)

for y in range(N):
    for x in range(N):
        eqn = linindex(x,y)
        if eqn > 0:
            Laplace[eqn, linindex(x,y)] = -4/h**2
            Laplace[eqn, linindex(x+1,y)] = 1/h**2
            Laplace[eqn, linindex(x-1,y)] = 1/h**2
            Laplace[eqn, linindex(x,y+1)] = 1/h**2
            Laplace[eqn, linindex(x,y-1)] = 1/h**2
        else:
            # Normierungselement
            Laplace[0, linindex(x, y)] = 1

# Normierung und flachdruecken
rho_fd_flat = np.empty(N*N)
for y in range(N):
    for x in range(N):
        if x == 0 and y == 0:
            rho_fd_flat[linindex(x,y)] = 0
        else:
            rho_fd_flat[linindex(x,y)] = rho[x,y]

psi_fd_flat = np.linalg.solve(Laplace, -rho_fd_flat/eps)

psi_fd = np.empty((N, N))
for y in range(N):
    for x in range(N):
        psi_fd[x,y] = psi_fd_flat[linindex(x,y)]

# Loesung per Fouriertrafo
##############################################

rho_fft = np.fft.fft2(rho)

for kx in range(N//2 + 1):
    for ky in range(N//2+1):
        if kx == 0 and ky == 0:
            rho_fft[kx,ky] = 0
        else:
            k2 = (kx**2 + ky**2)*(2*np.pi/L)**2
            rho_fft[      kx,   ky] /= k2
            if kx > 0 and kx < N/2:
                rho_fft[N-kx,   ky] /= k2
            if ky > 0 and ky < N/2:
                rho_fft[  kx, N-ky] /= k2
            if kx > 0 and ky > 0 and kx < N/2 and ky < N/2:
                rho_fft[N-kx, N-ky] /= k2

psi_fft = np.fft.ifft2(rho_fft/eps)

print("Imaginaerteil, sollte verschwinden:", np.amax(np.imag(psi_fft)))
psi_fft = np.real(psi_fft)

# Ausgabe
#############################################

figure = plt.figure(figsize=(8,8))
figure.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

# rechts oben: rho
#############################################
graph = figure.add_subplot(222)

# relevanter Bereich fuer den Colorbar:
rng = (-3.5, +3.5)

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

im = graph.imshow(rho, interpolation="bilinear", origin="lower",
                  cmap = cm_densities,
                  extent=(0,L,0,L), norm = colors.Normalize(rng[0], rng[1]))
figure.colorbar(im,shrink=0.8)

# links oben: FD
#############################################
graph = figure.add_subplot(221)

print("Reichweite im Potential", np.amin(psi_fd), np.amax(psi_fd))

graph.imshow(psi_fd, interpolation="bilinear", origin="lower",
             cmap = cm_densities,
             extent=(0,L,0,L), norm = colors.Normalize(rng[0], rng[1]))

# links unten: FFT
#############################################
graph = figure.add_subplot(223)

graph.imshow(psi_fft, interpolation="bilinear", origin="lower",
                  cmap = cm_densities,
                  extent=(0,L,0,L), norm = colors.Normalize(rng[0], rng[1]))

# rechts unten: Differenz
#############################################
graph = figure.add_subplot(224)

maxpot = np.amax(abs(psi_fft))

im = graph.imshow(abs(psi_fft - psi_fft[0] - psi_fd + psi_fd[0])/maxpot, interpolation="bilinear",
                  origin="lower",
                  cmap = cm.gray, extent=(0,L,0,L))
figure.colorbar(im,shrink=0.8)

figure.savefig("poisson.pdf")
