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
# Test des Poisson-Loesers
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import poisson

L = 5.0
N = 200

# Test Ladungsdichte
##############################################

rho = poisson.charge_density(0.0001, L, N)

# Dichte ist neutral
total = 0
for i in range(rho.shape[0]):
    for k in range(rho.shape[1]):
        total += rho[i, k]

np.testing.assert_allclose(total, 0, atol=1e-7)

# Test Loeser an Plattenkondensator
##############################################

rho = np.zeros((N, N))
for i in range(N):
    rho[0, i] = 1 / N
    rho[N // 2, i] = -1 / N

psi = poisson.solve(rho, L, N)

import matplotlib.pyplot as plt
ausgabe = plt.figure(figsize=(8, 4))

x = np.arange(0, L, L / N)

np.testing.assert_allclose(
     psi[:, 0], 0.5 / L *(abs(x - 0.5*L) -  0.25*L), atol=0.01)
