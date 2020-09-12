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
# Test des Poisson-Boltzmann-Loesers
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

import pb

rho_fix, chi = pb.fixed_charge_density()

# Ladung ist nicht zugänglich
total = 0
for i in range(pb.N):
    for k in range(pb.N):
        idx = pb.linindex(i, k)
        if chi[idx]:
            np.testing.assert_allclose(rho_fix[idx], 0)
        total += rho_fix[idx]*pb.h**2

np.testing.assert_allclose(total, 1, atol=0.1)

psi = pb.solve(rho_fix, chi)

for i in range(pb.N):
    idx = pb.linindex(i, 0)
    np.testing.assert_allclose(psi[idx], 0, atol=0.1)
    idx = pb.linindex(i, pb.N - 1)
    np.testing.assert_allclose(psi[idx], 0, atol=0.1)
    idx = pb.linindex(0, i)
    np.testing.assert_allclose(psi[idx], 0, atol=0.1)
    idx = pb.linindex(pb.N - 1, i)
    np.testing.assert_allclose(psi[idx], 0, atol=0.1)
