# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Demonstration QR-Algorithmus
######################################

import math
import sys
import numpy as np

# Fibonacci-Matrix
A = np.array(((0, 1), (1, 1)))

# QR-Algorithmus
################################
print("QR-Algorithmus")

Ak = A.copy()
I = np.identity(Ak.shape[0])
for i in range(6):
    print("Iteration {}, Ak ist".format(i))
    print(Ak)

    shift = Ak[-1, -1]
    Q, R = np.linalg.qr(Ak - shift * I)
    Ak = np.dot(R, Q) + shift * I

print("Eigenwerte sollten {} und {} sein.".format(
    0.5 * (1 - math.sqrt(5)), 0.5 * (1 + math.sqrt(5))))

# Inverse Iteration
################################

# Iteration ueber die Eigenwerte
for l in np.diag(Ak):
    print("Inverse Iteration zum Eigenwert {}".format(l))
    Ashift = A - l * I
    # Startwert
    xk = np.ones(A.shape[0])
    for i in range(2):
        xk = np.linalg.solve(Ashift, xk)
        xk = xk / np.linalg.norm(xk)
        print("Iteration {}, xk ist".format(i))
        print(xk)

    print("Der Eigenvektor sollte sein:")
    v = np.array([1.0, l])
    print(v / np.linalg.norm(v))
