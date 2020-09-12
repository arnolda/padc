# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
import numpy as np


def damped_newton(f, fprime, x0, epsilon):
    """Gedaempftes Newtonverfahren."""
    xi = x0
    while True:
        # Newton-Korrektur
        di = np.linalg.solve(fprime(xi), -f(xi))
        if np.linalg.norm(di) < epsilon:
            break
        # Schrittweitendaempfung
        lmbda = 1.0
        while True:
            # neue Naeherung
            xneu = xi + lmbda * di
            if np.linalg.norm(f(xneu)) < np.linalg.norm(f(xi)):
                break
            lmbda *= 0.5
        xi = xneu

    return xi
