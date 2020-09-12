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


def armijo_steepest_descent(f, gradf, x0, alpha=0.1, rho=0.5,
                            tol=1e-5, maxiter=1000):
    x = x0.copy()
    fx = f(x)
    grad = gradf(x)
    step = 0
    while np.linalg.norm(grad) > tol:
        d = -grad
        # Armjio-Schrittweite berechnen
        lmbda = 1
        xneu = x + d
        grad2 = np.dot(grad, d)
        while f(xneu) > fx + alpha * lmbda * grad2:
            lmbda = rho * lmbda
            xneu = x + lmbda * d
            step += 1
            if step > maxiter:
                # Abbruch, letzte Naeherung zurueckgegeben
                return x
        # mit dem gefundenen Schritt vorwaerts
        x = xneu
        fx = f(x)
        grad = gradf(x)

    return x
