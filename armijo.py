# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
from scipy import *
from scipy.linalg import *

def armijo_steepest_descent(f, gradf, x0, alpha=0.1, rho=0.5, tol=1e-5):
    x = x0.copy()
    fx = f(x)
    grad = gradf(x)
    while norm(grad) > tol:
        d = -grad
        # Armjio-Schrittweite berechnen
        # Mindestschrittweite
        minlmbda = tol**2/norm(d)
        lmbda = 1
        xneu = x + d
        grad2 = dot(grad, d)
        while f(xneu) > fx + alpha*lmbda*grad2 and lmbda > minlmbda:
            lmbda = rho*lmbda
            xneu = x + lmbda*d
        # Abbruch, wenn kein Abstieg mehr
        if f(xneu) >= fx:
            break
        # mit dem gefundenen Schritt vorwaerts
        x = xneu
        fx = f(x)
        grad = gradf(x)

    return x

