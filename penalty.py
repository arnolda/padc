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

from armijo import armijo_steepest_descent




# Das sollte grade nicht in optimierung.tex auftauchen
# naechste ist Zeile 20, und das sollte so bleiben
def penalty(f, gradf, g, gradg, x0, sigma_exp=2, eps_exp=-2, tol=1e-5):
    x = x0.copy()
    converged = False
    step = 1.0
    while not converged:
        sigma = step**sigma_exp
        eps   = step**eps_exp

        # Funktion + aktuelle Straffunktionen
        def q(x):
            res = f(x)
            for lg in g(x):
                res += minimum(0, sigma*lg - eps)**2
            return res
        def gradq(x):
            res = gradf(x)
            for lg, lgg in zip(g(x), gradg(x)):
                res += 2*minimum(0, sigma*lg - eps)*sigma*lgg
            return res

        # steilster Abstieg mit Armijo-Schrittweite
        x = armijo_steepest_descent(q, gradq, x, tol)
        if min(g(x)) > -tol:
            converged = True
        step += 1.0

    return x
