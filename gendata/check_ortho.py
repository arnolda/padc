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
# Test verschiedene Orthogonalisierungsverfahren
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

np.random.seed(123)

# Erzeugen einer quadratischen Zufallsmatrix
# da quadratisch, geben alle Verfahren bis auf Vz dasselbe
n = 10
a = np.random.uniform(0, 1, n * n)
a = a.reshape((n, n))

from gramschmidt import gramschmidt
from householder import householder
from givens import givens

# Normierung auf r_ii > 0, nur fuer Householder noetig


def normalize(q, r):
    for i in range(min(r.shape)):
        if r[i, i] < 0:
            r[i, :] = -r[i, :]
            q[:, i] = -q[:, i]


def check(q, r, a, method, qref=None, rref=None):
    tol = 1e-10
    np.testing.assert_almost_equal(
        np.dot(
            q,
            r),
        a,
        decimal=10,
        err_msg=f"{method}: q*r != a")
    np.testing.assert_allclose(
        np.identity(
            q.shape[0]),
        np.dot(
            q.transpose().conj(),
            q),
        atol=tol,
        err_msg=f"{method}: q^Hq != I")
    np.testing.assert_allclose(
        np.tril(r, -1), 0.0, atol=tol, err_msg=f"{method}: r hat Subdiagonalelemente")

    if qref is not None and rref is not None:
        np.testing.assert_allclose(
            q, qref, atol=tol, err_msg=f"{method}: q not like reference")
        np.testing.assert_allclose(
            r, rref, atol=tol, err_msg=f"{method}: r not like reference")


q, r = gramschmidt(a)
qref = q
rref = r
check(q, r, a, "gramschmidt")

q, r = householder(a)
check(q, r, a, "householder")
normalize(q, r)
check(q, r, a, "householder", qref, rref)

q, r = givens(a)
check(q, r, a, "givens", qref, rref)
