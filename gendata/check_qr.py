# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Test QR-Code
######################################

from scipy import *
from scipy.linalg import *
from numpy.random import *
import sys
sys.path.append("..")

from qr import qr_eigenwerte
from inverse_iteration import inverse_iteration

seed(123)

A = uniform(0,1,10*10)
A = A.reshape((10, 10))
# symmetrisieren fuer reelle Eigenwerte
for i in range(10):
    for k in range(10):
        A[i, k] = A[k,i]

ews = qr_eigenwerte(A, 1e-5)

for ew in ews:
    x = inverse_iteration(A, ew, 1e-5) 
    if norm(ew*x - dot(A, x))/norm(x) > 1e-5:
        raise Exception("qr oder inverse_iteration haben ein Problem")

print "Alle Eigenvektoren und Eigenwerte passen"
