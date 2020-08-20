# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Test SOR-Code
######################################

import sys
import numpy as np
sys.path.append("..")

from sor import sor_step

np.random.seed(123)

n = 20

A = np.random.uniform(0,1,n*n)
A = A.reshape((n, n))
# diagonal dominant machen
for i in range(n):
    A[i,i] = np.sum([abs(A[i,k]) for k in range(n)]) - abs(A[i,i])

# zufaelliger Zielvektor/Ladungsdichte
b = np.random.normal(0,1,n)

x = np.zeros(n)
cnt = 1
while cnt < 100 and np.linalg.norm(np.dot(A,x) - b) > 1e-5:
    sor_step(A, b, 1.5, x)
    cnt += 1

if cnt >= 100:
    raise Exception("SOR ist nicht konvergiert")
else:
    print(f"Konvergiert nach {cnt} Schritten")
