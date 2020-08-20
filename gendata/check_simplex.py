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
# Test des Simplexalgorithmus
##############################################

import sys
import numpy as np
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")

from simplex import simplex

# Abbildung 9.3.1
# die zweite Bedingung ist offenbar sinnlos
A = np.array(((-1, -1, -1),
           (2, 2, 2)))
b = np.array((-1, 2))
for c, res in (((0.5, 0.5, 0), (0.0,0.0,1.0)),
               ((-0.5, 0.5, 0), (1.0,0.0,0.0)),
               ((-0.1, -0.5, 0), (0.0,1.0,0.0))):
    answer = simplex(np.array(c), A, b)
    if not np.allclose(answer, np.array(res)):
        raise Exception("Simplex hat das Minimum nicht gefunden")

# die guten, alten Beispiele aus den Numerik-Uebungen
# Bauer 2a)
A = np.array(((1,   1, 1,0,0),
           (40,120, 0,1,0),
           ( 7, 12, 0,0,1)),dtype=float)
b = np.array((40, 2400, 312), dtype=float)
c = np.array((-100, -250, 0, 0, 0), dtype=float)

answer = simplex(np.array(c), A, b)
if not np.allclose(answer, np.array((24, 12, 4, 0, 0), dtype=float)):
    raise Exception("Simplex hat das Minimum nicht gefunden")

# Bauer 2b)
A = np.array(((1,   1,  1, 1,0,0),
           (40,120, 10, 0,1,0),           
           ( 7, 12,  7, 0,0,1),
           ( 1,  0, -1, 0,0,0),
           ( 7,  0, -7, 0,0,0),
           ), dtype=float)
b = np.array((40, 2400, 312, 0, 0), dtype=float)
c = np.array((-100, -250, -120, 0, 0, 0), dtype=float)

answer = simplex(np.array(c), A, b)
if not np.allclose(answer, np.array((8, 16 + 2./3., 8, 7 + 1./3., 0, 0), dtype=float)):
    raise Exception("Simplex hat das Minimum nicht gefunden")

# 2c), unbeschraenkte zulaessige Menge
A = np.array(((1, -1, -1, 0),
           (2, -1,  0,-1)), dtype=float)
b = np.array((-5, 10), dtype=float)
c = np.array((-1, -1, 0, 0), dtype=float)

try:
    answer = simplex(np.array(c), A, b)
except Exception as e:
    if str(e) != "zulaessige Menge unbeschraenkt!":
        raise Exception("Simplex hat unbeschraenkte zulaessige Menge nicht erkannt")
else:
    raise Exception("Simplex hat unbeschraenkte zulaessige Menge nicht erkannt")

# 2d), leere zulaessige Menge
A = np.array((( 1, 1, -1, 0, 0),
           (-1, 1,  0,-1, 0),
           (-2, 2,  0, 0, 1)), dtype=float)
b = np.array((1, 2, 3), dtype=float)
c = np.array((1, 1, 0, 0,0), dtype=float)

try:
    answer = simplex(np.array(c), A, b)
except Exception as e:
    if str(e) != "zulaessige Menge ist leer!":
        raise Exception("Simplex hat leere zulaessige Menge nicht erkannt")
else:
    raise Exception("Simplex hat leere zulaessige Menge nicht erkannt")
