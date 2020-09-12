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


def meanvar(x):
    summe = sum(x)
    summe2 = sum(v * v for v in x)
    mittel = summe / len(x)
    sigma2 = (summe2 / len(x) - mittel**2) * len(x) / (len(x) - 1)
    fehler = np.sqrt(sigma2 / len(x))
    return mittel, sigma2, fehler
