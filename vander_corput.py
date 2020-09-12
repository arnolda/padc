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


def vanderCorput(N, p):
    # zu wandelnde Zahlen
    numbers = np.arange(1, int(N) + 1)
    # bitumgekehrtes Ergebnis
    result = np.zeros(N)
    # Wert der aktuellen, inversen Stelle
    frac = 1.0 / p

    # solange die groesste Zahl noch Stellen hat
    while numbers[-1] > 0:
        # unterste Stelle abschneiden
        digit = numbers % p
        numbers = numbers // p
        # ... und zum Ergebnis hinzufuegen
        result += frac * digit
        frac /= p

    return np.array(result)
