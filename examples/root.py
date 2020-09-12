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
"""
Newtonverfahren zur Wurzelnaeherung
"""


def main(a, t):
    "Berechne k-te Wurzel aus a mit Genauigkeit t"
    k = 2

    xcurr = a
    xlast = 0.0
    step = 1

    print(f"x[{step}] = {xcurr:.15f}")

    # Schleife, bis das Ergebnis konvergiert ist
    while abs(xcurr - xlast) > t:
        xlast = xcurr
        # Newtonschritt
        xcurr = (1 - 1.0 / k) * xlast + a * 1.0 / (k * xlast**(k - 1))
        print(f"x[{step}] = {xcurr:.15f}")
        step += 1


if __name__ == "__main__":
    main(2.0, 1e-10)
