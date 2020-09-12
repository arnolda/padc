# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Haar-Wavelet-Transformation
#############################
import numpy as np


def haar_trafo(data):
    """Diskrete Wavelettransformation fuer das Haar-Wavelet."""
    # Daten mit kleinstem Integrationsschritt multiplizieren
    c = data.copy() / len(data)
    # Temporaerer Puffer, um benoetigte Werte nicht zu ueberschreiben
    ctmp = np.zeros(c.shape)
    width = len(c) // 2
    while width >= 1:
        for n in range(width):
            tmp1 = c[2 * n]
            tmp2 = c[2 * n + 1]
            # Detail
            ctmp[width + n] = tmp1 - tmp2
            # Downsampling
            ctmp[n]         = tmp1 + tmp2
        # Puffer zurueckschreiben
        c[:2 * width] = ctmp[:2 * width]
        width = width // 2
    return c


def inverse_haar_trafo(c):
    """Inverse Diskrete Wavelettransformation fuer das Haar-Wavelet."""
    # Rueckgabewerte
    data = np.zeros(len(c))
    # phi mitnehmen auf der niedrigsten Stufe
    data[0] = c[0]
    width = 1
    cstart = 1
    while width <= len(c) // 2:
        for n in range(width - 1, -1, -1):
            tmp = data[n]
            data[2 * n]     = tmp + width * c[cstart + n]
            data[2 * n + 1] = tmp - width * c[cstart + n]
        cstart += width
        width = width * 2
    return data


def example():
    """Anwendungsbeispiel."""
    x = np.linspace(0, 1, 256)
    y = np.cos(x)
    yrueck = inverse_haar_trafo(haar_trafo(y))
    np.testing.assert_allclose(yrueck, y)


if __name__ == "__main__":
    example()
