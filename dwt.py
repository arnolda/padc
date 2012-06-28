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
from scipy import *

def haar_trafo(data):
    "Diskrete Wavelettransformation mit Hilfe des Haar-Wavelets."
    # Daten mit kleinstem Integrationsschritt multiplizieren
    c = data.copy() / len(data)
    # Temporaerer Puffer, um benoetigte Werte nicht zu ueberschreiben
    ctmp = zeros(c.shape)
    width = len(c)/2
    while width >= 1:
        for n in range(width):
            tmp1 = c[2*n]
            tmp2 = c[2*n+1]
            # Detail
            ctmp[width + n] = tmp1 - tmp2
            # Downsampling
            ctmp[n]         = tmp1 + tmp2
        # Puffer zurueckschreiben
        c[:2*width] = ctmp[:2*width]
        width = width / 2
    return c

def inverse_haar_trafo(c):
    "Inverse Diskrete Wavelettransformation mit Hilfe des Haar-Wavelets"
    # Rueckgabewerte
    data = zeros(len(c))
    # phi mitnehmen auf der niedrigsten Stufe
    data[0] = c[0]
    width = 1
    cstart = 1
    while width <= len(c)/2:
        for n in range(width-1, -1, -1):
            tmp = data[n]
            data[2*n]     = tmp + width*c[cstart + n]
            data[2*n + 1] = tmp - width*c[cstart + n]
        cstart += width
        width = width * 2
    return data

# Anwendungsbeispiel

x = linspace(0,1,256)
y = cos(x)
coeff = haar_trafo(y)
yrueck = inverse_haar_trafo(coeff)
print max(abs(yrueck - y))
