# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Spline-Interpolation der Rungefunktion
#
############################################
import scipy.interpolate as ip
import numpy as np
import matplotlib.pyplot as pyplot

# die zu interpolierende Funktion
def runge(x): return 1/(1+x**2)

# interpolating splines

x = np.linspace(-5,5,11)
y = runge(x)

cubsp11 = ip.interp1d(x, y, kind = "cubic", bounds_error=False)

x = np.linspace(-5,5,7)
y = runge(x)

linsp7 = ip.interp1d(x, y, kind = "linear", bounds_error=False)
cubsp7 = ip.interp1d(x, y, kind = "cubic", bounds_error=False)

# Ausgabe
#################################
px = np.linspace(-5.1, 5.1, 200)

figure = pyplot.figure(figsize=(4,4))

graph = figure.add_subplot(111)
graph.plot(px, runge(px), "k", linewidth=0.5)
graph.plot(px, cubsp11(px), "r--",linewidth=1.5)
graph.plot(px, cubsp7(px), "g:")
graph.plot(px, linsp7(px), "b-.")
graph.axis([-5.1,5.1,0,1.1])

figure.savefig("splines.pdf")
