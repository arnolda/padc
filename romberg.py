# Romberg-Integration
############################################
from scipy import *
from scipy.special import *
from scipy.interpolate import *
import matplotlib.pyplot as pyplot


def romberg(f, a, b, kmax):
    """Romberg-Integration der Funktion f von a bis b.
    Die Schrittweiten reichen von (b-a) bis (b-a)*2**-kmax.
    """

    # Romberg-Stuetzpunkte erzeugen
    ############################################

    # aktuelle Schrittweite
    h = b - a
    # aktuelle Naeherung (fuer den Anfang Trapezregel auf 2 Punkten)
    tf = h*(0.5*f(a) + 0.5*f(b))

    # Stuetzpunkte fuer die Interpolation, verwendete Schrittweiten
    # und abgeschaetzte Integrale
    h_list  = [ h ]
    tf_list = [ tf ]

    for i in range(0, kmax):
        # Schaetzung mit der halben Schrittweite
        tf = 0.5*tf + 0.5*h*sum(f(a + (arange(0,2**i) + 0.5)*h))
        # Schrittweite nachziehen
        h *= 0.5
        # und alles speichern
        tf_list.append(tf)
        h_list.append(h)

    tf = array(tf_list)
    h = array(h_list)

    # Auswerten des interpolierenden Polynoms im Nullpunkt
    ######################################################

    return lagrange(h**2, tf)(0)

print "pi = %.16f" % (4*romberg(lambda x: 1.0/(1.0 + x**2), 0.0, 1.0, 8))
