# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
from numpy import *

def rk_explicit(verfahren, f, y0, tmax, h):
    def step(hc, hA, hb, f, yn, tn):
        "ein einzelner Schritt"
        k = []
        for i in range(len(hc)):
            k.append(f(tn + hc[i], yn + dot(hA[i,:i], k[:i])))
        return yn + dot(hb, k)

    # Zur Beschleunigung skalierte Parameter vorberechnen
    hA = h*verfahren['A']
    hc = h*verfahren['c']
    hb = h*verfahren['b']

    # Startwert und -zeit
    tn = 0.0
    yn = y0.copy()
    # Ergebnisvektor mit Zeit und Punkten
    result = [ concatenate(((tn,), yn.copy())) ]
    while tn < tmax:
        yn = step(hc, hA, hb, f, yn, tn)
        tn += h
        result.append(concatenate(((tn,), yn.copy())))
    return array(result)

# Butchertableau fuer das explizite Eulerverfahren
euler = { 'c': array((0,)),
          'A': array(((0,),)),
          'b': array((1,)) }

# Butchertableau fuer das klassische Runge-Kutta-Verfahren
rk_klassisch = { 'c': array((0,0.5,0.5,1)),
                 'A': array(((0  ,  0, 0),
                             (0.5,  0, 0),
                             (0  ,0.5, 0),
                             (0  ,  0, 1),
                             )),
                 'b': array((1./6,1./3,1./3,1./6))}
