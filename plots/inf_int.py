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
# pi/2 = int_0^infty 1/(1+x**2) mit Romberg
#
############################################
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

def f(x):
    return 1.0/(1.0+x**2)

def romberg(f, a, b, kmax):
    # Romberg-Stuetzpunkte erzeugen
    ############################################

    h = b - a
    tf = h*(0.5*f(a) + 0.5*f(b))

    tf_list = [ tf ]
    h_list  = [ h ]

    for i in range(0, kmax):
        # print("Tf(%f) = " % h, tf
        tf = 0.5*tf + 0.5*h*sum(f(a + (np.arange(0,2**i) + 0.5)*h))
        h *= 0.5
        tf_list.append(tf)
        h_list.append(h)

    tf = np.array(tf_list)
    h = np.array(h_list)

    ############################################

    res = scipy.interpolate.lagrange(h**2, tf)(0)

    # print("Tf(0) = ", res

    return res

# inverses integral
inverse = 2*romberg(f, 0, 1, 6)
print("Integral kmax = ", 6, "invers ", inverse, " Fehler = ", abs(inverse - np.pi/2))

# integral 0 bis gross direkt
b_list = np.linspace(1,1000,100)
r_list_6 = []
r_list_8 = []
for b in b_list:
    result = romberg(f, 0, b, 6)
    r_list_6.append(result)
    result = romberg(f, 0, b, 8)
    r_list_8.append(result)

print("min Fehler kmax = 6 direkt", min(abs(np.array(r_list_6) - np.pi/2)))
print("min Fehler kmax = 8 direkt", min(abs(np.array(r_list_8) - np.pi/2)))

figure = plt.figure(figsize=(4,4))
graph = figure.add_subplot(111)

graph.plot(b_list, r_list_6, "r-")
graph.plot(b_list, r_list_8, "b:")
graph.plot((min(b_list), max(b_list)), (np.pi/2, np.pi/2), "k--")
graph.axis((min(b_list), max(b_list), np.pi/2-.5, np.pi/2+.5))
graph.set_xlabel("b")

figure.savefig("inf_int.pdf")
