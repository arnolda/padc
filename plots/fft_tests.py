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
# Einfacher Fouriertest fuer RNGs
#
############################################
from scipy import *
from numpy.random import *
from numpy.fft import *

import matplotlib.pyplot as pyplot

from rng_tests import *

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.15,wspace=0.3, left=0.1,right=0.95)

rng1 = Rand()
data1 = array([ rng1.next() for x in range(1000000) ])
rng2 = Randu()
data2 = array([ rng2.next() for x in range(1000000) ])

##########################################

ft1 = rfft(data1)/len(data1)
ft2 = rfft(data2)/len(data2)

graph = figure.add_subplot(121)

rng=30

graph.bar(arange(0,rng), abs(ft1[:rng]), width=0.5, bottom=1e-5,
          color="blue", edgecolor="blue", linewidth=0)

graph.bar(0.5 + arange(0,rng), abs(ft2[:rng]), width=0.5, bottom=1e-5,
          color="red", edgecolor="white", hatch="////", linewidth=0)

graph.set_yscale("log")

##########################################

def autocorr(A):
    ft = fft(A)
    return real(ifft(ft*ft.conj()))/A.shape[0]

ac1 = autocorr(data1 - 0.5)
ac2 = autocorr(data2 - 0.5)

graph = figure.add_subplot(122)

rng=30

graph.plot(arange(0,rng), ac1[:rng], "b--",
           dashes=(0,9,2,2,2,1), dash_capstyle="butt", linewidth=2)
graph.plot(arange(0,rng), ac2[:rng], "r--",
           dashes=(8,8), dash_capstyle="butt", linewidth=2)

##########################################

figure.savefig("fft_tests.pdf")
