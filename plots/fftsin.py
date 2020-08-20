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
import matplotlib.pyplot as plt

# verrauschter Sinus
####################

N=200

x = np.linspace(0, 2*np.pi, N+1)

np.random.seed(1)
y = np.sin(x) + 0.1*np.sin(10*x) + np.random.normal(0, 0.1, x.shape)

dft = np.fft.fft(y[0:N])/N
# periodischer Ringschluss fuer die Ausgabe
dft = np.concatenate((dft, (dft[0],)))

# Ausgabe
#################################

figure = plt.figure(figsize=(8,4))

graph = figure.add_subplot(121)
graph.plot(x, y, "ko", markersize=1)
graph.axis([0, 2*np.pi,-1.1,1.1])

graph = figure.add_subplot(122)
graph.set_yscale("log")
graph.bar(np.arange(N//2) - 0.5, abs(dft[0:N//2])**2, color="k", linewidth=0)
graph.axis([-2, N/2,0.00001,0.3])

figure.savefig("fftsin.pdf")
