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
# Optische Tests, Quadrat- und Wuerfeltest
#
############################################
from scipy import *
from numpy.random import *
import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D

from rng_tests import *

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(bottom=0.15,wspace=0.3, left=0.1,right=0.95)

##########################################

rng = Randu()
data = [ rng.next() for x in range(2000) ]

graph = figure.add_subplot(121)

graph.plot(data[:-1], data[1:], ".", markersize=2)

##########################################

rng = Randu()
data = [ rng.next() for x in range(2000) ]

rect = figure.add_subplot(1, 2, 2).get_position()
graph = Axes3D(figure, rect)

graph.scatter(data[:-2], data[2:], data[1:-1], s=1, marker="o",
              edgecolors="blue", facecolors="blue")

figure.savefig("optical_test.pdf")
