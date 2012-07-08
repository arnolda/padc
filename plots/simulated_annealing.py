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
# LJ-Simulated Annealing
##############################################

from scipy import *
from scipy.linalg import norm
import matplotlib.pyplot as pyplot
import matplotlib.colors as colors
import matplotlib.cm as cm
import sys
sys.path.append("..")

# Quadratkantenlaenge
L=12
# Faktoren
epsilons = ("0.1", "0.05", "0.01", "0.001")

path="./"

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.13,wspace=0.2,hspace=0.2)

posplus  = {}
posminus = {}
E        = {}
for eps in epsilons:
    posplus[eps]  = loadtxt(path + "lj_snapshot_A_%s.coord" % eps).transpose()
    posminus[eps] = loadtxt(path + "lj_snapshot_B_%s.coord" % eps).transpose()
    E[eps] = loadtxt(path + "lj_E_%s.data" % eps).transpose()
    # Zeitschritte statt Zeit
    E[eps][0] /= 0.005

# Energien
##################################

graph_E = figure.add_subplot(122)

# LJ-Potential
##################################

for eps, place, style in zip(epsilons, (241,242,245,246),
                             ("k:", "g-", "b--", "m-.")):
    graph = figure.add_subplot(place)

    x, y = posplus[eps]
    graph.plot(x, y, ".",
               markerfacecolor="#ff0000",
               markeredgecolor="#ff0000", markersize=4)
    x, y = posminus[eps]
    graph.plot(x, y, ".",
               markerfacecolor="#000060", 
               markeredgecolor="#000060", 
               markersize=4)
    graph.set_xticks(())
    graph.set_yticks(())
    graph.set_title("$\epsilon=%s$" % eps, position=(0.5,-0.2))
    t, energy = E[eps]
    graph_E.plot(t, energy, style)
    graph_E.plot(t[-1], energy[-1], "ko", markersize=4)

graph_E.xaxis.set_label_text("Simulationsschritt")
graph_E.set_xscale("log")
graph_E.axis((1e3,1e6, -75, 0))

figure.savefig("simulated_annealing.pdf")
