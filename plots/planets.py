# -*- coding: utf-8 -*-
# 
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
# Himmelsmechanik
##############################################
import matplotlib.pyplot as plt
import numpy as np

import sys
# da liegt der Integrator da er Teil des Skripts ist
sys.path.append("..")

# Einheiten:
#   Laenge 1 AE (Abstand Sonne-Erde),
#   Zeit 1 Jahr
#   Erdenmassen
G= 1.1858e-4

# Erde
Merde = 1
xyerde = [0, 1]
vxyerde = [6.278, 0]
# Sonne
Msonne  = 333000
xysonne  = [0,0]
vxysonne = [0, 0]

# Mond
Mmond = 0.0123
xymond = [0, 1 + 0.00257]
vxymond = [6.278+0.2154, 0]

# alles zusammengepackt und flachgedrueckt
x0 = np.array([xysonne, xyerde, xymond]).ravel()
v0 = np.array([vxysonne, vxyerde, vxymond]).ravel()
m  = np.array([Msonne,  Merde,  Mmond])

# Gravitationsgesetz
def F(dx, m):
    r = np.linalg.norm(dx)
    return -G*m/r**3 * dx
def E(dx, m):
    r = np.linalg.norm(dx)
    return -G*m/r

# Beschleunigungen auf Objekte
def acc(t, x):
    N = len(x)//2
    x = x.reshape((N, 2))
    a = np.zeros((N, 2))
    for i in range(N):
        for k in range(N):
            if i != k:
                dx = x[i] - x[k]
                a[i] += F(dx, m[k])
    return a.ravel()

def tot_en(x, v):
    N = len(x)//2
    x = x.reshape((N, 2))
    v = v.reshape((N, 2))
    en = 0.0
    for i in range(N):
        en += 0.5*np.linalg.norm(v[i])**2*m[i]
        for k in range(N):
            if i != k:
                dx = x[i] - x[k]
                en += E(dx, m[i]*m[k])
    return en

# Loeser
#############################################

from vv import velocity_verlet

from rk import rk_explicit, rk_klassisch

# decomposing
def f(t, xv):
    n = len(xv)//2
    x = xv[:n]
    v = xv[n:]
    return np.concatenate((v, acc(t, x)))

def velocity_verlet_w_energy(acc, x0, v0, tmax, h):
    N = len(x0)
    result = velocity_verlet(acc, x0, v0, tmax, h)
    tns = result[:,0]
    xns = result[:,1:N+1]
    vns = result[:,N+1:2*N+1]

    ens = [ tot_en(xn, vn) for xn, vn in zip(xns, vns) ]
    return tns, xns, ens

def rk_explicit_w_energy(acc, x0, v0, tmax, h):
    N = len(x0)
    result = rk_explicit(rk_klassisch, f, np.concatenate((x0, v0)), tmax, h)
    tns = result[:,0]
    xns = result[:,1:N+1]
    vns = result[:,N+1:2*N+1]

    ens = [ tot_en(xn, vn) for xn, vn in zip(xns, vns) ]
    return tns, xns, ens

h = 1.0/365
tmax = 10

def unpack(tns, xns, ens):
    xerde, yerde = [], []
    xmond, ymond = [], []
    xsonne, ysonne = [], []
    for pos in xns:
        xsonne.append(pos[0])
        ysonne.append(pos[1])
        xerde.append(pos[2])
        yerde.append(pos[3])
        xmond.append(pos[4])
        ymond.append(pos[5])

    return np.array(xerde), np.array(yerde), \
        np.array(xsonne), np.array(ysonne), \
        np.array(xmond), np.array(ymond)

tns_vv, xns_vv, ens_vv = velocity_verlet_w_energy(acc, x0, v0, tmax, h)
xerde_vv, yerde_vv, \
    xsonne_vv, ysonne_vv, \
    xmond_vv, ymond_vv = unpack(tns_vv, xns_vv, ens_vv)

tns_rk, xns_rk, ens_rk = rk_explicit_w_energy(acc, x0, v0, tmax, h)
xerde_rk, yerde_rk, \
    xsonne_rk, ysonne_rk, \
    xmond_rk, ymond_rk = unpack(tns_rk, xns_rk, ens_rk)

# Ausgabe
#############################################

figure = plt.figure(figsize=(8,8))
figure.subplots_adjust(left=0.15, right=0.95,wspace=0.3)

# links, Erde um Sonne
#############################################
graph = figure.add_subplot(221)

tgtx = xerde_rk - xsonne_rk
tgty = yerde_rk - ysonne_rk
graph.plot(tgtx, tgty, "r--", dashes=(4,60), linewidth=1)
graph.plot(tgtx[20::365], tgty[20::365], "r*",
           markeredgecolor="r", markersize=6)

tgtx = xerde_vv - xsonne_vv
tgty = yerde_vv - ysonne_vv
graph.plot(tgtx, tgty, "b--", dashes=(1,10), linewidth=1)
graph.plot(tgtx[::365], tgty[::365], "bD",
           markeredgecolor="b", markersize=4)


graph.plot((0,),
           (0,), "yo", markeredgecolor="y")
graph.axis((-1.05,1.05,-1.05,1.05))

# rechts, Mond um Erde
#############################################
graph = figure.add_subplot(222)

tgtx = xmond_vv - xerde_vv
tgty = ymond_vv - yerde_vv
graph.plot(tgtx[:365//2], tgty[:365//2], "r--", linewidth=1)

tgtx = xmond_rk - xerde_rk
tgty = ymond_rk - yerde_rk
graph.plot(tgtx[:365//2], tgty[:365//2], "b:", linewidth=1)

graph.plot((0,),
           (0,), "bo", markeredgecolor="b")

graph.set_xticks((-0.002, 0, 0.002))

# unten, Energie
#############################################
graph = figure.add_subplot(212)

h = 20.0/365
tmax = 100

tns_vv2, xns_vv, ens_vv2 = velocity_verlet_w_energy(acc, x0, v0, tmax, h)
tns_rk2, xns_rk, ens_rk2 = rk_explicit_w_energy(acc, x0, v0, tmax, h)

graph.plot(tns_vv2, ens_vv2, "r-", linewidth=1)
graph.plot(tns_rk2, ens_rk2, "b--", linewidth=1)
graph.xaxis.set_label_text("Jahre")
graph.yaxis.set_label_text("Gesamtenergie")

graph.axis((0,100,-65,-55))

figure.savefig("planets.pdf")

