# -*- coding: utf-8 -*-
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
# Iterative Gleichungsloeser
#
##############################################
from scipy import *
from scipy.linalg import *
from numpy.random import *
import matplotlib.pyplot as pyplot

seed(123)

# DGL-Diskretisierungsschritte
N=5
# Matrixgroesse
n=N*N

# Diskretisierte 2d-Laplace-Gleichung
##############################################

def linindex(x, y):
    global N
    return (x % N) + N*(y % N)
# Schrittweite
h=1
# Matrix
Ap=zeros((N*N, N*N))

eqn=0
for y in range(N):
    for x in range(N):
        if eqn > 0:
            Ap[eqn, linindex(x,y)] = -4/h**2
            Ap[eqn, linindex(x+1,y)] = 1/h**2
            Ap[eqn, linindex(x-1,y)] = 1/h**2
            Ap[eqn, linindex(x,y+1)] = 1/h**2
            Ap[eqn, linindex(x,y-1)] = 1/h**2
        else:
            # Normierung
            Ap[0, 0] = 1
        eqn += 1

# zufaelliger Zielvektor/Ladungsdichte
bp = normal(0,1,N*N)

# zum Vergleich, strikt diagonal dominant
##############################################

A = normal(0,1,n*n)
A = A.reshape((n,n))
# strikt diagonal dominant machen
for i in range(n):
    A[i,i] = 1 + sum([abs(A[i,k]) for k in range(n)]) - abs(A[i,i])

# zufaelliger Zielvektor/Ladungsdichte
b = normal(0,1,n)

# Verfahren
##############################################

def jacobi_step(A, b, x):
    xnew = zeros_like(x)
    n = A.shape[0]
    for j in range(n):
        s = b[j]
        s -= sum([A[j,k]*x[k] for k in range(j)])
        s -= sum([A[j,k]*x[k] for k in range(j+1,n)])
        xnew[j] = s/A[j,j]
    return xnew

def sor_step(A, b, omega, x):
    n = A.shape[0]
    for j in range(n):
        s = b[j]
        # diese Werte sind schon neu
        s -= sum([A[j,k]*x[k] for k in range(j)])
        # und diese Werte sind noch alt
        s -= sum([A[j,k]*x[k] for k in range(j+1,n)])
        # Ueberschreiben von x[j]!
        x[j] = omega*s/A[j,j] + (1-omega)*x[j]

# Ausgabe
##########################################

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)

##########################################

graph = figure.add_subplot(121)

graph.set_yscale("log")

# Anzahl Jacobi/etc Schritte
steps = 40

# gutartige Matrix A
xjacobi  = zeros(n)
xjacobis  = [ xjacobi ]
xgs      = zeros(n)
xgss      = [ xgs.copy() ]
# nicht strikt diagdom Ap
xgsp = zeros(N*N)
xgssp = [ xgsp ]
xsorp    = zeros(N*N)
xsorsp    = [ xsorp.copy() ]

for step in range(steps):
    xjacobi = jacobi_step(A, b, xjacobi)
    xjacobis.append(xjacobi)
    sor_step(A, b, 1, xgs)
    xgss.append(xgs.copy())

    sor_step(Ap, bp, 1, xgsp)
    xgssp.append(xgsp)
    sor_step(Ap, bp, 1.63, xsorp)
    xsorsp.append(xsorp.copy())

jac_error  = [ norm(dot(A,x) - b) for x in xjacobis ]
gs_error   = [ norm(dot(A,x) - b) for x in xgss ]
gs_errorp =  [ norm(dot(Ap,x) - bp) for x in xgssp ]
sor_errorp = [ norm(dot(Ap,x) - bp) for x in xsorsp ]

graph.plot(range(steps+1), jac_error, "bo", clip_on=False)
graph.plot(range(steps+1), gs_errorp, "r+", clip_on=False)
graph.plot(range(steps+1), gs_error, "rD", clip_on=False)
graph.plot(range(steps+1), sor_errorp, "g*", clip_on=False)

graph.axis((0,steps+1,1e-17,10))
graph.xaxis.set_label_text("Schritte")
graph.yaxis.set_label_text("Fehler")

##########################################

graph = figure.add_subplot(122)

graph.set_yscale("log")

# Anzahl SOR-Schritte, nach der wir messen
steps = 20

omegas = linspace(0,1,10, endpoint=False)
omegas = concatenate((omegas, linspace(1,2,30)))
errors = []

for omega in omegas:
    xsorp = zeros(N*N)

    for step in range(steps):
        sor_step(Ap, bp, omega, xsorp)
    errors.append(norm(dot(Ap,xsorp) - bp))

graph.plot(omegas, errors, "g*", clip_on=False)
graph.axis((0,2,5e-3,20))
graph.xaxis.set_label_text(u"ω")

figure.savefig("iterative.pdf")

pyplot.show()
