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
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(123)

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
h=0.1
# Matrix
Ap=np.zeros((N*N, N*N))

eqn=0
for y in range(N):
    for x in range(N):
        if eqn > 0:
            Ap[eqn, linindex(x,y)]   = -4/h**2
            Ap[eqn, linindex(x+1,y)] = 1/h**2
            Ap[eqn, linindex(x-1,y)] = 1/h**2
            Ap[eqn, linindex(x,y+1)] = 1/h**2
            Ap[eqn, linindex(x,y-1)] = 1/h**2
        else:
            # Normierung
            Ap[0, 0] = 1
        eqn += 1

# zufaelliger Zielvektor/Ladungsdichte
bp = np.random.normal(0,1,N*N)

# zum Vergleich, strikt diagonal dominant
##############################################

A = np.random.normal(0,1,n*n)
A = A.reshape((n,n))
# strikt diagonal dominant machen
for i in range(n):
    A[i,i] = sum([abs(A[i,k]) for k in range(n)]) - abs(A[i,i])/0.9

# zufaelliger Zielvektor/Ladungsdichte
b = np.random.normal(0,1,n)

# Verfahren
##############################################

def jacobi_step(A, b, x):
    xnew = np.zeros_like(x)
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

figure = plt.figure(figsize=(8,4))
figure.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9)

##########################################

graph = figure.add_subplot(121)

graph.set_yscale("log")

# Anzahl Jacobi/etc Schritte
steps  = 25
isteps = 4

# gutartige Matrix A
xjacobi  = np.zeros(n)
xjacobis  = [ xjacobi ]
xgs      = np.zeros(n)
xgss      = [ xgs.copy() ]
# nicht strikt diagdom Ap
xjacobip = np.zeros(N*N)
xjacobisp = [ xjacobip ]
xgsp     = np.zeros(N*N)
xgssp     = [ xgsp ]
xsorp    = np.zeros(N*N)
xsorsp    = [ xsorp.copy() ]

for step in range(0, steps):
    for istep in range(isteps):
        xjacobi = jacobi_step(A, b, xjacobi)
        sor_step(A, b, 1, xgs)
        
        xjacobip = jacobi_step(Ap, bp, xjacobip)
        sor_step(Ap, bp, 1, xgsp)
        sor_step(Ap, bp, 1.66, xsorp)

    xjacobis.append(xjacobi.copy())
    xgss.append(xgs.copy())
    
    xjacobisp.append(xjacobip.copy())
    xgssp.append(xgsp.copy())
    xsorsp.append(xsorp.copy())

jac_error  = [ np.linalg.norm(np.dot(A,x) - b) for x in xjacobis ]
gs_error   = [ np.linalg.norm(np.dot(A,x) - b) for x in xgss ]
jac_errorp = [ np.linalg.norm(np.dot(Ap,x) - bp) for x in xjacobisp ]
gs_errorp  = [ np.linalg.norm(np.dot(Ap,x) - bp) for x in xgssp ]
sor_errorp = [ np.linalg.norm(np.dot(Ap,x) - bp) for x in xsorsp ]

graph.plot(range(0, isteps*(steps+1),isteps), jac_error, "bo", clip_on=False)
graph.plot(range(0, isteps*(steps+1),isteps), gs_error, "rD", clip_on=False)
graph.plot(range(0, isteps*(steps+1),isteps), jac_errorp, "b.", clip_on=False)
graph.plot(range(0, isteps*(steps+1),isteps), gs_errorp, "r+", clip_on=False)
graph.plot(range(0, isteps*(steps+1),isteps), sor_errorp, "g*", clip_on=False)

graph.xaxis.set_label_text("Schritte")
graph.yaxis.set_label_text("Fehler")

##########################################

graph = figure.add_subplot(122)

graph.set_yscale("log")

# Anzahl SOR-Schritte, nach der wir messen
steps = 20

omegas = np.linspace(0,1,10, endpoint=False)
omegas = np.concatenate((omegas, np.linspace(1,2,30)))
errors = []

for omega in omegas:
    xsorp = np.zeros(N*N)

    for step in range(steps):
        sor_step(Ap, bp, omega, xsorp)
    errors.append(np.linalg.norm(np.dot(Ap,xsorp) - bp))

graph.plot(omegas, errors, "g*", clip_on=False)
graph.xaxis.set_label_text(u"ω")

figure.savefig("iterative.pdf")

print("optimal omega:", omegas[np.argmin(errors)])
