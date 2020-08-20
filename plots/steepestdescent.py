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
# Steepest Descent
#
############################################
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

def armijo_steepest_descent(f, gradf, x0, alpha=0.1, rho=0.5,
                            tolerance=1e-10, maxsteps=None):
    ls = []
    x = np.array(x0)
    path = [ x.copy() ]
    fx = f(x)
    grad = gradf(x)

    step = 0
    while np.linalg.norm(grad) > tolerance and (not maxsteps or maxsteps > step):
        step += 1

        d = -grad

        # Armjio-Schrittweite
        if alpha != None:
            lmbda = 1
            xneu = x + d
            grad2 = np.dot(grad, d)
            while  f(xneu) > fx + alpha*lmbda*grad2:
                lmbda = rho*lmbda
                xneu = x + lmbda*d
        else:
            lmbda = rho
            xneu = x + lmbda*d

        # mit dem gefundenen Schritt vorwaerts
        x = xneu
        path.append(x.copy())
        ls.append(lmbda)
        fx = f(x)
        grad = gradf(x)

        # Abbruch, wenn wir kompletten Unsinn rechnen
        if min(x) < -5 or max(x) > 5:
            break
    return np.array(path).transpose(), np.array(ls)

############################################

figure = plt.figure(figsize=(8,8))
figure.subplots_adjust(left=0.05, right=0.95)
plt.gray()

def draw(f, gradf, x0, where, method, axis, levels, limit=None):
    graph = figure.add_subplot(where)

    print("method is", method)

    if method == "armijo":
        path, lambdas = armijo_steepest_descent(f, gradf, x0, tolerance=0.01)
    else:
        path, lambdas = armijo_steepest_descent(f, gradf, x0, alpha=None,
                                                rho=0.01, tolerance=0.01,
                                                maxsteps=limit)


    print("steps taken", len(lambdas))
    print("min step was", min(lambdas))
    print("max step was", max(lambdas))

    # contour
    n=100
    rx = np.linspace(axis[0], axis[1], n)
    ry = np.linspace(axis[2], axis[3], n)
    x, y = np.meshgrid(rx, ry)
    z = np.zeros_like(x)
    for k in range(n):
        for l in range(n):
            z[k, l] = max(1e-20, f([ x[k, l], y[k, l] ]))

    cont = graph.contour(x, y, z, levels=levels,
                         norm=colors.LogNorm(min(levels), max(levels)*2),linewidths=0.5)
    graph.clabel(cont, inline=1, fontsize=10, fmt="%g")

    # target cross
    graph.plot(1, 1 , "+", markeredgecolor="black",markersize=20)

    # path
    graph.plot(path[0], path[1] , "-",color="#ffa0a0",linewidth=1)
    graph.plot(path[0], path[1] , "o",
               markeredgecolor="red",markerfacecolor="red",markersize=3)
    graph.plot(path[0,::200], path[1,::200] , "o",
               markeredgecolor="#808080",markerfacecolor="#808080",markersize=4)

    graph.axis(axis)

# Oben, quadratische Funktion
############################################

A = np.array([[70,30],[30,50]])

print(A)

def f(x):
    return np.dot(x, np.dot(A,x))

def gradf(x):
    return 2*np.dot(A,x)

draw(f, gradf, (-0.5,-1), 221, "armijo", (-1,1,-1,1), (1,10,30,50,100))
draw(f, gradf, (-0.5,-1), 222, "fixed", (-1,1,-1,1), (1,10,30,50,100))

# Unten, Rosenbrockfunktion
############################################

def rosenbrock(x):
    return (1.0-x[0])**2 + 100.0*(x[1]-x[0]**2)**2

def gradrosenbrock(x):
    return np.array((2*(200*x[0]**3 - 200*x[0]*x[1] + x[0] - 1),
                  200*(x[1]-x[0]**2)))

draw(rosenbrock, gradrosenbrock, (0,0.2), 223, "armijo", (0,1,-0.25,1), (0.2,1,5,10,15,20))
draw(rosenbrock, gradrosenbrock, (0,0.2), 224, "fixed", (0,1,-0.25,1), (0.2,1,5,10,15,20), limit=20)

figure.savefig("steepestdescent.pdf")
