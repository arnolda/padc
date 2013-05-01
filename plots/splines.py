# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
# Spline-Interpolation der Rungefunktion
#
############################################
import scipy.interpolate as ip
import numpy as np
import matplotlib.pyplot as pyplot

# die zu interpolierende Funktion
def runge(x): return 1/(1+x**2)

# derivatives
def dd(f, x, step=0.001):
    "second derivative"
    return (-0.5*f(x + step) + f(x) - 0.5*f(x - step))/step**2

def d(f, x, step=0.001):
    "derivative"
    return (f(x + step) - f(x))/step

# interpolating splines
class SplineInterpolation:
    def __init__(self, x, y, leftddx = 0, rightddx = 0):
        n = len(x)
        if n < 3:
            raise Exception("need at least three support points to compute a spline")
        self.x = x.copy()
        self.y = y.copy()
        self.n = n

        # compute lambda, mu and S
        mu = np.zeros(n-1)
        mu[1:]=(x[1:-1]-x[:-2])/(x[2:]-x[:-2])

        l = np.zeros(n-1)
        l[1:]=(x[2:]-x[1:-1])/(x[2:]-x[:-2])

        S = np.zeros(n)
        S[1:-1]=6*(((y[2:]  -y[1:-1])/(x[2:]  -x[1:-1]))- \
                   ((y[1:-1]-y[:-2]) /(x[1:-1]-x[:-2]))) / (x[2:]-x[:-2])
        # second moment boundary conditions
        S[0]    = leftddx
        S[n-1]  = rightddx

        # construct coeff matrix
        tmp = np.zeros((n, n))
        tmp[0,0] = 1
        for i in range(1, n-1):
            tmp[i,i] = 2
            tmp[i,i-1] = mu[i]
            tmp[i,i+1] = l[i]
        tmp[n-1,n-1] = 1

        # compute coeff
        self.M = np.linalg.solve(tmp, S)

        # compute the dependent coeffs
        self.alpha = (self.M[1:]-self.M[:-1])/(self.x[1:]-self.x[:-1])

        self.m = (y[1:]-y[:-1])/(x[1:]-x[:-1]) - \
            1./6.*(x[1:]-x[:-1])*(2*self.M[:-1]+self.M[1:])

    def __call__(self, x):
        y = np.zeros_like(x)
        for i in range(0, len(x)):
            index=0
            while index < len(self.x) - 1 and self.x[index] <= x[i]:
                index += 1
            if index > 0: index -= 1

            y[i] = self.y[index] + \
                self.m[index]*(x[i] - self.x[index]) + \
                0.5*self.M[index]*(x[i] - self.x[index])**2 + \
                1./6.*self.alpha[index]*(x[i] - self.x[index])**3
        return np.array(y)

x = np.linspace(-5,5,11)
y = runge(x)

cubsp11 = SplineInterpolation(x, y)

x = np.linspace(-5,5,7)
y = runge(x)

linsp7 = ip.interp1d(x, y, kind = "linear", bounds_error=False)
cubsp7 = SplineInterpolation(x, y)

# Ausgabe
#################################
px = np.linspace(-5.1, 5.1, 200)

figure = pyplot.figure(figsize=(8,4))
figure.subplots_adjust(left=0.05, bottom=0.1, wspace=0.2,right=0.95, top=0.95)

# left: spline interpolation
graph = figure.add_subplot(121)
graph.plot(px, runge(px), "k", linewidth=0.5)
graph.plot(px, cubsp11(px), "r--",linewidth=1.5)
graph.plot(px, cubsp7(px), "g:",linewidth=1.5)
graph.plot(px, linsp7(px), "b-.")
graph.axis([-5.1,5.1,0,1.1])

# right: derivatives
graph = figure.add_subplot(122)
graph.plot(px, dd(runge, px), "k",linewidth=0.5)
graph.plot(px, d(cubsp11, px), "r:",linewidth=1.5)
graph.plot(px, dd(cubsp11, px), "b",linewidth=1.5)
graph.axis([-5.1,5.1,-0.7,0.7])

figure.savefig("splines.pdf")
