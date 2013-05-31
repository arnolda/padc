from numpy import *
from matplotlib.pyplot import *

T = 3.0
N = 32
sigma = 0.2

def gauss(x, mu=0.0, sigma=1.0):
    return 1.0/(sigma*sqrt(2.0*pi)) * exp(-0.5 * ((x-mu)/sigma)**2)

figure()

ixs = arange(-N/2, N/2)

# create and plot Gaussian
t = linspace(-T/2, T/2, 1000)
g = gauss(t, sigma=sigma)
subplot(321)
plot(t, g, '-')

# discretize Gaussian
t = linspace(-T/2, T/2, N)
g = gauss(t, sigma=sigma)

subplot(323)
plot(g, 'o-')
xlim((0, 31))

# compute and plot FT
n = arange(-N/2,N/2)
ghat = T/N * 1./sqrt(2*pi) * fft.fft(g)
#ghat = T/N * ((-1.)**n)/sqrt(2*pi) * fft.fft(g)
subplot(325)
plot(ixs, real(ghat[ixs]), 'o-')
xlim((-16, 15))

# discretize Gaussian and shift it
g = roll(g, N/2)
subplot(324)
plot(g, 'o-')
xlim((0, 31))

# compute and plot FT
ghat = T/N * 1./sqrt(2*pi) * fft.fft(g)
subplot(326)
plot(ixs, ghat[ixs], 'o-')
xlim((-16, 15))

show()
