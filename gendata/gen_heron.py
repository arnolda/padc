# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Demonstration Heron-Verfahren
######################################

import math
import numpy as np

a = 2
res = math.sqrt(a)


def digits(x):
    return int(math.floor(-math.log(abs(res - x)) / math.log(10)))


def heron(x):
    return 0.5 * (x + a / x)


x = 1
for step in range(8):
    print(f"{step} & {x:.15f} & {digits(x)} \\\\")
    x = heron(x)

# Alternativ: Polynom-Interpolation
############################################

import scipy.interpolate as inter

# Chebyshev
steps = 3 * 5 // 2

print(f"Polynome mit Koeffizienten: {steps}")

x = np.zeros(steps)
for i in range(steps):
    x[i] = 2.5 + 2.5 * math.cos((2. * i + 1) / (2. * steps) * np.pi)

val = inter.lagrange(x, np.sqrt(x))(2)
print(f"mit Lagrange und Chebyshev auf [0,5]: {val} {digits(val)}")

# Alternativ: Taylor um 1
############################################

val = 1
binom = 1
for n in range(1, steps + 1):
    binom *= (0.5 - n + 1) / n
    val += binom

print(f"mit Taylor um 1 {val} {digits(val)}")
