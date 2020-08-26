# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Demonstration Division
#########################

import math

a = 2.0
res = 1.0 / a


def digits(x):
    return int(math.floor(-math.log(abs(res - x)) / math.log(10)))


def div(x):
    return x * (2 - a * x)


x = 0.9
for step in range(6):
    print(f"{step} & {x:.15f} & {digits(x)} \\\\")
    x = div(x)
