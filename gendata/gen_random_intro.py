# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.

# Ausmessen der Geschwindigkeit von urandom und randint
#######################################################

import os
import time
import numpy as np

start = time.time()
np.random.uniform(0, 1, 1000000)
end = time.time()
ri_time = end - start
print(f"randint -> {ri_time*1000.0:.1f} ms")

start = time.time()
os.urandom(4000000)
end = time.time()
ur_time = end - start
print(f"urandom -> {ur_time*1000.0:.1f} ms")
print(f"factor {ur_time / ri_time:.1f}")

c = [0 for x in range(6)]
for l in range(3):
    s = ""
    for x in np.random.randint(1, 7, 60):
        c[x - 1] += 1
        s += str(x)
    print(s, "\\\\")
print(c)
