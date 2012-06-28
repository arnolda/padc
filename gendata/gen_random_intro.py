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

import numpy.random as random
import os
import time

start = time.time()
random.uniform(0,1, 1000000)
end = time.time()
ri_time = end-start
print "randint -> %.1f ms" % (ri_time*1000.0)

start = time.time()
os.urandom(4000000)
end = time.time()
ur_time = end-start
print "urandom -> %.1f ms" % (ur_time*1000.0)
print "factor %.1f" % (ur_time / ri_time)

c = [ 0 for x in range(6) ]
for l in range(3):
    s = ""
    for x in random.randint(1,7,60):
        c[x-1]+=1
        s += str(x)
    print s,"\\\\"

print c
