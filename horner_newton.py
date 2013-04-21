# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
def poly_newton_step(series, xn):
    n = len(series) - 1
    p  = series[n]
    dp = n*series[n]
    for i in range(n - 1,0,-1):
        p  =  p*xn +   series[i]
        dp = dp*xn + i*series[i]
    p = p*xn +   series[0]
    return xn - p/dp
