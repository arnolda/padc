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
    """Ein Newtonschritt zur Nullstellensuche."""
    p = series[-1]
    dp = (len(series) - 1) * series[-1]
    for i, coeff in reversed(list(enumerate(series[1:-1]))):
        p  =  p * xn +   coeff
        dp = dp * xn + i * coeff
    p = p * xn + series[0]
    return xn - p / dp
