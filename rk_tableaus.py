# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
import numpy as np

# Butchertableau fuer das explizite Eulerverfahren
euler = { 'c': np.array((0,)),
          'A': np.array(((0,),)),
          'b': np.array((1,))}

# Butchertableau fuer das klassische Runge-Kutta-Verfahren
rk_klassisch = { 'c': np.array((0, 0.5, 0.5, 1)),
                 'A': np.array(((0,     0, 0),
                                (0.5,   0, 0),
                                (0,   0.5, 0),
                                (0,     0, 1))),
                 'b': np.array((1. / 6, 1. / 3, 1. / 3, 1. / 6))}
