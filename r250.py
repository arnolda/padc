# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
class R250:
    m = (1 << 31) - 1

    def __init__(self, init_rng):
        # Zustand mit anderem Generator initialisieren
        self.state = []
        for i in range(250):
            self.state.append(init_rng())

    def __call__(self):
        newval = (self.state[0] + self.state[250 - 103]) % self.m
        # neuen Wert anhaengen
        self.state.append(newval)
        # ... und aeltesten rauswerfen
        del self.state[0]
        return self.state[-1]
