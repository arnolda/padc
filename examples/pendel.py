# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugänglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
#
"""
Simulation der Bahn eines Pendels
mit Hilfe des Velocity-Verlet-Integrators.
"""
import argparse
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


class Pendel:
    """Kraft und Energie des Pendels.

    Arguments:
        l: Laenge des Pendelarms.
        g: Erdbeschleunigung.

    """

    def __init__(self, g, l):
        self.g = g
        self.l = l

    def kraft(self, a):
        "Berechnet die Kraft auf das Pendel"
        return -self.g / self.l * np.sin(a)

    def energie(self, a, da):
        """Berechnet die potentielle Energie"""
        return 0.5 * (self.l * da)**2 + self.g * (self.l - self.l * np.cos(a))


class Propagator:
    """Propagiert das System in erster Ordnung.

    Arguments:
        dt: Zeitschritt
        start_a: Startposition
        start_da: Startwinkelgeschwindigkeit

    """

    def __init__(self, system, dt, start_a, start_da):
        self.system = system
        self.dt = dt
        self.a = start_a
        self.da = start_da
        self.t = 0

    def step(self):
        """Einen Zeitschritt vorwärts"""
        self.da += self.system.kraft(self.a) * self.dt
        self.a += self.da * self.dt

        self.t += self.dt

    def propagate(self, n):
        """n Schritte vorwärts"""
        for _ in range(n):
            self.step()

    def energie(self):
        """Energie im System"""
        return self.system.energie(self.a, self.da)


class VVPropagator(Propagator):
    "Propagiert das System mit Hilfe des VV"

    def step(self):
        """Einen Zeitschritt vorwärts"""
        self.da += 0.5 * self.system.kraft(self.a) * self.dt
        self.a += self.da * self.dt
        self.da += 0.5 * self.system.kraft(self.a) * self.dt
        self.t += self.dt


class Recorder:
    """Aufnahme der Observablen in Ringpuffern.

    Arguments:
        memory: Dauer (Breite) der zyklischen Graphen in Zeiteinheiten
    """

    def __init__(self, memory):
        self.memory = memory
        self.a = []
        self.v = []
        self.E = []
        self.t = []

    def save(self, t, a, da, E):
        """Werte im Ringpuffer sichern"""
        self.a.append(a)
        self.v.append(da)
        self.E.append(E)
        self.t.append(t)

        if self.t[-1] - self.t[0] > self.memory:
            del self.t[0], self.E[0], self.a[0], self.v[0]


def draw_pendulum(kurven, propagator, system):
    """Ausgabe des Pendels."""
    kurve_w, kurve_arm, *_ = kurven

    # Skizze der Pendelbewegung
    kurve_w.set_data(
        [system.l * np.sin(propagator.a)],
        [system.l - system.l * np.cos(propagator.a)])
    kurve_arm.set_data(
        [0, system.l * np.sin(propagator.a)],
        [system.l, system.l * (1 - np.cos(propagator.a))])


def draw_observables(kurven, recorder):
    """Ausgabe der Observablen in Graphen."""
    _, _, kurve_a, graph_a, kurve_v, graph_v, kurve_E, graph_E = kurven
    t_min, t_max = recorder.t[0], recorder.t[-1]

    kurve_a.set_data(recorder.t, recorder.a)
    graph_a.axis((t_min, t_max, np.min(recorder.a), np.max(recorder.a)))

    kurve_v.set_data(recorder.t, recorder.v)
    graph_v.axis((t_min, t_max, np.min(recorder.v), np.max(recorder.v)))
    # die Energie ist praktisch konstant, daher den Bereich etwas vergroessern
    E_min = 0.99 * np.min(recorder.E)
    E_max = 1.01 * np.max(recorder.E)
    kurve_E.set_data(recorder.t, recorder.E)
    graph_E.axis((t_min, t_max, E_min, E_max))


def init_curves(abb, system):
    """Initialisiert ein Graphikfenster mit 4 Subgraphen."""

    # Pendelskizze
    graph_pos = abb.add_subplot(221, title="Pendel")
    kurve_arm, kurve_w = graph_pos.plot([], [], "b", [], [], "ro")
    graph_pos.axis((-1.1 * system.l, 1.1 * system.l,
                    -0.1 * system.l, 2.1 * system.l))

    graph_v = abb.add_subplot(222, title="Winkelgeschwindigkeit")
    kurve_v, = graph_v.plot([], [])
    graph_a = abb.add_subplot(223, title="Winkel")
    kurve_a, = graph_a.plot([], [])
    graph_E = abb.add_subplot(224, title="Energie")
    kurve_E, = graph_E.plot([], [])
    return (kurve_w, kurve_arm,
            kurve_a, graph_a,
            kurve_v, graph_v,
            kurve_E, graph_E)


def animate(_frame, propagator, steps, recorder, kurven):
    """
    Die Hauptroutine, die staendig die Pendelposition in der Zeit propagiert.
    """
    # einige Schritte propagieren
    propagator.propagate(steps)

    recorder.save(propagator.t, propagator.a,
                  propagator.da, propagator.energie())

    draw_pendulum(kurven, propagator, propagator.system)
    draw_observables(kurven, recorder)


def main():
    """Hauptroutine"""
    parser = argparse.ArgumentParser()

    parser.add_argument("--dt", dest="dt",
                        help="Integrationszeitschritt", metavar="Zeitschritt",
                        type=float, default=0.01)
    parser.add_argument("--verlet", dest="method",
                        help="Velocity-Verlet-Integrator benutzen",
                        action="store_const", const="verlet",
                        default="simple")
    parser.add_argument("--start_a", dest="start_a",
                        help="Anfangsauslenkung", metavar="Winkel",
                        type=float, default=0.1)
    parser.add_argument("--start_da", dest="start_da",
                        help="Anfangsgeschwindigkeit",
                        metavar="Geschwindigkeit",
                        type=float, default=0)
    parser.add_argument("--l", dest="l",
                        help=u"Länge des Pendels", metavar=u"Länge",
                        type=float, default=1)
    parser.add_argument("--limit", dest="limit",
                        help=u"Zeitlimit in 0.1s", metavar=u"Limit",
                        type=int)

    args = parser.parse_args()

    system = Pendel(l=args.l, g=9.81)

    propagator_class = Propagator if args.method == "simple" else VVPropagator
    propagator = propagator_class(system, args.dt, args.start_a, args.start_da)

    recorder = Recorder(20)

    abb = plt.figure()

    kurven = init_curves(abb, system)

    anim = FuncAnimation(
        abb, animate, fargs=[propagator, 10, recorder, kurven],
        frames=None, interval=50, repeat=False)
    plt.show(block=False)
    for _ in range(args.limit) if args.limit else it.count():
        plt.pause(.1)
    return anim


if __name__ == "__main__":
    main()
