# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
#
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
#
# Test der Zufallszahlengeneratoren
##############################################

import sys
# da liegen die Methoden, da sie Teil des Skripts sind
sys.path.append("..")


def test_rng(name, rng, expected):
    print(f"testing {name}")

    for val in expected:
        check = rng()
        assert check == val, f"mismatch {check} != {val}"


# Test MINSTD
##################################
import minstd

rng = minstd.MinStd(1)
expected = [
    16807,
    282475249,
    1622650073,
    984943658,
    1144108930,
    470211272,
    101027544,
    1457850878,
    1458777923,
    2007237709]

test_rng("minstd", rng, expected)

# Test R250
##################################
import r250

rng = r250.R250(minstd.MinStd(1))
expected = [
    1213127486,
    823912584,
    524243432,
    1981441630,
    1060561381,
    740860367,
    530002863,
    2143434332,
    662639499,
    131866351
]

test_rng("r250", rng, expected)

# Test rand
##################################
import rand

rng = rand.Rand(1)
expected = [
    1103527590,
    377401575,
    662824084,
    1147902781,
    2035015474,
    368800899,
    1508029952,
    486256185,
    1062517886,
    267834847
]

test_rng("rand", rng, expected)

# Zeigt, dass bei rand das niedrigste Bit periodisch wechselt
##################################

rng = rand.Rand(1)
print([rng() % 2 for i in range(10)])
