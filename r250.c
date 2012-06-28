/* Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
   Axel Arnold, Universitaet Stuttgart.

   Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
   Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
   zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
   http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
   schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
   View, California, 94041, USA. */
/* Position des am weitesten verzoegerten Beitrags x_{n-250}
   Dieser wird im Ringspeicher durch den neuen Wert ersetzt,
   da nicht mehr gebraucht */
int position = 0;
// Der Ringspeicher
unsigned int state[250];

void r250_seed(unsigned int seed) {
  // etwa minstd zur Initialisierung
  minstd_seed(seed);
  for (int i = 0; i < 250; ++i)
    state[i] = minstd();
  position = 0;
}

unsigned int r250() {
  const unsigned int m = (1u << 31) - 1;
  unsigned int newval =
    (state[position] + state[(position + 250 - 103) % 250]) % m;
  state[position] = newval;
  // Ringspeicher weiterschieben
  position = (position + 1) % 250;
  return newval;
}
