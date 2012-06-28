/* Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
   Axel Arnold, Universitaet Stuttgart.

   Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
   Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
   zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
   http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
   schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
   View, California, 94041, USA.*/
int state = 1;

void minstd_seed(unsigned int seed)
{
  state = seed;
}

unsigned int minstd()
{
  const int m = (1u << 31) - 1, a = 16807;
  state = ((long int)state)*a % m;
  return state;
}
