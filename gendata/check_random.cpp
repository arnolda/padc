/* Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
   Axel Arnold, Universitaet Stuttgart.

   Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
   Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
   zugaenglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
   http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
   schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
   View, California, 94041, USA.

   Austesten der RNGs-Bespielcodes
*/
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

/* Namespaces, um die C-RNGs aus der Vorlesung getrennt zu halten
 */
namespace Rand {
  #include "../rand.c"
}

namespace Minstd {
  #include "../minstd.c"
}

namespace R250 {
  void minstd_seed(int _s) { Minstd::minstd_seed(_s); }
  int minstd()             { return Minstd::minstd(); }
  #include "../r250.c"
}

void test_rng(const char *name, const gsl_rng_type *T,
	      unsigned int seed, unsigned int (* myrng)())
{
  printf("testing %s\n", name);
  gsl_rng *r = gsl_rng_alloc(T);
  gsl_rng_set (r, seed);
  for (int i = 0; i < 10; i++) {
    int a = myrng(), b = gsl_rng_get(r);
    if (a != b) {
      printf("%d-%d  ", a, b);
    }
  }
  printf("\n");
  
  gsl_rng_free (r);
}

void run_rng(const char *name, unsigned int (* myrng)())
{
  printf("running %s\n", name);
  for (int i = 0; i < 10; i++) {
    printf("%d  ", myrng());
  }
  printf("\n");
}

int main()
{
  gsl_rng_env_setup();
     
  // show that rand48 is 2^17-periodic in the lowest bit
  for (int i = 0; i < 10; i++) {
    for (int i = 0; i < (1<<17)-1; i++) {
      lrand48();
    }
    printf("%ld ", lrand48() % 2);
  }
  printf("\n");

  // standard rand
  test_rng("rand", gsl_rng_rand, Rand::state, Rand::rand);

  // minstd
  test_rng("minstd", gsl_rng_minstd, Minstd::state, Minstd::minstd);

  // r250
  R250::r250_seed(1);
  run_rng("r250", R250::r250);

  return 0;
}
