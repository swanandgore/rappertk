
#include <math.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#include "RanGen.h"

long RanGen::seed = -1973;

RanGen::RanGen() { ran2(&seed); }

void RanGen::seedme(int r) { seed = r ; ran2(&seed); }

double RanGen::ran_uni_0_1() { return ran2(&seed); }

RanGen & RanGen::instance() {
    static RanGen rangen;
    return rangen;
}

//  adapted from Rapper
//
//  FROM NUMERICAL RECIPES IN C
//
//  Long period (>2 E18) random number generator of
//  L'Ecuyer with Bays-Durham suffle and added
//  safeguards. Returns a uniform random deviate between
//  0.0 and 1.0 (exclusive of the endpoint values).
//  Call with idum a negative integer to initialise;
//  thereafter, do not alter idum between successive
//  deviates in a sequence. RNMX should approximate the
//  largest floating value that is less than 1.
double RanGen::ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum = 1;
    else *idum = -(*idum);
    idum2 = (*idum);
    for (j = NTAB+7; j >= 0; j--) {
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0) *idum += IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp = AM * iy) > RNMX) return RNMX;
  else return temp;
} // ran2

//double ran01() {
  //  return rand() / (RAND_MAX+0.);
  //  return RanGen::instance().ran_uni_0_1();
//}

#ifdef DONT_DEFINE_ME
#include <iostream>
using namespace std;

main() {
    RanGen & r = RanGen::instance();
    for(int i=0; i < 1000; i++)
        cout << r.ran01() << endl;
}
#endif // DONT_DEFINE_ME
