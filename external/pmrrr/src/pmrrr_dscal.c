#include <stdlib.h>
#include "global.h"

/* non-optimized, non-threaded DSCAL replacement */
void pmrrr_dscal(int *n, double *alpha, double *restrict x, int *incx)
{
  int i;
  int stride = *incx;
  int size   = *n;
  double s   = *alpha;

  for (i=0; i<size; i++)
    x[i * stride] *= s; 

  return;
}
