/* Added by Jack Poulson to ease the integration with Elemental */
void
pmrrr_dscal( int* n, double* alpha, double* x, int* incx )
{
  int i;
  int stride = *incx;
  int height = *n;
  double scale = *alpha;

  for( i=0; i<height; ++i )
    x[i*stride] *= scale;
}

