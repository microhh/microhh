#include <cstdlib>
#include <cstdio>

int func(double * __restrict__ a, double * __restrict__ b, const int isize, const int jsize, const int ksize)
{
  const int c0 = 1.;
  const int c1 = 2.;
  const int c2 = 3.;
  const int c3 = 4.;
  const int c4 = 5.;
  const int c5 = 6.;
  const int c6 = 7.;

  int ijk,ii,jj,kk;
  ii = 1;
  jj = isize+6;
  kk = (isize+6)*(jsize+6);

  for(int k=3; k<ksize+3; k++)
    for(int j=3; j<jsize+3; j++)
      for(int i=3; i<isize+3; i++)
      {
        ijk = i + j*jj + k*kk;
        b[ijk] += c0*a[ijk-3*ii] + c1*a[ijk-2*ii] + c2*a[ijk-1*ii] + c3*a[ijk] + c4*a[ijk+1*ii] + c5*a[ijk+2*ii] + c6*a[ijk+3*ii]
               +  c0*a[ijk-3*jj] + c1*a[ijk-2*jj] + c2*a[ijk-1*jj] + c3*a[ijk] + c4*a[ijk+1*jj] + c5*a[ijk+2*jj] + c6*a[ijk+3*jj]
               +  c0*a[ijk-3*kk] + c1*a[ijk-2*kk] + c2*a[ijk-1*kk] + c3*a[ijk] + c4*a[ijk+1*kk] + c5*a[ijk+2*kk] + c6*a[ijk+3*kk];
      }

  return 0;
}

int main()
{
  int isize = 2048;
  int jsize = 8;
  int ksize = 8;

  int nsize = (isize+6)*(jsize+6)*(ksize+6);

  double *a, *b;
  a = new double[(isize+6)*(jsize+6)*(ksize*6)];
  b = new double[(isize+6)*(jsize+6)*(ksize*6)];

  std::srand(0);
  for(int n=0; n<nsize; n++)
  {
    a[n] = (double)(std::rand());
    b[n] = 0.;
  }

  for(int t=0; t<200; t++)
    func(a, b, isize, jsize, ksize);

  // for(int n=0; n<nsize; n++)
  //   std::printf("%d, %E\n", n, b[n]);

  delete[] a;
  delete[] b;
}

