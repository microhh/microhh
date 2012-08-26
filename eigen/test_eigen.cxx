#include <cstdlib>
#include <cstdio>
#include <Eigen/Dense>
using namespace Eigen;

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
    {
      ijk = 3+j*jj + k*kk;
      Map<VectorXd> vcc0(&a[ijk], isize);
      Map<VectorXd> sum (&b[ijk], isize);

      Map<VectorXd> vim3(&a[ijk-3*ii], isize);
      Map<VectorXd> vim2(&a[ijk-2*ii], isize);
      Map<VectorXd> vim1(&a[ijk-1*ii], isize);
      Map<VectorXd> vip1(&a[ijk+1*ii], isize);
      Map<VectorXd> vip2(&a[ijk+2*ii], isize);
      Map<VectorXd> vip3(&a[ijk+3*ii], isize);

      Map<VectorXd> vjm3(&a[ijk-3*jj], isize);
      Map<VectorXd> vjm2(&a[ijk-2*jj], isize);
      Map<VectorXd> vjm1(&a[ijk-1*jj], isize);
      Map<VectorXd> vjp1(&a[ijk+1*jj], isize);
      Map<VectorXd> vjp2(&a[ijk+2*jj], isize);
      Map<VectorXd> vjp3(&a[ijk+3*jj], isize);

      Map<VectorXd> vkm3(&a[ijk-3*kk], isize);
      Map<VectorXd> vkm2(&a[ijk-2*kk], isize);
      Map<VectorXd> vkm1(&a[ijk-1*kk], isize);
      Map<VectorXd> vkp1(&a[ijk+1*kk], isize);
      Map<VectorXd> vkp2(&a[ijk+2*kk], isize);
      Map<VectorXd> vkp3(&a[ijk+3*kk], isize);

      sum += c0*vim3 + c1*vim2 + c2*vim1 + c3*vcc0 + c4*vip1 + c5*vip2 + c6*vip3;
      sum += c0*vjm3 + c1*vjm2 + c2*vjm1 + c3*vcc0 + c4*vjp1 + c5*vjp2 + c6*vjp3;
      sum += c0*vkm3 + c1*vkm2 + c2*vkm1 + c3*vcc0 + c4*vkp1 + c5*vkp2 + c6*vkp3;
    }

  return 0;
}

int main()
{
  int isize = 8;
  int jsize = 8;
  int ksize = 2048;

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

