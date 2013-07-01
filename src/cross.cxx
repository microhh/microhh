#include <cstdio>
#include <cmath>
#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "cross.h"
#include "defines.h"

ccross::ccross(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

ccross::~ccross()
{
}

int ccross::readinifile(cinput *inputin)
{
  int n = 0;

  // optional, by default switch cross off
  n += inputin->getItem(&swcross, "cross", "swcross", "", "0");

  if(swcross == "1")
  {
    // get the list of indices at which to take cross sections
    n += inputin->getList(&jxz, "cross", "jxz", "");
    n += inputin->getList(&kxy, "cross", "kxy", "");

    // get the list of variables per type of cross
    n += inputin->getList(&simple, "cross", "simple", "");
    n += inputin->getList(&lngrad, "cross", "lngrad", "");
  }

  if(n > 0)
    return 1;

  return 0;
}

int ccross::exec(int iteration)
{
  if(swcross == "0")
    return 0;

  if(mpi->mpiid == 0) std::printf("Saving cross sections for iteration %d\n", iteration);

  // loop over chosen indices in order to make xz cross sections
  for(std::vector<int>::iterator itn = jxz.begin(); itn < jxz.end(); ++itn)
  {
    for(std::vector<std::string>::iterator itc = simple.begin(); itc < simple.end(); ++itc)
      crosssimple(fields->s[*itc]->data, fields->s["tmp1"]->data, fields->s[*itc]->name, *itn, iteration);

    for(std::vector<std::string>::iterator itc = lngrad.begin(); itc < lngrad.end(); ++itc)
      crosslngrad(fields->s[*itc]->data, fields->s["tmp1"]->data, fields->s["tmp2"]->data, grid->dzi4, fields->s[*itc]->name + "lngrad", *itn, iteration);
  }
  
  return 0;
}

int ccross::crosssimple(double * restrict data, double * restrict tmp, std::string name, int index, int iteration)
{
  // define the file name
  char filename[256];
  std::sprintf(filename, "%s.%s.%07d", name.c_str(), "xzcross", iteration);
  if(mpi->mpiid == 0) std::printf("Saving \"%s\"\n", filename);

  // save the slice with the correct index
  grid->savexzslice(data, tmp, index, filename);

  return 0;
}
 
int ccross::crosslngrad(double * restrict a, double * restrict lngrad, double * restrict tmp, double * restrict dzi4, std::string name, int index, int iteration)
{
  int ijk,ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  double dxi,dyi;
  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // calculate the log of the gradient
  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        lngrad[ijk] = std::log( 
                        std::pow( ( cg0*(ci0*a[ijk-ii3] + ci1*a[ijk-ii2] + ci2*a[ijk-ii1] + ci3*a[ijk    ])
                                  + cg1*(ci0*a[ijk-ii2] + ci1*a[ijk-ii1] + ci2*a[ijk    ] + ci3*a[ijk+ii1])
                                  + cg2*(ci0*a[ijk-ii1] + ci1*a[ijk    ] + ci2*a[ijk+ii1] + ci3*a[ijk+ii2])
                                  + cg3*(ci0*a[ijk    ] + ci1*a[ijk+ii1] + ci2*a[ijk+ii2] + ci3*a[ijk+ii3]) ) * cgi*dxi, 2.)

                      + std::pow( ( cg0*(ci0*a[ijk-jj3] + ci1*a[ijk-jj2] + ci2*a[ijk-jj1] + ci3*a[ijk    ])
                                  + cg1*(ci0*a[ijk-jj2] + ci1*a[ijk-jj1] + ci2*a[ijk    ] + ci3*a[ijk+jj1])
                                  + cg2*(ci0*a[ijk-jj1] + ci1*a[ijk    ] + ci2*a[ijk+jj1] + ci3*a[ijk+jj2])
                                  + cg3*(ci0*a[ijk    ] + ci1*a[ijk+jj1] + ci2*a[ijk+jj2] + ci3*a[ijk+jj3]) ) * cgi*dyi, 2.)

                      + std::pow( ( cg0*(ci0*a[ijk-kk3] + ci1*a[ijk-kk2] + ci2*a[ijk-kk1] + ci3*a[ijk    ])
                                  + cg1*(ci0*a[ijk-kk2] + ci1*a[ijk-kk1] + ci2*a[ijk    ] + ci3*a[ijk+kk1])
                                  + cg2*(ci0*a[ijk-kk1] + ci1*a[ijk    ] + ci2*a[ijk+kk1] + ci3*a[ijk+kk2])
                                  + cg3*(ci0*a[ijk    ] + ci1*a[ijk+kk1] + ci2*a[ijk+kk2] + ci3*a[ijk+kk3]) ) * dzi4[k], 2.) );
      }

  // define the file name
  char filename[256];
  std::sprintf(filename, "%s.%s.%07d", name.c_str(), "xzcross", iteration);
  if(mpi->mpiid == 0) std::printf("Saving \"%s\"\n", filename);

  // save the slice with the correct index
  grid->savexzslice(lngrad, tmp, index, filename);

  return 0;
}
 
