/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BOUNDARY_SURFACE
#define BOUNDARY_SURFACE

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "boundary.h"
#include "thermo.h"

class cboundary_surface : public cboundary
{
  public:
    cboundary_surface(cgrid *, cfields *, cmpi *);
    ~cboundary_surface();

    int readinifile(cinput *);
    int init();
    int setvalues();
    int exec();

    int setdepends(cthermo *);

    int save(int);
    int load(int);

    double *obuk;
    double *ustar;

  private:
    // surface scheme
    int bcvalues();
    int stability(double *, double *, double *,
                  double *, double *, double *,
                  double *, double *, double *,
                  double *, double *);
    int stability_neutral(double *, double *,
                          double *, double *,
                          double *, double *,
                          double *, double *);
    int surfm(double *, double *,
              double *, double *, double *, double *,
              double *, double *, double *, double *,
              double, int);
    int surfs(double *, double *, double *,
              double *, double *, double *,
              double, int);
    double calcobuk_noslip_flux     (double, double, double, double);
    double calcobuk_noslip_dirichlet(double, double, double, double);
    inline double fm(double, double, double);
    inline double fh(double, double, double);
    inline double psim(double);
    inline double psih(double);
    inline double phim(double);
    inline double phih(double);
    double ustarin;
    double z0m;
    double z0h;

    bool allocated;

    typedef std::map<std::string, int> bcbotmap;
    int surfmbcbot;
    bcbotmap surfsbcbot;

    cthermo *thermo;
};
#endif
