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

#ifndef BOUNDARY
#define BOUNDARY

// forward declarations to speed up build time
class cmodel;
class cgrid;
class cfields;
class cmpi;

struct field3dbc
{
  double bot;
  double top;
  int bcbot;
  int bctop;
};

class cboundary
{
  public:
    cboundary(cmodel *);
    virtual ~cboundary();

    virtual int readinifile(cinput *);
    virtual int init();
    virtual int setvalues();

    virtual int save(int);
    virtual int load(int);

    int exec();

  protected:
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    std::string swspatialorder;

    int mbcbot;
    int mbctop;

    typedef std::map<std::string, field3dbc *> bcmap;
    bcmap sbc;

    int setbc(double *, double *, double *, int, double, double, double);
    
  private:
    virtual int bcvalues();

    int setgcbot_2nd(double *, double *, int, double *, double *);
    int setgctop_2nd(double *, double *, int, double *, double *);
    int setgcbot_4th(double *, double *, int, double *, double *);
    int setgctop_4th(double *, double *, int, double *, double *);

    int setgcbotw_4th(double *);
    int setgctopw_4th(double *);

    inline double grad4x(const double, const double, const double, const double);
};
#endif
