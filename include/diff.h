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

#ifndef DIFF
#define DIFF

// forward declaration to speed up build time
class cmodel;
class cgrid;
class cfields;
class cmpi;

class cdiff
{
  public:
    cdiff(cmodel *);
    virtual ~cdiff();

    virtual int readinifile(cinput *);
    virtual int setvalues();
    virtual int execvisc();
    virtual int exec();

    virtual unsigned long gettimelim(unsigned long, double);
    virtual double getdn(double);

    double dnmax;

  protected:
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
