/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
class cmaster;

class cdiff
{
  public:
    cdiff(cmodel *, cinput *);
    virtual ~cdiff();
    static cdiff* factory(cmaster *, cinput *, cmodel *, const std::string); ///< Factory function for diff class generation.

    virtual void setvalues();
    virtual int execvisc();
    virtual int exec();

    std::string getname();
    virtual unsigned long gettimelim(unsigned long, double);
    virtual double getdn(double);

    double dnmax;

    // GPU functions and variables
    virtual int prepareDevice(); 

  protected:
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
    cmaster *master;

    std::string swdiff;
};
#endif
