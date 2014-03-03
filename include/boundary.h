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
class cmaster;
class cmodel;
class cgrid;
class cfields;

/**
 * Structure containing the boundary options and values per 3d field.
 */
struct field3dbc
{
  double bot; ///< Value of the bottom boundary.
  double top; ///< Value of the top boundary.
  int bcbot;  ///< Switch for the bottom boundary.
  int bctop;  ///< Switch for the top boundary.
};

/**
 * Base class for the advection scheme.
 * This class handles the case when advection is turned off. Derived classes are
 * implemented that handle different advection schemes.
 */
class cboundary
{
  public:
    cboundary(cmodel *);  ///< Constuctor of the boundary class.
    virtual ~cboundary(); ///< Destructor of the boundary class.

    virtual int readinifile(cinput *); ///< Process the data from the input file.
    virtual int init();                ///< Initialize the fields.
    virtual int setvalues();           ///< Set all 2d fields to the prober BC value.

    virtual int save(int); ///< Save boundary conditions related fields for restarts.
    virtual int load(int); ///< Load boundary conditions related fields for restarts.

    int exec();              ///< Update the boundary conditions.
    virtual int execcross(); ///< Execute cross sections of surface

  protected:
    cmaster *master; ///< Pointer to master class.
    cmodel  *model;  ///< Pointer to model class.
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.

    int mbcbot;
    int mbctop;

    typedef std::map<std::string, field3dbc *> bcmap;
    bcmap sbc;

    int processbcs(cinput *); ///< Process the boundary condition settings from the ini file.
    int setbc(double *, double *, double *, int, double, double, double); ///< Set the values for the boundary fields.

  private:
    virtual int bcvalues(); ///< Update the boundary values.


    int setgcbot_2nd(double *, double *, int, double *, double *); ///< Set the bottom ghost cells with 2nd-order accuracy.
    int setgctop_2nd(double *, double *, int, double *, double *); ///< Set the top ghost cells with 2nd-order accuracy.
    int setgcbot_4th(double *, double *, int, double *, double *); ///< Set the bottom ghost cells with 4th-order accuracy.
    int setgctop_4th(double *, double *, int, double *, double *); ///< Set the top ghost cells with 4th-order accuracy.

    int setgcbotw_4th(double *); ///< Set the bottom ghost cells for the vertical velocity with 4th order accuracy.
    int setgctopw_4th(double *); ///< Set the top ghost cells for the vertical velocity with 4th order accuracy.

    inline double grad4x(const double, const double, const double, const double); ///< Calculate a 4th order gradient.
};
#endif
