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

#ifndef BOUNDARY
#define BOUNDARY

// forward declarations to speed up build time
class Master;
class Model;
class Grid;
class Fields;
struct Mask;

/**
 * Base class for the advection scheme.
 * This class handles the case when advection is turned off. Derived classes are
 * implemented that handle different advection schemes.
 */
class Boundary
{
  public:
    Boundary(Model *, Input *); ///< Constuctor of the boundary class.
    virtual ~Boundary();          ///< Destructor of the boundary class.
    static Boundary* factory(Master *, Input *, Model *); ///< Factory function for boundary class generation.

    virtual void init(Input *);   ///< Initialize the fields.
    virtual void create(Input *); ///< Create the fields.
    virtual void setTimeDep();
    virtual void setValues();      ///< Set all 2d fields to the prober BC value.

    virtual void save(int); ///< Save boundary conditions related fields for restarts.
    virtual void load(int); ///< Load boundary conditions related fields for restarts.

    virtual void exec();              ///< Update the boundary conditions.
    virtual void execCross(); ///< Execute cross sections of surface
    virtual void execStats(Mask *); ///< Execute statistics of surface

    enum BoundaryType {DirichletType, NeumannType, FluxType, UstarType};

    // GPU functions and variables
    virtual void prepareDevice(); 
    virtual void forwardDevice(); 
    virtual void backwardDevice(); 

  protected:
    Master *master; ///< Pointer to master class.
    Model  *model;  ///< Pointer to model class.
    Grid   *grid;   ///< Pointer to grid class.
    Fields *fields; ///< Pointer to fields class.

    BoundaryType mbcbot;
    BoundaryType mbctop;

    /**
     * Structure containing the boundary options and values per 3d field.
     */
    struct field3dbc
    {
      double bot; ///< Value of the bottom boundary.
      double top; ///< Value of the top boundary.
      BoundaryType bcbot; ///< Switch for the bottom boundary.
      BoundaryType bctop; ///< Switch for the top boundary.
    };

    typedef std::map<std::string, field3dbc *> bcmap;
    bcmap sbc;

    // time dependent variables
    std::string swtimedep;
    std::vector<double> timedeptime;
    std::vector<std::string> timedeplist;
    std::map<std::string, double *> timedepdata;

    void processbcs(Input *); ///< Process the boundary condition settings from the ini file.
    void processtimedep(Input *); ///< Process the time dependent settings from the ini file.
    void setbc(double *, double *, double *, BoundaryType, double, double, double); ///< Set the values for the boundary fields.

    // GPU functions and variables
    void setbc_g(double *, double *, double *, BoundaryType, double, double, double); ///< Set the values for the boundary fields.

  private:
    virtual void bcvalues(); ///< Update the boundary values.

    void setgcbot_2nd(double *, double *, BoundaryType, double *, double *); ///< Set the bottom ghost cells with 2nd-order accuracy.
    void setgctop_2nd(double *, double *, BoundaryType, double *, double *); ///< Set the top ghost cells with 2nd-order accuracy.
    void setgcbot_4th(double *, double *, BoundaryType, double *, double *); ///< Set the bottom ghost cells with 4th-order accuracy.
    void setgctop_4th(double *, double *, BoundaryType, double *, double *); ///< Set the top ghost cells with 4th-order accuracy.

    void setgcbotw_4th(double *); ///< Set the bottom ghost cells for the vertical velocity with 4th order accuracy.
    void setgctopw_4th(double *); ///< Set the top ghost cells for the vertical velocity with 4th order accuracy.

    inline double grad4x(const double, const double, const double, const double); ///< Calculate a 4th order gradient.

  protected:
    static const double noVelocity = 0.;
    static const double noOffset = 0.;
};
#endif
