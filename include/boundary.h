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

class Master;
class Model;
class Input;
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
    virtual ~Boundary();        ///< Destructor of the boundary class.

    static Boundary* factory(Master *, Input *, Model *); ///< Factory function for boundary class generation.

    virtual void init(Input *);   ///< Initialize the fields.
    virtual void create(Input *); ///< Create the fields.
    virtual void updateTimeDep(); ///< Update the time dependent parameters.
    virtual void setValues();     ///< Set all 2d fields to the prober BC value.

    virtual void save(int); ///< Save boundary conditions related fields for restarts.
    virtual void load(int); ///< Load boundary conditions related fields for restarts.

    virtual void exec(); ///< Update the boundary conditions.

    virtual void execStats(Mask *); ///< Execute statistics of surface
    virtual void execCross();       ///< Execute cross sections of surface

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
    struct Field3dBc
    {
      double bot; ///< Value of the bottom boundary.
      double top; ///< Value of the top boundary.
      BoundaryType bcbot; ///< Switch for the bottom boundary.
      BoundaryType bctop; ///< Switch for the top boundary.
    };

    typedef std::map<std::string, Field3dBc *> BcMap;
    BcMap sbc;

    // time dependent variables
    std::string swtimedep;
    std::vector<double> timedeptime;
    std::vector<std::string> timedeplist;
    std::map<std::string, double *> timedepdata;

    void processBcs(Input *);     ///< Process the boundary condition settings from the ini file.
    void processTimeDep(Input *); ///< Process the time dependent settings from the ini file.

    void setBc(double *, double *, double *, BoundaryType, double, double, double); ///< Set the values for the boundary fields.

    // GPU functions and variables
    void setBc_g(double *, double *, double *, BoundaryType, double, double, double); ///< Set the values for the boundary fields.

  private:
    virtual void updateBcs(); ///< Update the boundary values.

    void calcGhostCellsBot_2nd(double *, double *, BoundaryType, double *, double *); ///< Calculate the bottom ghost cells with 2nd-order accuracy.
    void calcGhostCellsTop_2nd(double *, double *, BoundaryType, double *, double *); ///< Calculate the top ghost cells with 2nd-order accuracy.
    void calcGhostCellsBot_4th(double *, double *, BoundaryType, double *, double *); ///< Calculate the bottom ghost cells with 4th-order accuracy.
    void calcGhostCellsTop_4th(double *, double *, BoundaryType, double *, double *); ///< Calculate the top ghost cells with 4th-order accuracy.

    void calcGhostCellsBotw_4th(double *); ///< Calculate the bottom ghost cells for the vertical velocity with 4th order accuracy.
    void calcGhostCellsTopw_4th(double *); ///< Calculate the top ghost cells for the vertical velocity with 4th order accuracy.

    inline double grad4x(const double, const double, const double, const double); ///< Calculate a 4th order gradient.

    /*
  protected:
    static const double noVelocity = 0.;
    static const double noOffset = 0.;
    */
};
#endif
