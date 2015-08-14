/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
 * Base class for the boundary scheme.
 * This class handles the case when the boundary is turned off. Derived classes are
 * implemented that handle different boundary schemes.
 */
class Boundary
{
    public:
        Boundary(Model*, Input*); ///< Constuctor of the boundary class.
        virtual ~Boundary();      ///< Destructor of the boundary class.

        static Boundary* factory(Master*, Input*, Model*); ///< Factory function for boundary class generation.

        virtual void init(Input*);   ///< Initialize the fields.
        virtual void create(Input*); ///< Create the fields.

        virtual void update_time_dependent(); ///< Update the time dependent parameters.

        virtual void set_values(); ///< Set all 2d fields to the prober BC value.

        virtual void exec(); ///< Update the boundary conditions.
        virtual void set_ghost_cells_w_diff(bool); ///< Update the boundary conditions.

        virtual void exec_stats(Mask*); ///< Execute statistics of surface
        virtual void exec_cross();       ///< Execute cross sections of surface

        enum Boundary_type {Dirichlet_type, Neumann_type, Flux_type, Ustar_type};

        // GPU functions and variables
        virtual void prepare_device();
        virtual void forward_device();
        virtual void backward_device();

    protected:
        Master* master; ///< Pointer to master class.
        Model*  model;  ///< Pointer to model class.
        Grid*   grid;   ///< Pointer to grid class.
        Fields* fields; ///< Pointer to fields class.

        Boundary_type mbcbot;
        Boundary_type mbctop;

        /**
         * Structure containing the boundary options and values per 3d field.
         */
        struct Field3dBc
        {
            double bot; ///< Value of the bottom boundary.
            double top; ///< Value of the top boundary.
            Boundary_type bcbot; ///< Switch for the bottom boundary.
            Boundary_type bctop; ///< Switch for the top boundary.
        };

        typedef std::map<std::string, Field3dBc*> BcMap;
        BcMap sbc;

        // Variables to handle time dependency.
        std::string swtimedep;
        std::vector<double> timedeptime;
        std::vector<std::string> timedeplist;
        std::map<std::string, double*> timedepdata;

        void process_bcs(Input *); ///< Process the boundary condition settings from the ini file.

        void process_time_dependent(Input *); ///< Process the time dependent settings from the ini file.

        void set_bc(double*, double*, double*, Boundary_type, double, double, double); ///< Set the values for the boundary fields.

        // GPU functions and variables
        void set_bc_g(double*, double*, double*, Boundary_type, double, double, double); ///< Set the values for the boundary fields.

    private:
        virtual void update_bcs();       ///< Update the boundary values.
        virtual void update_slave_bcs(); ///< Update the slave boundary values.

        void calc_ghost_cells_bot_2nd(double*, double*, Boundary_type, double*, double*); ///< Calculate the bottom ghost cells with 2nd-order accuracy.
        void calc_ghost_cells_top_2nd(double*, double*, Boundary_type, double*, double*); ///< Calculate the top ghost cells with 2nd-order accuracy.
        void calc_ghost_cells_bot_4th(double*, double*, Boundary_type, double*, double*); ///< Calculate the bottom ghost cells with 4th-order accuracy.
        void calc_ghost_cells_top_4th(double*, double*, Boundary_type, double*, double*); ///< Calculate the top ghost cells with 4th-order accuracy.

        void calc_ghost_cells_botw_4th(double*); ///< Calculate the bottom ghost cells for the vertical velocity with 4th order accuracy.
        void calc_ghost_cells_topw_4th(double*); ///< Calculate the top ghost cells for the vertical velocity with 4th order accuracy.

        void calc_ghost_cells_botw_diff_4th(double*); ///< Calculate the bottom ghost cells for the vertical velocity with 4th order accuracy.
        void calc_ghost_cells_topw_diff_4th(double*); ///< Calculate the top ghost cells for the vertical velocity with 4th order accuracy.
};
#endif
