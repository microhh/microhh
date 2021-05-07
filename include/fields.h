/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef FIELDS_H
#define FIELDS_H

#include <memory>
#include <map>
#include <vector>
#include <mutex>
#include "field3d.h"
#include "field3d_io.h"
#include "field3d_operators.h"
#include "boundary_cyclic.h"

class Master;
class Input;
class Netcdf_file;
class Netcdf_handle;

template<typename> class Grid;
template<typename> class Soil_grid;
template<typename> class Stats;
template<typename> class Advec;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Soil_field3d;
template<typename> class Field3d_io;
template<typename> class Field3Field3d_operators;
template<typename> struct Mask;

template<typename TF>
using Field_map = std::map<std::string, std::shared_ptr<Field3d<TF>>>;
template<typename TF>
using Soil_field_map = std::map<std::string, std::shared_ptr<Soil_field3d<TF>>>;
template<typename TF>
using Field_2d_map = std::map<std::string, std::vector<TF>>;

enum class Fields_mask_type {Wplus, Wmin};

template<typename TF>
class Fields
{
    public:
        Fields(Master&, Grid<TF>&, Soil_grid<TF>&, Input&); ///< Constructor of the fields class.
        ~Fields(); ///< Destructor of the fields class.

        void init(Input&, Dump<TF>&, Cross<TF>&, const Sim_mode);  ///< Initialization of the field arrays.
        void create(Input&, Netcdf_file&); ///< Initialization of the fields (random perturbations, vortices).
        void create_stats(Stats<TF>&);    ///< Initialization of the fields statistics.
        void create_column(Column<TF>&);  ///< Initialization of the single column output.
        void create_dump(Dump<TF>&);      ///< Initialization of the single column output.
        void create_cross(Cross<TF>&);    ///< Initialization of the single column output.

        void exec();
        void get_mask(Stats<TF>&, std::string);
        void exec_stats(Stats<TF>&);   ///< Calculate the statistics
        void exec_column(Column<TF>&);   ///< Output the column

        void reset_tendencies();

        void init_momentum_field(
                const std::string&, const std::string&,
                const std::string&, const std::string&,
                const std::array<int,3>&);

        void init_prognostic_field(
                const std::string&, const std::string&,
                const std::string&, const std::string&,
                const std::array<int,3>&);

        void init_diagnostic_field(
                const std::string&, const std::string&,
                const std::string&, const std::string&,
                const std::array<int,3>&);

        void init_prognostic_soil_field(const std::string&, const std::string&, const std::string&);
        void init_prognostic_2d_field(const std::string&);

        std::string simplify_unit(const std::string, const std::string, const int = 1, const int = 1);
        void init_tmp_field();

        #ifdef USECUDA
        void init_tmp_field_g();
        #endif

        void save(int);
        void load(int);

        TF check_momentum();
        TF check_tke();
        TF check_mass();

        bool has_mask(std::string);

        void set_calc_mean_profs(bool);

        void exec_cross(Cross<TF>&, unsigned long);
        void exec_dump(Dump<TF>&, unsigned long);

        Field_map<TF> a;  ///< Map containing all field3d instances.
        Field_map<TF> ap; ///< Map containing all prognostic field3d instances.
        Field_map<TF> at; ///< Map containing all tendency field3d instances.

        Field_map<TF> mp; ///< Map containing all momentum field3d instances.
        Field_map<TF> mt; ///< Map containing all momentum tendency field3d instances.

        Field_map<TF> sd; ///< Map containing all diagnostic scalar field3d instances.
        Field_map<TF> sp; ///< Map containing all prognostic scalar field3d instances.
        Field_map<TF> st; ///< Map containing all prognostic scalar tendency field3d instances.

        Soil_field_map<TF> sps; ///< Map containing all prognostic soil scalar fields.
        Soil_field_map<TF> sts; ///< Map containing all prognostic soil scalar tendencies.

        Field_2d_map<TF> ap2d; ///< Map containing all prognostic 2D fields.
        Field_2d_map<TF> at2d; ///< Map containing all prognostic 2D field tendencies.

        std::shared_ptr<Field3d<TF>> get_tmp();
        void release_tmp(std::shared_ptr<Field3d<TF>>&);

        #ifdef USECUDA
        std::shared_ptr<Field3d<TF>> get_tmp_g();
        void release_tmp_g(std::shared_ptr<Field3d<TF>>&);
        #endif

        std::vector<TF> rhoref;  ///< Reference density at full levels
        std::vector<TF> rhorefh; ///< Reference density at half levels

        // TODO remove these to and bring them to diffusion model
        TF visc;

        /*
         *Device (GPU) functions and variables
         */
        void prepare_device();  ///< Allocation of all fields at device
        void forward_device();  ///< Copy of all fields from host to device
        void backward_device(); ///< Copy of all fields required for statistics and output from device to host
        void clear_device();    ///< Deallocation of all fields at device

        void forward_field_device_3d (TF*, TF*);       ///< Copy of a single 3d field from host to device
        void forward_field_device_2d (TF*, TF*);       ///< Copy of a single 2d field from host to device
        void forward_field_device_1d (TF*, TF*, int);  ///< Copy of a single array from host to device
        void backward_field_device_3d(TF*, TF*);       ///< Copy of a single 3d field from device to host
        void backward_field_device_2d(TF*, TF*);       ///< Copy of a single 2d field from device to host
        void backward_field_device_1d(TF*, TF*, int);  ///< Copy of a single array from device to host

        TF* rhoref_g;  ///< Reference density at full levels at device
        TF* rhorefh_g; ///< Reference density at half levels at device


    private:
        Master& master;
        Grid<TF>& grid;
        Soil_grid<TF>& soil_grid;
        Field3d_io<TF> field3d_io;
        Field3d_operators<TF> field3d_operators;
        Boundary_cyclic<TF> boundary_cyclic;

        bool calc_mean_profs;

        int n_tmp_fields;   ///< Number of temporary fields.

        std::vector<std::shared_ptr<Field3d<TF>>> atmp;
        std::vector<std::shared_ptr<Field3d<TF>>> atmp_g;

        std::mutex tmp_fld_mutex;

        // cross sections
        std::vector<std::string> crosslist; ///< List with all crosses from the ini file.
        std::vector<std::string> dumplist;  ///< List with all 3d dumps from the ini file.

        // Cross sections split per type.
        std::vector<std::string> cross_simple;
        std::vector<std::string> cross_lngrad;
        std::vector<std::string> cross_bot;
        std::vector<std::string> cross_top;
        std::vector<std::string> cross_fluxbot;
        std::vector<std::string> cross_fluxtop;
        std::vector<std::string> cross_path;

        void check_added_cross(
                const std::string&,
                const std::string&,
                std::vector<std::string>&,
                std::vector<std::string>&);

        // // masks
        std::vector<std::string> available_masks;   // Vector with the masks that fields can provide
        // void calc_mask_wplus(double*, double*, double*, int*, int*, int*, double*);
        // void calc_mask_wmin (double*, double*, double*, int*, int*, int*, double*);

        // perturbations
        TF rndamp;
        TF rndz;
        TF rndexp;
        TF vortexamp;
        int vortexnpair;
        std::string vortexaxis;

        void add_mean_profs(Netcdf_handle&);
        // int add_mean_prof(Input*, std::string, double*, double);
        void randomize(Input&, std::string, TF* const restrict);
        void add_vortex_pair(Input&);

        // statistics
        std::vector<TF> umodel;
        std::vector<TF> vmodel;

        // double* umodel;
        // double* vmodel;

        /*
         *Device (GPU) functions and variables
         */
        void forward_field3d_device(Field3d<TF> *);  ///< Copy of a complete Field3d instance from host to device
        void backward_field3d_device(Field3d<TF> *); ///< Copy of a complete Field3d instance from device to host
};
#endif
