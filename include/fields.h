/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#ifndef FIELDS
#define FIELDS
#include <memory>
#include <map>
#include <vector>
#include "field3d.h"

class Master;
class Input;
template<typename> class Grid;

//class Stats;
//struct Mask;

// typedef std::map<std::string, Field3d*> Field_map;
template<typename TF>
using Field_map = std::map<std::string, std::shared_ptr<Field3d<TF>>>;

template<typename TF>
class Fields
{
    public:
        Fields(Master&, Grid<TF>&, Input&); ///< Constructor of the fields class.
        ~Fields(); ///< Destructor of the fields class.

        void init();         ///< Initialization of the field arrays.
        void create(Input&, Data_block&); ///< Initialization of the fields (random perturbations, vortices).
        // void create_stats(); ///< Initialization of the fields statistics.

        // void exec();
        // void get_mask(Field3d*, Field3d*, Mask*);
        // void exec_stats(Mask*);

        void init_momentum_field  (std::string, std::string, std::string);
        void init_prognostic_field(std::string, std::string, std::string);
        void init_diagnostic_field(std::string, std::string, std::string);
        void init_tmp_field       (std::string, std::string, std::string);

        void save(int);
        void load(int);

        double check_momentum();
        double check_tke();
        double check_mass();

        // void set_calc_mean_profs(bool);
        // void set_minimum_tmp_fields(int);

        // void exec_cross();
        // void exec_dump();

        Field_map<TF> a;  ///< Map containing all field3d instances
        Field_map<TF> ap; ///< Map containing all prognostic field3d instances
        Field_map<TF> at; ///< Map containing all tendency field3d instances

        Field_map<TF> mp; ///< Map containing all momentum field3d instances
        Field_map<TF> mt; ///< Map containing all momentum tendency field3d instances

        Field_map<TF> sd; ///< Map containing all diagnostic scalar field3d instances
        Field_map<TF> sp; ///< Map containing all prognostic scalar field3d instances
        Field_map<TF> st; ///< Map containing all prognostic scalar tendency field3d instances

        Field_map<TF> atmp; ///< Map containing all temporary field3d instances

        std::vector<TF> rhoref;  ///< Reference density at full levels 
        std::vector<TF> rhorefh; ///< Reference density at half levels

        // TODO remove these to and bring them to diffusion model
        TF visc;

        /* 
         *Device (GPU) functions and variables
         */
        /*
        enum Offset_type {Offset, No_offset};

        void prepare_device();  ///< Allocation of all fields at device 
        void forward_device();  ///< Copy of all fields from host to device
        void backward_device(); ///< Copy of all fields required for statistics and output from device to host
        void clear_device();    ///< Deallocation of all fields at device

        void forward_field_device_3d (double*, double*, Offset_type); ///< Copy of a single 3d field from host to device
        void forward_field_device_2d (double*, double*, Offset_type); ///< Copy of a single 2d field from host to device
        void forward_field_device_1d (double*, double*, int);         ///< Copy of a single array from host to device
        void backward_field_device_3d(double*, double*, Offset_type); ///< Copy of a single 3d field from device to host
        void backward_field_device_2d(double*, double*, Offset_type); ///< Copy of a single 2d field from device to host
        void backward_field_device_1d(double*, double*, int);         ///< Copy of a single array from device to host

        double* rhoref_g;  ///< Reference density at full levels at device
        double* rhorefh_g; ///< Reference density at half levels at device
        */

    private:
        Master& master;
        Grid<TF>& grid;
        // Stats*  stats;

        bool calc_mean_profs;

        int n_tmp_fields;   ///< Number of temporary fields.

        // cross sections
        // std::vector<std::string> crosslist; ///< List with all crosses from the ini file.
        // std::vector<std::string> dumplist;  ///< List with all 3d dumps from the ini file.

        // Cross sections split per type.
        // std::vector<std::string> crosssimple;
        // std::vector<std::string> crosslngrad;   
        // std::vector<std::string> crossbot;
        // std::vector<std::string> crosstop;
        // std::vector<std::string> crossfluxbot;
        // std::vector<std::string> crossfluxtop;

        // void check_added_cross(std::string, std::string, std::vector<std::string>*, std::vector<std::string>*);

        // // masks
        // void calc_mask_wplus(double*, double*, double*, int*, int*, int*, double*);
        // void calc_mask_wmin (double*, double*, double*, int*, int*, int*, double*);

        // perturbations
        TF rndamp;
        TF rndz;
        TF rndexp;
        TF vortexamp;
        int vortexnpair;
        std::string vortexaxis;

        // Kernels for the check functions.
        // double calc_momentum_2nd(double*, double*, double*, double*);
        // double calc_tke_2nd     (double*, double*, double*, double*);
        // double calc_mass        (double*, double*);

        void add_mean_profs(Data_block&);
        // int add_mean_prof(Input*, std::string, double*, double);
        // int randomize    (Input*, std::string, double*);
        // int add_vortex_pair(Input*);

        // statistics
        // double* umodel;
        // double* vmodel;

        /* 
         *Device (GPU) functions and variables
         */
        // void forward_field3d_device(Field3d *);  ///< Copy of a complete Field3d instance from host to device
        // void backward_field3d_device(Field3d *); ///< Copy of a complete Field3d instance from device to host
};
#endif
