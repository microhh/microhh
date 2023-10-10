/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef BACKGROUND_PROFS_H
#define BACKGROUND_PROFS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <memory>
#include "types.h"
#include "Source_functions.h"
#include "Gas_concs.h"

using Aerosol_concs = Gas_concs;

class Master;
class Input;
class Netcdf_handle;
class Netcdf_file;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timedep;
template<typename> class Timeloop;
template<typename> class Stats;
template<typename> class Thermo;
template<typename> class Field3d;

template<typename TF>
class Background
{
public:
    Background(Master&, Grid<TF>&, Fields<TF>&, Input&);
    ~Background();

    void init(Netcdf_handle&, Timeloop<TF>&);
    void create(Input&, Netcdf_handle&, Stats<TF>&);
    void exec_stats(Stats<TF>&);
    void update_time_dependent(Timeloop<TF>&);

    void get_tpm(Array<Float,2>&, Array<Float,2>&, Array<Float,2>&, Array<Float,2>&, Gas_concs&);
    void get_gasses(Gas_concs&);
    void get_aerosols(Aerosol_concs&);

    int get_n_era_levels() const { return n_era_levels; }

private:
    Master& master;
    Grid<TF>& grid;
    Fields<TF>& fields;

    // Case switches
    bool sw_update_background;
    bool sw_aerosol;
    bool sw_aerosol_timedep;
    double dt_rad;
    unsigned long idt_rad;

    TF n_era_layers;
    TF n_era_levels;

    std::map<std::string, Timedep<TF>*> tdep_gases;
    std::vector<std::string> gaslist;        ///< List of gases that have timedependent background profiles.
    std::map<std::string, std::vector<TF>> gasprofs; ///< Map of profiles with gases stored by its name.

    // Arrays
    // to fill with input values
    // temperature, pressure and moisture
    std::vector<TF> t_lay;
    std::vector<TF> t_lev;
    std::vector<TF> p_lay;
    std::vector<TF> p_lev;
    std::vector<TF> h2o;
    //aerosols
    std::vector<TF> aermr01;
    std::vector<TF> aermr02;
    std::vector<TF> aermr03;
    std::vector<TF> aermr04;
    std::vector<TF> aermr05;
    std::vector<TF> aermr06;
    std::vector<TF> aermr07;
    std::vector<TF> aermr08;
    std::vector<TF> aermr09;
    std::vector<TF> aermr10;
    std::vector<TF> aermr11;

    std::unique_ptr<Timedep<TF>> tdep_t_lay;
    std::unique_ptr<Timedep<TF>> tdep_t_lev;
    std::unique_ptr<Timedep<TF>> tdep_p_lay;
    std::unique_ptr<Timedep<TF>> tdep_p_lev;
    std::unique_ptr<Timedep<TF>> tdep_h2o;
    std::unique_ptr<Timedep<TF>> tdep_o3;
    std::unique_ptr<Timedep<TF>> tdep_aermr01;
    std::unique_ptr<Timedep<TF>> tdep_aermr02;
    std::unique_ptr<Timedep<TF>> tdep_aermr03;
    std::unique_ptr<Timedep<TF>> tdep_aermr04;
    std::unique_ptr<Timedep<TF>> tdep_aermr05;
    std::unique_ptr<Timedep<TF>> tdep_aermr06;
    std::unique_ptr<Timedep<TF>> tdep_aermr07;
    std::unique_ptr<Timedep<TF>> tdep_aermr08;
    std::unique_ptr<Timedep<TF>> tdep_aermr09;
    std::unique_ptr<Timedep<TF>> tdep_aermr10;
    std::unique_ptr<Timedep<TF>> tdep_aermr11;
};

#endif //BACKGROUND_PROFS_H
