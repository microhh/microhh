//
// Created by Mirjam Tijhuis on 27/03/2023.
//

#ifndef BACKGROUND_PROFS_H
#define BACKGROUND_PROFS_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <memory>
#include "Types.h"
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
    void exec(Thermo<TF>&);
    void exec_stats(Stats<TF>&);
    void update_time_dependent(Timeloop<TF>&);

    void get_tp(Array<Float,2>&, Array<Float,2>&, Array<Float,2>&, Array<Float,2>&);
    void get_gasses(Gas_concs gas_concs_col);
    void get_aerosols(Aerosol_concs aerosol_concs_col);

private:
    Master& master;
    Grid<TF>& grid;
    Fields<TF>& fields;

    // Case switches
    bool sw_update_background;
    bool sw_aerosol;
    double dt_rad;
    unsigned long idt_rad;

//    const TF n_era_layers = 136;
//    const TF n_era_levels = 137;
    TF n_era_layers;
    TF n_era_levels;

    // Arrays
    // to fill with input values
    // temperature and pressure
    std::vector<TF> t_lay;
    std::vector<TF> t_lev;
    std::vector<TF> p_lay;
    std::vector<TF> p_lev;
    //gasses
    std::vector<TF> h2o;
    std::vector<TF> o3;
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
