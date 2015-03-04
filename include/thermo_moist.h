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

#ifndef THERMO_MOIST
#define THERMO_MOIST

#include "thermo.h"

class Master;
class Grid;
class Fields;
class Stats;
struct Mask;

class Thermo_moist : public Thermo
{
    public:
        Thermo_moist(Model*, Input*);
        virtual ~Thermo_moist();

        void init();
        void create(Input*);
        void exec();
        void get_mask(Field3d*, Field3d*, Mask*);
        void exec_stats(Mask*);
        void exec_cross();
        void exec_dump();

        // functions to retrieve buoyancy properties, to be called from other classes
        bool check_thermo_field(std::string name);
        void get_thermo_field(Field3d*, Field3d*, std::string name);
        void get_buoyancy_surf(Field3d*);
        void get_buoyancy_fluxbot(Field3d*);
        void get_prog_vars(std::vector<std::string>*); ///< Retrieve a list of prognostic variables.

#ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
#endif

    private:
        void initStat();  ///< Initialize the thermo statistics
        void initCross(); ///< Initialize the thermo cross-sections
        void initDump();  ///< Initialize the thermo field dumps

        int swupdatebasestate;
        std::string thvar; ///< Name of prognostic potential temperature variable

        // cross sections
        std::vector<std::string> crosslist;        ///< List with all crosses from ini file
        std::vector<std::string> allowedcrossvars; ///< List with allowed cross variables
        std::vector<std::string> dumplist;         ///< List with all 3d dumps from the ini file.

        Stats *stats;

        // masks
        void calcMask_ql    (double*, double*, double*, int *, int *, int *, double *);
        void calcMask_qlcore(double*, double*, double*, int *, int *, int *, double*, double*, double *);

        void calcBuoyancyTend_2nd(double*, double*, double*, double*, double*, double*, double*, double *);
        void calcBuoyancyTend_4th(double*, double*, double*, double*, double*, double*, double*, double *);

        void calcBuoyancy(double*, double*, double*, double*, double*, double *);
        void calcN2(double*, double*, double*, double *); ///< Calculation of the Brunt-Vaissala frequency.
        void calcBaseState(double*, double*, double*, double*, double*, double*, double*, double*, double*, double *);

        void calcLiquidWater(double*, double*, double*, double *);
        void calcBuoyancyBot(double*, double*,
                double*, double*,
                double*, double*,
                double*, double *);
        void calcBuoyancyFluxBot(double*, double*, double*, double*, double*, double *);

        std::string swbasestate;
        double pbot;
        double thvref0; ///< Reference virtual potential temperature in case of Boussinesq

        // REFERENCE PROFILES
        double *thl0;    // Initial thl profile 
        double *qt0;     // Initial qt profile
        double *thvref; 
        double *thvrefh;
        double *exnref;
        double *exnrefh;
        double *pref;
        double *prefh;

        // GPU functions and variables
        double *thvref_g; 
        double *thvrefh_g;
        double *exnref_g;
        double *exnrefh_g;
        double *pref_g;
        double *prefh_g;
};
#endif
