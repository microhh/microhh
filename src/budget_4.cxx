/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include <cstdio>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fast_math.h"
#include "finite_difference.h"
#include "model.h"
#include "thermo.h"
#include "diff.h"
#include "advec.h"
#include "force.h"
#include "stats.h"
#include "field3d_operators.h"

#include "budget.h"
#include "budget_4.h"

namespace
{
    template<typename TF>
    void calc_ke(TF* restrict ke, TF* restrict tke,
                 const TF* restrict u, const TF* restrict v, const TF* restrict w,
                 const TF* restrict umodel, const TF* restrict vmodel,
                 const TF utrans, const TF vtrans,
                 const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                 const int icells, const int ijcells)
    {
        using Fast_math::pow2;

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        const TF ci0 = Finite_difference::O4::ci0<TF>;
        const TF ci1 = Finite_difference::O4::ci1<TF>;
        const TF ci2 = Finite_difference::O4::ci2<TF>;
        const TF ci3 = Finite_difference::O4::ci3<TF>;
    
        for (int k=kstart; k<kend; ++k)
        {
            // ke [k] = TF(0.);
            // tke[k] = TF(0.);
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    const TF u2 = ci0*pow2(u[ijk-ii1] + utrans) + ci1*pow2(u[ijk    ] + utrans)
                                + ci2*pow2(u[ijk+ii1] + utrans) + ci3*pow2(u[ijk+ii2] + utrans);
                    const TF v2 = ci0*pow2(v[ijk-jj1] + vtrans) + ci1*pow2(v[ijk    ] + vtrans)
                                + ci2*pow2(v[ijk+jj1] + vtrans) + ci3*pow2(v[ijk+jj2] + vtrans);
                    const TF w2 = ci0*pow2(w[ijk-kk1]) + ci1*pow2(w[ijk]) + ci2*pow2(w[ijk+kk1]) + ci3*pow2(w[ijk+kk2]);
                    ke[ijk] = TF(0.5)*(u2 + v2 + w2);
                }
    
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    const TF u2 = ci0*pow2(u[ijk-ii1] - umodel[k]) + ci1*pow2(u[ijk    ] - umodel[k])
                                + ci2*pow2(u[ijk+ii1] - umodel[k]) + ci3*pow2(u[ijk+ii2] - umodel[k]);
                    const TF v2 = ci0*pow2(v[ijk-jj1] - vmodel[k]) + ci1*pow2(v[ijk    ] - vmodel[k])
                                + ci2*pow2(v[ijk+jj1] - vmodel[k]) + ci3*pow2(v[ijk+jj2] - vmodel[k]);
                    const TF w2 = ci0*pow2(w[ijk-kk1]) + ci1*pow2(w[ijk]) + ci2*pow2(w[ijk+kk1]) + ci3*pow2(w[ijk+kk2]);
                    tke[ijk] = TF(0.5)*(u2 + v2 + w2);
                }
        }
    }
}

template<typename TF>
Budget_4<TF>::Budget_4(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Input& inputin) :
    Budget<TF>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin),
    field3d_operators(masterin, gridin, fieldsin)
{}

template<typename TF>
Budget_4<TF>::~Budget_4()
{}

template<typename TF>
void Budget_4<TF>::init()
{
    auto& gd = grid.get_grid_data();

    umodel.resize(gd.kcells);
    vmodel.resize(gd.kcells);
}

template<typename TF>
void Budget_4<TF>::create(Stats<TF>& stats)
{
    // Add the profiles for the kinetic energy to the statistics.
    stats.add_prof("ke" , "Kinetic energy" , "m2 s-2", "z");
    stats.add_prof("tke", "Turbulent kinetic energy" , "m2 s-2", "z");

    // Add the profiles for the kinetic energy budget to the statistics.
    stats.add_prof("u2_shear" , "Shear production term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_shear" , "Shear production term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("tke_shear", "Shear production term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_shear" , "Shear production term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_turb" , "Turbulent transport term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_visc" , "Viscous transport term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_diss" , "Dissipation term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_pres" , "Pressure transport term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" );
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" );
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh");
    stats.add_prof("uw_rdstr", "Pressure redistribution term in UW budget", "m2 s-3", "zh");

    if (thermo.get_switch() != "0")
    {
        stats.add_prof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh");
        stats.add_prof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" );
        stats.add_prof("uw_buoy" , "Buoyancy production/destruction term in UW budget" , "m2 s-3", "zh");

        stats.add_prof("b2_shear", "Shear production term in B2 budget"   , "m2 s-5", "z");
        stats.add_prof("b2_turb" , "Turbulent transport term in B2 budget", "m2 s-5", "z");
        stats.add_prof("b2_visc" , "Viscous transport term in B2 budget"  , "m2 s-5", "z");
        stats.add_prof("b2_diss" , "Dissipation term in B2 budget"        , "m2 s-5", "z");

        stats.add_prof("bw_shear", "Shear production term in BW budget"   , "m2 s-4", "zh");
        stats.add_prof("bw_turb" , "Turbulent transport term in BW budget", "m2 s-4", "zh");
        stats.add_prof("bw_visc" , "Viscous transport term in BW budget"  , "m2 s-4", "zh");
        stats.add_prof("bw_rdstr", "Redistribution term in BW budget"     , "m2 s-4", "zh");
        stats.add_prof("bw_buoy" , "Buoyancy term in BW budget"           , "m2 s-4", "zh");
        stats.add_prof("bw_diss" , "Dissipation term in BW budget"        , "m2 s-4", "zh");
        stats.add_prof("bw_pres" , "Pressure transport term in BW budget" , "m2 s-4", "zh");
    }

    if (thermo.get_switch() != "0")
    {
        // Add the profiles for the potential energy budget to the statistics.
        stats.add_prof("bsort", "Sorted buoyancy", "m s-2", "z");
        stats.add_prof("zsort", "Height diff buoyancy and sorted buoyancy", "m", "z");
        stats.add_prof("pe"   , "Total potential energy", "m2 s-2", "z");
        stats.add_prof("ape"  , "Available potential energy", "m2 s-2", "z");
        stats.add_prof("bpe"  , "Background potential energy", "m2 s-2", "z");

        // Add the budget terms for the potential energy.
        stats.add_prof("pe_turb", "Turbulent transport term in potential energy budget", "m2 s-3", "z");
        stats.add_prof("pe_visc", "Viscous transport term in potential energy budget", "m2 s-3", "z");
        stats.add_prof("pe_bous", "Boussinesq term in potential energy budget", "m2 s-3", "z");

        // add the budget terms for the background potential energy
        // stats.add_prof("bpe_turb", "Turbulent transport term in background potential energy budget", "m2 s-3", "z");
        // stats.add_prof("bpe_visc", "Viscous transport term in background potential energy budget", "m2 s-3", "z");
        // stats.add_prof("bpe_diss", "Dissipation term in background potential energy budget", "m2 s-3", "z");
    }
}

template<typename TF>
void Budget_4<TF>::exec_stats(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Calculate the mean of the fields.
    field3d_operators.calc_mean_profile(umodel.data(), fields.mp.at("u")->fld.data());
    field3d_operators.calc_mean_profile(vmodel.data(), fields.mp.at("v")->fld.data());

    // Calculate the TKE budget.
    auto ke  = fields.get_tmp();
    auto tke = fields.get_tmp();

    const TF no_offset = 0.;
    const int no_threshold = 0;

    calc_ke(ke->fld.data(), tke->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            umodel.data(), vmodel.data(),
            grid.utrans, grid.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_stats("ke" , *ke , no_offset, no_threshold, {"mean"});
    stats.calc_stats("tke", *tke, no_offset, no_threshold, {"mean"});

    /*
    calc_tke_budget_shear_turb(fields.u->data, fields.v->data, fields.w->data,
                               fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data,
                               umodel, vmodel,
                               m->profs["u2_shear"].data, m->profs["v2_shear"].data, m->profs["tke_shear"].data, m->profs["uw_shear"].data,
                               m->profs["u2_turb"].data, m->profs["v2_turb"].data, m->profs["w2_turb"].data, m->profs["tke_turb"].data, m->profs["uw_turb"].data,
                               grid.dzi4, grid.dzhi4);

    calc_tke_budget(fields.u->data, fields.v->data, fields.w->data, fields.sd["p"]->data,
                    fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data,
                    umodel, vmodel,
                    m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data, m->profs["uw_visc"].data,
                    m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data, m->profs["uw_diss"].data,
                    m->profs["w2_pres"].data, m->profs["tke_pres"].data, m->profs["uw_pres"].data,
                    m->profs["u2_rdstr"].data, m->profs["v2_rdstr"].data, m->profs["w2_rdstr"].data, m->profs["uw_rdstr"].data,
                    grid.dzi4, grid.dzhi4, fields.visc);

    // calculate the buoyancy term of the TKE budget
    if (thermo.get_switch() != "0")
    {
        // store the buoyancy in the tmp1 field
        thermo.get_thermo_field(fields.atmp["tmp1"], fields.atmp["tmp2"], "b", true);

        grid.calc_mean(fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->data, grid.kcells);
        grid.calc_mean(fields.sd["p"]->datamean, fields.sd["p"]->data, grid.kcells);

        calc_tke_budget_buoy(fields.u->data, fields.w->data, fields.atmp["tmp1"]->data,
                             umodel, fields.atmp["tmp1"]->datamean,
                             m->profs["w2_buoy"].data, m->profs["tke_buoy"].data, m->profs["uw_buoy"].data);

        calc_b2_budget(fields.w->data, fields.atmp["tmp1"]->data,
                       fields.atmp["tmp1"]->datamean,
                       m->profs["b2_shear"].data, m->profs["b2_turb"].data, m->profs["b2_visc"].data, m->profs["b2_diss"].data,
                       grid.dzi4, grid.dzhi4,
                       fields.visc);

        calc_bw_budget(fields.w->data, fields.sd["p"]->data, fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data,
                       fields.sd["p"]->datamean, fields.atmp["tmp1"]->datamean,
                       m->profs["bw_shear"].data, m->profs["bw_turb"].data, m->profs["bw_visc"].data,
                       m->profs["bw_buoy"].data, m->profs["bw_rdstr"].data, m->profs["bw_diss"].data, m->profs["bw_pres"].data,
                       grid.dzi4, grid.dzhi4,
                       fields.visc);
    }

    // calculate the potential energy budget
    if (thermo.get_switch() != "0")
    {
        // calculate the sorted buoyancy profile, tmp1 still contains the buoyancy
        stats.calc_sorted_prof(fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, m->profs["bsort"].data);

        // calculate the potential energy back, tmp1 contains the buoyancy, tmp2 will contain height that the local buoyancy
        // will reach in the sorted profile
        calc_pe(fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp2"]->databot, fields.atmp["tmp2"]->datatop,
                grid.z,
                m->profs["bsort"].data,
                m->profs["pe"].data, m->profs["ape"].data, m->profs["bpe"].data,
                m->profs["zsort"].data);


        // calculate the budget of background potential energy, start with this one, because tmp2 contains the needed height
        // which will be overwritten inside of the routine
        // calcBpeBudget(fields.w->data, fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp2"]->databot, fields.atmp["tmp2"]->datatop,
        //               m->profs["bpe_turb"].data, m->profs["bpe_visc"].data, m->profs["bpe_diss"].data,
        //               // TODO put the correct value for visc here!!!!!
        //               m->profs["bsort"].data,
        //               grid.z, grid.dzi4, grid.dzhi4,
        //               fields.visc);

        // calculate the budget of potential energy
        calc_pe_budget(fields.w->data, fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp2"]->datatop,
                       m->profs["pe_turb"].data, m->profs["pe_visc"].data, m->profs["pe_bous"].data,
                       // TODO put the correct value for visc here!!!!!
                       grid.z, grid.zh, grid.dzi4, grid.dzhi4,
                       fields.visc);
    }
    */
}

template class Budget_4<double>;
template class Budget_4<float>;
