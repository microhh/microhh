
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

#include <iostream>

#include "chemistry.h"
#include "tools.h"
#include "fields.h"
#include "thermo.h"
#include "constants.h"
#include "chemistry_plume_kernels.h"

namespace
{
    template<typename TF>
    __global__ void pss_g(
            TF* const restrict thno3,
            TF* const restrict th2o2,
            TF* const restrict tco,
            TF* const restrict thcho,
            TF* const restrict trooh,
            TF* const restrict tc3h6,
            TF* const restrict to3,
            TF* const restrict tno,
            TF* const restrict tno2,
            const TF* const restrict hno3,
            const TF* const restrict h2o2,
            const TF* const restrict co,
            const TF* const restrict hcho,
            const TF* const restrict rooh,
            const TF* const restrict c3h6,
            const TF* const restrict o3,
            const TF* const restrict no,
            const TF* const restrict no2,
            const TF* const restrict jval,
            const TF* const restrict vdo3,
            const TF* const restrict vdno,
            const TF* const restrict vdno2,
            const TF* const restrict vdhno3,
            const TF* const restrict vdh2o2,
            const TF* const restrict vdrooh,
            const TF* const restrict vdhcho,
            const TF* const restrict tprof,
            const TF* const restrict qprof,
            const TF* const restrict dzi,
            const TF* const restrict rhoref,
            const TF sdt,
            const int istart,
            const int iend,
            const int jstart,
            const int jend,
            const int kstart,
            const int kend,
            const int jstride,
            const int kstride)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jstride + k*kstride;
            const int ij = i + j*jstride;

            const TF c_m = TF(1e-3) * rhoref[k] * Constants::Na<TF> * Constants::xmair_i<TF>;   // molecules/cm3 for chemistry!

            // From ppb (units mixing ratio) to molecules/cm3 --> changed: now mol/mol unit for transported tracers:
            const TF cfactor = c_m;
            const TF cfac_i = TF(1) / cfactor;
            const TF c_h2 = TF(500e-9) * cfactor ; // 500 ppb --> #/cm3

            // rate constants on horizontal average quantities,
            const TF temp = tprof[k];
            const TF c_h2o = std::max(qprof[k] * Constants::xmair<TF> * c_m * Constants::xmh2o_i<TF>, TF(1));

            const TF rconst0  = arr3(TF(1.7E-12), TF(-940), temp);
            const TF rconst1  = arr3(TF(1.E-14), TF(-490), temp);
            const TF rconst2  = arr3(TF(4.8E-11), TF(250), temp);
            const TF rconst3  = epr(TF(3.E-13), TF(460), TF(2.1E-33), TF(920), TF(1.4E-21), TF(2200), c_m, c_h2o, temp);
            const TF rconst4  = arr3(TF(2.9E-12), TF(-160), temp);
            const TF rconst5  = arr3(TF(2.8E-12), TF(-1800), temp) * c_h2;
            const TF rconst6  = arr3(TF(3.E-12), TF(-1500), temp);
            const TF rconst7  = arr3(TF(1.4E-13), TF(-2470), temp);
            const TF rconst8  = arr3(TF(1.8E-11), TF(110), temp);
            const TF rconst9  = troe_ifs(TF(3.6E-30), TF(4.1), TF(1.9E-12), TF(-0.2), TF(10), c_m, temp);
            const TF rconst10 = troe_ifs2(TF(1.3E-3), TF(-3.5), TF(9.7E14), TF(0.1), TF(10), c_m, TF(-11000), TF(-11080), temp);
            const TF rconst11 = arr3(TF(3.3E-12), TF(270), temp);
            const TF rconst12 = troe_no2oh(TF(3.2E-30), TF(4.5), TF(3.E-11), TF(10), c_m, temp);
            const TF rconst14 = rk28(TF(2.4E-14), TF(460), TF(6.51E-34), TF(1335), TF(2.69E-17), TF(2199), c_m, temp);
            const TF rconst16 = arr3(TF(2.45E-12), TF(-1775), temp);
            const TF rconst17 = arr3(TF(3.8E-13), TF(780), temp) * (TF(1)-(TF(1)/(TF(1) + arr3(TF(498), TF(-1160), temp))));
            const TF rconst18 = arr3(TF(3.8E-13), TF(780), temp) * (TF(1)/(TF(1.) + arr3(TF(498.), TF(-1160), temp)));
            const TF rconst19 = arr3(TF(2.8E-12), TF(300), temp);
            const TF rconst21 = arr3(TF(3.8E-12), TF(200), temp);
            const TF rconst22 = arr3(TF(5.5E-12), TF(125), temp);
            const TF rconst24 = troe_cooh(TF(5.9E-33), TF(1.4), TF(1.1E-12), TF(-1.3), TF(1.5E-13), TF(-0.6), TF(2.1E9), TF(-6.1), TF(0.6), c_m, temp);
            const TF rconst25 = arr3(TF(9.5E-14), TF(390), temp);
            const TF rconst26 = arr3(TF(5.5E-15), TF(-1880), temp);
            const TF rconst27 = k3rd_iupac(TF(8.6E-27), TF(3.5), TF(3.E-11), TF(1), TF(0.6), c_m, TF(0.5), temp);
            const TF rconst28 = arr3(TF(4.6E-13), TF(-1155), temp);
            const TF rconst29 = usr_o3_hv_h2o(temp, c_m, c_h2o, jval[Jval::o31d]);
            const TF rconst37 = TF(0);
            const TF rconst38 = TF(0);

            const TF fix_ch4 = TF(1800e-9) * cfactor;   // methane concentration

            // Convert to molecules per cm3 and add tendencies of other processes.
            const TF var_hno3 = std::max((hno3[ijk] + thno3[ijk] * sdt) * cfactor, TF(0));
            const TF var_h2o2 = std::max((h2o2[ijk] + th2o2[ijk] * sdt) * cfactor, TF(0));
            const TF var_co   = std::max((co[ijk]   + tco[ijk]   * sdt) * cfactor, TF(0));
            const TF var_hcho = std::max((hcho[ijk] + thcho[ijk] * sdt) * cfactor, TF(0));
            const TF var_rooh = std::max((rooh[ijk] + trooh[ijk] * sdt) * cfactor, TF(0));
            const TF var_rh   = std::max((c3h6[ijk] + tc3h6[ijk] * sdt) * cfactor, TF(0));
            const TF var_o3   = std::max((o3[ijk]   + to3[ijk]   * sdt) * cfactor, TF(0));
            const TF var_no   = std::max((no[ijk]   + tno[ijk]   * sdt) * cfactor, TF(0));
            const TF var_no2  = std::max((no2[ijk]  + tno2[ijk]  * sdt) * cfactor, TF(0));

            const TF rconst39 = (k == kstart) ? vdo3[ij]   * dzi[k] : TF(0);
            const TF rconst40 = (k == kstart) ? vdno[ij]   * dzi[k] : TF(0);
            const TF rconst41 = (k == kstart) ? vdno2[ij]  * dzi[k] : TF(0);
            const TF rconst42 = (k == kstart) ? vdhno3[ij] * dzi[k] : TF(0);
            const TF rconst43 = (k == kstart) ? vdh2o2[ij] * dzi[k] : TF(0);
            const TF rconst44 = (k == kstart) ? vdhcho[ij] * dzi[k] : TF(0);
            const TF rconst45 = (k == kstart) ? vdrooh[ij] * dzi[k] : TF(0);

            // QSSA Auxiliary Variables
            TF p_oh, l_oh;
            TF p_ho2, l_ho2;
            TF p_ro2, l_ro2;
            TF p_no3, l_no3;
            TF p_n2o5, l_n2o5;

            // Reset QSSA iteration for each grid point. Otherwise this causes a race condition in vectorised or GPU code.
            TF fix_oh   = TF(0);
            TF fix_ho2  = TF(0);
            TF fix_ro2  = TF(0);
            TF fix_no3  = TF(0);
            TF fix_n2o5 = TF(0);

            // =========================================================================
            // QSSA Iteration: Solve P = L * F for OH, HO2, RO2, NO3, N2O5
            // =========================================================================
            for (int i_qssa=0; i_qssa<2; i_qssa++)
            {
                // ---------------------------------------------------------------------
                // 1. Solve N2O5 (Index fix_n2o5)
                // ---------------------------------------------------------------------
                // Prod: NO2 + NO3 (A9)
                p_n2o5 = rconst9*var_no2*fix_no3;               // RF[9]: NO2(V6) + NO3(F6)

                // Loss: N2O5 -> NO2+NO3 (A10), N2O5 -> 2HNO3 (A15), Loss(A31)
                l_n2o5 = rconst10                         // RF[10]
                       + TF(0.0004)                      // RF[15] (Hardcoded k)
                       + jval[Jval::n2o5];                   // RF[31]

                fix_n2o5 = (l_n2o5 > TF(1e-30)) ? p_n2o5 / l_n2o5 : TF(0);

                // ---------------------------------------------------------------------
                // 2. Solve NO3 (Index fix_no3)
                // ---------------------------------------------------------------------
                // Prod: NO2+O3 (A7), N2O5->NO2+NO3 (A10), HNO3+OH (A14)
                p_no3 = rconst7*var_no2*var_o3                 // RF[7]: NO2(V6) + O3(V7)
                      + rconst10*fix_n2o5                     // RF[10]: N2O5(F7)
                      + rconst14*var_hno3*fix_oh;               // RF[14]: HNO3(V0) + OH(F3)

                // Loss: NO(A8), NO2(A9), HO2(A13), RO2(A20), HCHO(A23), RH(A28), Loss(A32)
                l_no3 = rconst8*var_no                      // RF[8]: NO(V8)
                      + rconst9*var_no2                      // RF[9]: NO2(V6) (To N2O5)
                      + TF(4.0e-12)*fix_ho2                 // RF[13]: HO2(F4) (Hardcoded k)
                      + TF(1.2e-12)*fix_ro2                 // RF[20]: RO2(F5) (Hardcoded k)
                      + TF(5.8e-16)*var_hcho                 // RF[23]: HCHO(V3) (Hardcoded k)
                      + rconst28*var_rh                     // RF[28]: RH(V5)
                      + jval[Jval::no3];                     // RF[32]: Loss

                fix_no3 = (l_no3 > TF(1e-30)) ? p_no3 / l_no3 : TF(0);

                // ---------------------------------------------------------------------
                // 3. Solve RO2 (Index fix_ro2)
                // ---------------------------------------------------------------------
                // Prod: CH4+OH(A16), ROOH+OH(0.6*A21), RH+O3(0.31*A26), RH+OH(A27)
                p_ro2 = rconst16*fix_oh*fix_ch4                // RF[16]: OH(F3) + CH4(F0)
                      + TF(0.6)*rconst21*var_rooh*fix_oh        // RF[21]: ROOH(V4) + OH(F3)
                      + TF(0.31)*rconst26*var_rh*var_o3       // RF[26]: RH(V5) + O3(V7)
                      + rconst27*var_rh*fix_oh;               // RF[27]: RH(V5) + OH(F3)

                // Loss: HO2(A17, A18), NO(A19), NO3(A20), RO2(2*A25)
                l_ro2 = rconst17*fix_ho2                     // RF[17]: HO2(F4)
                      + rconst18*fix_ho2                     // RF[18]: HO2(F4)
                      + rconst19*var_no                     // RF[19]: NO(V8)
                      + TF(1.2e-12)*fix_no3                 // RF[20]: NO3(F6) (Uses F6 now)
                      + TF(2.0)*rconst25*fix_ro2;            // RF[25]: RO2(F5) - Quadratic

                fix_ro2 = (l_ro2 > TF(1e-30)) ? p_ro2 / l_ro2 : TF(0);

                // ---------------------------------------------------------------------
                // 4. Solve HO2 and OH (Index fix_ho2 fix_oh)
                // ---------------------------------------------------------------------
                // HO2 Prod: O3+OH(A0), H2O2+OH(A4), OH+M(A5), RO2+NO(A19), RO2+NO3(A20),
                //       HCHO+OH(A22), CO+OH(A24), RO2+RO2(0.74*A25), RH+O3(0.19*A26),
                //       ROOH+hv(A33), HCHO+hv(2*A35)
                p_ho2 = rconst0*var_o3*fix_oh                 // RF[0]: O3(V7) + OH(F3)
                      + rconst4*var_h2o2*fix_oh                 // RF[4]: H2O2(V1) + OH(F3)
                      + rconst5*fix_oh                      // RF[5]: OH(F3)
                      + rconst19*fix_ro2*var_no                // RF[19]: RO2(F5) + NO(V8)
                      + TF(1.2e-12)*fix_no3*fix_ro2            // RF[20]: NO3(F6) + RO2(F5)
                      + rconst22*var_hcho*fix_oh                // RF[22]: HCHO(V3) + OH(F3)
                      + rconst24*var_co*fix_oh                  // RF[24]: CO(V2) + OH(F3)
                      + TF(0.74)*rconst25*fix_ro2*fix_ro2       // RF[25]: RO2(F5)^2
                      + TF(0.19)*rconst26*var_rh*var_o3         // RF[26]: RH(V5) + O3(V7)
                      + jval[Jval::ch3o2h]*var_rooh               // RF[33]: ROOH(V4)
                      + TF(2.0)*jval[Jval::ch2or]*var_hcho;      // RF[35]: HCHO(V3)
                // Loss: O3(A1), OH(A2), HO2(2*A3), NO(A11), NO3(A13), RO2(A17, A18)
                l_ho2 = rconst1*var_o3                      // RF[1]: O3(V7)
                      + rconst2*fix_oh                      // RF[2]: OH(F3)
                      + TF(2.0)*rconst3*fix_ho2               // RF[3]: HO2(F4)
                      + rconst11*var_no                     // RF[11]: NO(V8)
                      + TF(4.0e-12)*fix_no3                 // RF[13]: NO3(F6) (Uses F6 now)
                      + rconst17*fix_ro2                     // RF[17]: RO2(F5)
                      + rconst18*fix_ro2;                    // RF[18]: RO2(F5)

                fix_ho2 = (l_ho2 > TF(1e-30)) ? p_ho2 / l_ho2 : TF(0);

                // OH Prod: HO2+O3(A1), HO2+NO(A11), RH+O3(0.33*A26), O3+hv(2*A29),
                //       ROOH+hv(A33), H2O2+hv(2*A36)
                p_oh = rconst1*var_o3*fix_ho2                  // RF[1]: O3(V7) + HO2(F4)
                     + rconst11*fix_ho2*var_no                 // RF[11]: HO2(F4) + NO(V8)
                     + TF(0.33)*rconst26*var_rh*var_o3         // RF[26]: RH(V5) + O3(V7)
                     + TF(2.0)*rconst29*var_o3               // RF[29]: O3(V7)
                     + jval[Jval::ch3o2h]*var_rooh                // RF[33]: ROOH(V4)
                     + TF(2.0)*jval[Jval::h2o2]*var_h2o2;     // RF[36]: H2O2(V1)


                // Loss: O3(A0), HO2(A2), H2O2(A4), M(A5), NO2(A12), HNO3(A14),
                //       CH4(A16), ROOH(0.6*A21), HCHO(A22), CO(A24), RH(A27)
                l_oh = rconst0*var_o3                       // RF[0]: O3(V7)
                     + rconst2*fix_ho2                       // RF[2]: HO2(F4)
                     + rconst4*var_h2o2                       // RF[4]: H2O2(V1)
                     + rconst5                            // RF[5]
                     + rconst12*var_no2                      // RF[12]: NO2(V6)
                     + rconst14*var_hno3                      // RF[14]: HNO3(V0)
                     + rconst16*fix_ch4                      // RF[16]: CH4(F0)
                     + TF(0.6)*rconst21*var_rooh               // RF[21]: ROOH(V4)
                     + rconst22*var_hcho                      // RF[22]: HCHO(V3)
                     + rconst24*var_co                      // RF[24]: CO(V2)
                     + rconst27*var_rh;                     // RF[27]: RH(V5)

                fix_oh = (l_oh > TF(1e-30)) ? p_oh / l_oh : TF(0);
            }

            // Computation of equation rates.
            const TF rf0  = rconst0*var_o3*fix_oh;
            const TF rf1  = rconst1*var_o3*fix_ho2;
            const TF rf3  = rconst3*fix_ho2*fix_ho2;
            const TF rf4  = rconst4*var_h2o2*fix_oh;
            const TF rf6  = rconst6*var_o3*var_no;
            const TF rf7  = rconst7*var_no2*var_o3;
            const TF rf8  = rconst8*var_no*fix_no3;
            const TF rf9  = rconst9*var_no2*fix_no3;
            const TF rf10 = rconst10*fix_n2o5;
            const TF rf11 = rconst11*var_no*fix_ho2;
            const TF rf12 = rconst12*var_no2*fix_oh;
            const TF rf13 = TF(4e-12)*fix_ho2*fix_no3;
            const TF rf14 = rconst14*var_hno3*fix_oh;
            const TF rf15 = TF(0.0004)*fix_n2o5;
            const TF rf17 = rconst17*fix_ho2*fix_ro2;
            const TF rf18 = rconst18*fix_ho2*fix_ro2;
            const TF rf19 = rconst19*var_no*fix_ro2;
            const TF rf20 = TF(1.2e-12)*fix_ro2*fix_no3;
            const TF rf21 = rconst21*var_rooh*fix_oh;
            const TF rf22 = rconst22*var_hcho*fix_oh;
            const TF rf23 = TF(5.8e-16)*var_hcho*fix_no3;
            const TF rf24 = rconst24*var_co*fix_oh;
            const TF rf25 = rconst25*fix_ro2*fix_ro2;
            const TF rf26 = rconst26*var_rh*var_o3;
            const TF rf27 = rconst27*var_rh*fix_oh;
            const TF rf28 = rconst28*var_rh*fix_no3;
            const TF rf29 = rconst29*var_o3;
            const TF rf30 = jval[Jval::no2]*var_no2;
            const TF rf31 = jval[Jval::n2o5]*fix_n2o5;
            const TF rf32 = jval[Jval::no3]*fix_no3;
            const TF rf33 = jval[Jval::ch3o2h]*var_rooh;
            const TF rf34 = jval[Jval::ch2om]*var_hcho;
            const TF rf35 = jval[Jval::ch2or]*var_hcho;
            const TF rf36 = jval[Jval::h2o2]*var_h2o2;
            const TF rf37 = rconst37;
            const TF rf38 = rconst38;
            const TF rf39 = rconst39*var_o3;
            const TF rf40 = rconst40*var_no;
            const TF rf41 = rconst41*var_no2;
            const TF rf42 = rconst42*var_hno3;
            const TF rf43 = rconst43*var_h2o2;
            const TF rf44 = rconst44*var_hcho;
            const TF rf45 = rconst45*var_rooh;

            // Aggregate function
            const TF vdot_hno3 = rf12+rf13-rf14+TF(2)*rf15+rf23-rf42;
            const TF vdot_h2o2 = rf3-rf4-rf36-rf43;
            const TF vdot_co   = rf22-rf24+TF(0.56)*rf26+rf34+rf35;
            const TF vdot_hcho = rf18+rf19+rf20+TF(0.4)*rf21-rf22-rf23+TF(1.37)*rf25+TF(1.04)*rf26+TF(1.5)*rf27+rf33-rf34-rf35-rf44;
            const TF vdot_rooh = rf17-rf21-rf33-rf45;
            const TF vdot_rh   = -rf26-rf27-rf28+rf38;
            const TF vdot_no2  = rf6-rf7+TF(2)*rf8-rf9+rf10+rf11-rf12+rf19+rf20-rf30+rf31+rf32-rf41;
            const TF vdot_o3   = -rf0-rf1-rf6-rf7-rf26-rf29+rf30+rf32-rf39;
            const TF vdot_no   = -rf6-rf8-rf11-rf19+rf30+rf37-rf40;

            thno3[ijk] += vdot_hno3 * cfac_i;
            th2o2[ijk] += vdot_h2o2 * cfac_i;
            tco[ijk]   += vdot_co   * cfac_i;
            thcho[ijk] += vdot_hcho * cfac_i;
            trooh[ijk] += vdot_rooh * cfac_i;
            tc3h6[ijk] += vdot_rh   * cfac_i;
            to3[ijk]   += vdot_o3   * cfac_i;
            tno[ijk]   += vdot_no   * cfac_i;
            tno2[ijk]  += vdot_no2  * cfac_i;
        }
    }
}

#ifdef USECUDA
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo, const double sdt, const double dt)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    // Calculate the mean temperature profile.
    auto temperature = fields.get_tmp_g();
    thermo.get_thermo_field_g(*temperature, "T", true);
    field3d_operators.calc_mean_profile_g(temperature->fld_mean_g, temperature->fld_g);

    pss_g<TF><<<gridGPU, blockGPU>>>(
        fields.st.at("hno3")->fld_g,
        fields.st.at("h2o2")->fld_g,
        fields.st.at("co")  ->fld_g,
        fields.st.at("hcho")->fld_g,
        fields.st.at("rooh")->fld_g,
        fields.st.at("c3h6")->fld_g,
        fields.st.at("o3")  ->fld_g,
        fields.st.at("no")  ->fld_g,
        fields.st.at("no2") ->fld_g,
        fields.sp.at("hno3")->fld_g,
        fields.sp.at("h2o2")->fld_g,
        fields.sp.at("co")  ->fld_g,
        fields.sp.at("hcho")->fld_g,
        fields.sp.at("rooh")->fld_g,
        fields.sp.at("c3h6")->fld_g,
        fields.sp.at("o3")  ->fld_g,
        fields.sp.at("no")  ->fld_g,
        fields.sp.at("no2") ->fld_g,
        jval_g,
        vdo3_g,
        vdno_g,
        vdno2_g,
        vdhno3_g,
        vdh2o2_g,
        vdrooh_g,
        vdhcho_g,
        temperature->fld_mean_g,
        fields.sp.at("qt")->fld_mean_g,
        gd.dzi_g,
        fields.rhoref_g,
        sdt,
        gd.istart,
        gd.iend,
        gd.jstart,
        gd.jend,
        gd.kstart,
        gd.kend,
        gd.icells,
        gd.ijcells);

    fields.release_tmp_g(temperature);
}

template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary)
{
    if (!sw_chemistry)
        return;

    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);

    // There is nothing the GPU can speedup here, simply calculate at CPU and memcpy to device.
    jval[Jval::o31d]   = ifac.fac0 * jo31d[ifac.index0]   + ifac.fac1 * jo31d[ifac.index1];
    jval[Jval::h2o2]   = ifac.fac0 * jh2o2[ifac.index0]   + ifac.fac1 * jh2o2[ifac.index1];
    jval[Jval::no2]    = ifac.fac0 * jno2[ifac.index0]    + ifac.fac1 * jno2[ifac.index1];
    jval[Jval::no3]    = ifac.fac0 * jno3[ifac.index0]    + ifac.fac1 * jno3[ifac.index1];
    jval[Jval::n2o5]   = ifac.fac0 * jn2o5[ifac.index0]   + ifac.fac1 * jn2o5[ifac.index1];
    jval[Jval::ch2or]  = ifac.fac0 * jch2or[ifac.index0]  + ifac.fac1 * jch2or[ifac.index1];
    jval[Jval::ch2om]  = ifac.fac0 * jch2om[ifac.index0]  + ifac.fac1 * jch2om[ifac.index1];
    jval[Jval::ch3o2h] = ifac.fac0 * jch3o2h[ifac.index0] + ifac.fac1 * jch3o2h[ifac.index1];

    cuda_safe_call(cudaMemcpy(jval_g, jval.data(), n_jval*sizeof(TF), cudaMemcpyHostToDevice));

    //deposition->update_time_dependent(
    //        timeloop,
    //        boundary,
    //        vdo3.data(),
    //        vdno.data(),
    //        vdno2.data(),
    //        vdhno3.data(),
    //        vdh2o2.data(),
    //        vdrooh.data(),
    //        vdhcho.data());
}

template<typename TF>
void Chemistry<TF>::prepare_device()
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    const int time_size = time.size();

    // Allocate GPU arrays.
    // Only `jval_g` is enough; interpolation is performed at CPU, and memcpy'd to device.
    jval_g.allocate(n_jval);

    vdo3_g.allocate(gd.ijcells);
    vdno_g.allocate(gd.ijcells);
    vdno2_g.allocate(gd.ijcells);
    vdhno3_g.allocate(gd.ijcells);
    vdh2o2_g.allocate(gd.ijcells);
    vdrooh_g.allocate(gd.ijcells);
    vdhcho_g.allocate(gd.ijcells);

    // Copy data from host to device.
    cuda_safe_call(cudaMemcpy(jval_g,    jval.data(),    n_jval*sizeof(TF),     cudaMemcpyHostToDevice));

    const int memsize_ij = gd.ijcells * sizeof(TF);
    cuda_safe_call(cudaMemcpy(vdo3_g,    vdo3.data(),    memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdno_g,    vdno.data(),    memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdno2_g,   vdno2.data(),   memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdhno3_g,  vdhno3.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdh2o2_g,  vdh2o2.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdrooh_g,  vdrooh.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdhcho_g,  vdhcho.data(),  memsize_ij, cudaMemcpyHostToDevice));
}

template<typename TF>
void Chemistry<TF>::clear_device()
{
    if (!sw_chemistry)
        return;
}
#endif

#ifdef FLOAT_SINGLE
template class Chemistry<float>;
#else
template class Chemistry<double>;
#endif