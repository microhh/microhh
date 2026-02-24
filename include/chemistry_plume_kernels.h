/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#ifndef CHEMISTRY_PLUME_KERNELS_H
#define CHEMISTRY_PLUME_KERNELS_H

#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

template<typename TF>
CUDA_MACRO inline TF epr(
        const TF a1, const TF c1,
        const TF a2, const TF c2,
        const TF a3, const TF c3,
        const TF mmult, const TF ch2o, const TF temp)
{
    const TF k1 = a1 * exp(c1/temp);
    const TF k2 = a2 * exp(c2/temp) * mmult;
    const TF k3 = a3 * exp(c3/temp) * ch2o;
    return (k1 + k2) * (TF(1.0) + k3);
}

template<typename TF>
CUDA_MACRO inline TF arr3(
        const TF a0, const TF b0, const TF temp)
{
    return a0 * exp(b0 / temp);
}

template<typename TF>
CUDA_MACRO inline TF troe_no2oh(
        const TF kzero, const TF mzero, const TF kinf,
        const TF fmulti, const TF mn2, const TF temp)
{
    const TF k0t   = (kzero * pow((TF(300)/temp), mzero)) * mn2;
    const TF kinft = kinf;
    const TF znn   = TF(0.75) - (TF(1.27) * log10(TF(0.41)));
    return (k0t * kinft) / (k0t + kinft) * pow(fmulti, (log10(TF(0.41)) / (TF(1) + pow((log10(k0t/kinft))/znn, 2))));
}

template<typename TF>
CUDA_MACRO inline TF troe_cooh(
        const TF kzero, const TF mzero, const double kinf,
        const TF minf, const TF k3, const TF c3,
        const TF k4, const TF c4, const TF fmulti,
        const TF mn2, const TF temp)
{
    const TF k0t   = (kzero * pow((TF(300)/temp), mzero)) * mn2;
    const TF kinft = kinf * pow((TF(300)/temp), minf);
    const TF kval3 = k3 * pow((TF(300)/temp), c3);
    const TF kval4 = (k4 * pow((TF(300)/temp), c4)) / mn2;
    const TF kcooh = k0t / (TF(1) + k0t / kinft) * pow(fmulti, TF(1) / (TF(1) + pow(log10(k0t/kinft), 2)));
    return kcooh + (kval3 / (TF(1) + kval3 / kval4) * pow(fmulti, TF(1) / (TF(1) + pow(log10(kval3/kval4), 2))));
}

template<typename TF>
CUDA_MACRO inline TF troe_ifs(
        const TF kzero, const TF mzero, const TF kinf,
        const TF minf, const TF fmulti, const TF mn2, const TF temp)
{
    const TF k0t   = (kzero * pow((temp/TF(300)), mzero)) * mn2;
    const TF kinft = kinf * pow((temp/TF(300)), minf);
    const TF znn_b = TF(0.75) - (TF(1.27) * log10(TF(0.35)));
    return (k0t * kinft) / (k0t + kinft) * pow(fmulti, log10(TF(0.35)) / (TF(1) + pow(log10(k0t/kinft)/znn_b, 2)));
}

template<typename TF>
CUDA_MACRO inline TF troe_ifs2(
        const TF kzero, const TF mzero, const TF kinf,
        const TF minf, const TF fmulti, const TF mn2,
        const TF c1, const TF c2, const TF temp)
{
    const TF k0t   = (kzero * pow((temp/TF(300)), mzero)) * mn2 * exp(c1 / temp);
    const TF kinft = kinf * pow((temp/TF(300)), minf) * exp(c2 / temp);
    const TF znn_b = TF(0.75) - (TF(1.27) * log10(TF(0.35)));
    return (k0t * kinft) / (k0t + kinft) * pow(fmulti, log10(TF(0.35)) / (TF(1) + pow(log10(k0t/kinft)/znn_b, 2)));
}

template<typename TF>
CUDA_MACRO inline TF k3rd_iupac(
        const TF kzero, const TF mzero, const TF kinf,
        const TF minf, const TF fmulti, const TF mn2,
        const TF nm, const TF temp)
{
    const TF k0t   = (kzero * pow((temp/TF(300)), mzero)) * mn2;
    const TF kinft = kinf * pow((temp/TF(300)), minf);
    return (k0t / (TF(1) + (k0t/kinft))) * pow(fmulti, TF(1) / (TF(1) + pow((log10(k0t/kinft)/nm), 2)));
}

template<typename TF>
CUDA_MACRO inline TF usr_o3_hv_h2o(
        const TF temp, const TF c_m, const TF c_h2o, const TF j_o1d)
{
    const TF kh2o = TF(1.63e-10) * c_h2o * exp(TF(60)/temp);
    const TF kn2  = TF(2.15e-11) * exp(TF(110.0)/temp) * TF(0.79) * c_m;
    const TF ko2  = TF(3.30e-11) * exp(TF(55.0) /temp) * TF(0.21) * c_m;
    return (kh2o * j_o1d) / (kh2o + kn2 + ko2);
}

template<typename TF>
CUDA_MACRO inline TF rk28(
        const TF k0a, const TF k0ea,
        const TF k2a, const TF k2ea,
        const TF k3a, const TF k3ea,
        const TF mn2, const TF temp)
{
    const TF k0 = k0a * exp(k0ea/temp);
    const TF k2 = (k2a * exp(k2ea/temp)) * mn2;
    const TF k3 = k3a * exp(k3ea/temp);
    return k0 + k3 * k2 / (k3 + k2);
}

#endif // CHEMISTRY_PLUME_KERNELS_H
