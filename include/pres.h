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

#ifndef PRES_H
#define PRES_H

#ifdef USECUDA
#include <cufft.h>
#endif

#include "field3d_operators.h"
#include "fft.h"

class Master;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class FFT;
template<typename> class Stats;

template<typename TF>
class Pres
{
    public:
        Pres(Master&, Grid<TF>&, Fields<TF>&, FFT<TF>&, Input&);
        virtual ~Pres();

        static std::shared_ptr<Pres> factory(Master&, Grid<TF>&, Fields<TF>&, FFT<TF>&, Input&);

        virtual void init() = 0;
        virtual void set_values() = 0;
        virtual void create(Stats<TF>&) = 0;

        virtual void exec(double, Stats<TF>&) = 0;
        virtual TF check_divergence() = 0;

        #ifdef USECUDA
        virtual void prepare_device() = 0;
        virtual void clear_device() = 0;
        #endif

    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        FFT<TF>& fft;

        Field3d_operators<TF> field3d_operators;

        #ifdef USECUDA
        void make_cufft_plan();
        void fft_forward (TF*, TF*, TF*);
        void fft_backward(TF*, TF*, TF*);

        bool FFT_per_slice;
        bool force_FFT_per_slice;
        cufftHandle iplanf;
        cufftHandle jplanf;
        cufftHandle iplanb;
        cufftHandle jplanb;
        #endif

    private:
        #ifdef USECUDA
        void check_cufft_memory();
        #endif
};
#endif
