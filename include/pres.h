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

#ifndef PRES
#define PRES

class Master;
template<typename> class Grid;
template<typename> class Fields;

#ifdef USECUDA
#include <cufft.h>
#endif

template<typename TF>
class Pres
{
    public:
        Pres(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Pres();

        static std::shared_ptr<Pres> factory(Master&, Grid<TF>&, Fields<TF>&, Input&, const std::string); ///< Factory function for pres class generation.

        virtual void init();
        virtual void set_values();

        virtual void exec(double);
        virtual TF check_divergence();

        virtual void prepare_device();

    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

    private:
        #ifdef USECUDA
        void make_cufft_plan();
        void fft_forward (double*, double*, double*);
        void fft_backward(double*, double*, double*);

        bool FFTPerSlice;
        cufftHandle iplanf;
        cufftHandle jplanf;
        cufftHandle iplanb;
        cufftHandle jplanb;
        #endif

#ifdef USECUDA
        void check_cufft_memory();
#endif
};
#endif
