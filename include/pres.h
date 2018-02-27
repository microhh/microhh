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

        virtual void init() = 0;
        virtual void set_values() = 0;

        virtual void exec(double) = 0;
        virtual TF check_divergence() = 0;

        #ifdef USECUDA
        virtual void prepare_device() = 0;
        virtual void clear_device() = 0;
        #endif

    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        #ifdef USECUDA
        void make_cufft_plan();
        void fft_forward (TF*, TF*, TF*);
        void fft_backward(TF*, TF*, TF*);

        bool FFT_per_slice;
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
