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

#include <stdio.h>
#include <exception>
#include <string>

#ifndef TOOLS_H
#define TOOLS_H

/* CUDA error checking, from: http://choorucode.com/2011/03/02/how-to-do-error-checking-in-cuda/
   In debug mode, CUDACHECKS is defined and all kernel calls are checked with cudaCheckError().
   All CUDA api calls are always checked with cudaSafeCall() */
#define cuda_safe_call(err) Tools_g::__cuda_safe_call(err, __FILE__, __LINE__)
#define cuda_check_error()  Tools_g::__cuda_check_error(__FILE__, __LINE__)

namespace Tools_g
{
    enum Reduce_type {Sum_type, Max_type}; ///< Enumerator holding the different reduction types
    const int reduce_max_threads = 512;    ///< Maximum number of threads used in reduce algorithms

    template<typename TF>
    void reduce_interior(const TF *, TF *, int, int, int, int, int, int, int, int, int, int, Reduce_type);
    template<typename TF>
    void reduce_all(const TF *, TF *, int, int, int, Reduce_type, TF);
    template<typename TF> __global__
    void set_to_val(TF* __restrict__, int, TF);
    template<typename TF> __global__
    void mult_by_val(TF* __restrict__, int, TF);

    struct cuda_exception : public std::exception
    {
        cuda_exception(cudaError err, const char *file, const int line);
        cuda_exception(cudaError err, std::string msg);
        const char *what() const throw();
        cudaError error() const;

        private:
            const cudaError err_;
            const std::string message_;
    };

    // Wrapper to check for errors in CUDA api calls (e.g. cudaMalloc)
    inline void __cuda_safe_call(cudaError err, const char *file, const int line)
    {
        if (cudaSuccess != err)
        {
            printf("cudaSafeCall() failed at %s:%i : %s\n", file, line, cudaGetErrorString(err));
            throw cuda_exception(err, file, line);
        }
    }

    // Function to check for errors in CUDA kernels. Call directly after kernel.
    inline void __cuda_check_error(const char *file, const int line)
    {
        #ifdef CUDACHECKS
        cudaError err = cudaGetLastError();
        if (cudaSuccess != err)
        {
            printf("cudaCheckError() failed at %s:%i : %s\n", file, line, cudaGetErrorString( err ) );
            throw cuda_exception(err, file, line);
        }

        err = cudaDeviceSynchronize();
        if(cudaSuccess != err)
        {
            printf("cudaCheckError() with sync failed at %s:%i : %s\n", file, line, cudaGetErrorString( err ) );
            throw cuda_exception(err, file, line);
        }
        #endif
    }
}

#endif
