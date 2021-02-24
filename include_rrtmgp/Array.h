/*
 * This file is part of a C++ interface to the Radiative Transfer for Energetics (RTE)
 * and Rapid Radiative Transfer Model for GCM applications Parallel (RRTMGP).
 *
 * The original code is found at https://github.com/earth-system-radiation/rte-rrtmgp.
 *
 * Contacts: Robert Pincus and Eli Mlawer
 * email: rrtmgp@aer.com
 *
 * Copyright 2015-2020,  Atmospheric and Environmental Research and
 * Regents of the University of Colorado.  All right reserved.
 *
 * This C++ interface can be downloaded from https://github.com/earth-system-radiation/rte-rrtmgp-cpp
 *
 * Contact: Chiel van Heerwaarden
 * email: chiel.vanheerwaarden@wur.nl
 *
 * Copyright 2020, Wageningen University & Research.
 *
 * Use and duplication is permitted under the terms of the
 * BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
 *
 */

#ifndef ARRAY_H
#define ARRAY_H

#include <array>
#include <vector>
#include <algorithm>
#include <iostream>

template<int N>
inline std::array<int, N> calc_strides(const std::array<int, N>& dims)
{
    std::array<int, N> strides;
    strides[0] = 1;
    for (int i=1; i<N; ++i)
        strides[i] = strides[i-1]*dims[i-1];

    return strides;
}

template<int N>
inline int calc_index(
        const std::array<int, N>& indices,
        const std::array<int, N>& strides,
        const std::array<int, N>& offsets)
{
    int sum = 0;
    for (int i=0; i<N; ++i)
        sum += (indices[i]-offsets[i]-1)*strides[i];

    return sum;
}

template<int N>
inline std::array<int, N> calc_indices(
        int index, const std::array<int, N>& strides, const std::array<int, N>& offsets)
{
    std::array<int, N> indices;

    for (int i=N-1; i>=1; --i)
    {
        indices[i] = index / strides[i];
        index %= strides[i];
    }
    indices[0] = index;

    for (int i=0; i<N; ++i)
        indices[i] += offsets[i] + 1;

    return indices;
}

template<int N>
inline int product(const std::array<int, N>& array)
{
    int product = array[0];
    for (int i=1; i<N; ++i)
        product *= array[i];

    return product;
}

template<typename T, int N>
class Array
{
    public:
        // Create an empty array, without dimensions.
        Array() :
            dims({}),
            ncells(0)
        {}

        // Create an array of zeros with given dimensions.
        Array(const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(ncells),
            strides(calc_strides<N>(dims)),
            offsets({})
        {}

        // Create an array from copying the contents of an std::vector.
        Array(const std::vector<T>& data, const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(data),
            strides(calc_strides<N>(dims)),
            offsets({})
        {} // CvH Do we need to size check data?

        // Create an array from moving the contents of an std::vector.
        Array(std::vector<T>&& data, const std::array<int, N>& dims) :
            dims(dims),
            ncells(product<N>(dims)),
            data(data),
            strides(calc_strides<N>(dims)),
            offsets({})
        {} // CvH Do we need to size check data?

        // Define the default copy constructor and assignment operator.
        // Array(const Array<T, N>&) = default;
        // Array<T,N>& operator=(const Array<T, N>&) = default; // CvH does this one need empty checking?

        // Array(Array<T, N>&& array) :
        //     dims(std::exchange(array.dims, {})),
        //     ncells(std::exchange(array.ncells, 0)),
        //     data(std::move(array.data)),
        //     strides(std::exchange(array.strides, {})),
        //     offsets(std::exchange(array.offsets, {}))
        // {}

        inline void set_offsets(const std::array<int, N>& offsets)
        {
            this->offsets = offsets;
        }

        inline std::array<int, N> get_dims() const { return dims; }

        inline void set_dims(const std::array<int, N>& dims)
        {
            if (this->ncells != 0)
                throw std::runtime_error("Only arrays of size 0 can be resized");

            this->dims = dims;
            ncells = product<N>(dims);
            data.resize(ncells);
            strides = calc_strides<N>(dims);
            offsets = {};
        }

        inline std::vector<T>& v() { return data; }
        inline const std::vector<T>& v() const { return data; }

        inline T* ptr() { return data.data(); }
        inline const T* ptr() const { return data.data(); }

        inline int size() const { return ncells; }

        // inline std::array<int, N> find_indices(const T& value) const
        // {
        //     int pos = std::find(data.begin(), data.end(), value) - data.begin();
        //     return calc_indices<N>(pos, strides, offsets);
        // }

        inline T max() const
        {
            return *std::max_element(data.begin(), data.end());
        }

        inline T min() const
        {
            return *std::min_element(data.begin(), data.end());
        }

        inline void operator=(std::vector<T>&& data)
        {
            // CvH check size.
            this->data = data;
        }

        inline T& operator()(const std::array<int, N>& indices)
        {
            const int index = calc_index<N>(indices, strides, offsets);
            return data[index];
        }

        inline T operator()(const std::array<int, N>& indices) const
        {
            const int index = calc_index<N>(indices, strides, offsets);
            return data[index];
        }

        inline int dim(const int i) const { return dims[i-1]; }
        inline bool is_empty() const { return ncells == 0; }

        inline Array<T, N> subset(
                const std::array<std::pair<int, int>, N> ranges) const
        {
            // Calculate the dimension sizes based on the range.
            std::array<int, N> subdims;
            std::array<bool, N> do_spread;

            for (int i=0; i<N; ++i)
            {
                subdims[i] = ranges[i].second - ranges[i].first + 1;
                // CvH how flexible / tolerant are we?
                do_spread[i] = (dims[i] == 1);
            }

            // Create the array and fill it with the subset.
            Array<T, N> a_sub(subdims);
            for (int i=0; i<a_sub.ncells; ++i)
            {
                std::array<int, N> index;
                int ic = i;
                for (int n=N-1; n>0; --n)
                {
                    index[n] = do_spread[n] ? 1 : ic / a_sub.strides[n] + ranges[n].first;
                    ic %= a_sub.strides[n];
                }
                index[0] = do_spread[0] ? 1 : ic + ranges[0].first;
                a_sub.data[i] = (*this)(index);
            }

            return a_sub;
        }

        inline void fill(const T value)
        {
            std::fill(data.begin(), data.end(), value);
        }

    private:
        std::array<int, N> dims;
        int ncells;
        std::vector<T> data;
        std::array<int, N> strides;
        std::array<int, N> offsets;
};

template<typename T, int N>
bool any_vals_outside(const Array<T, N>& array, const T lower_limit, const T upper_limit)
{
    return std::any_of(
            array.v().begin(),
            array.v().end(),
            [lower_limit, upper_limit](T val){ return (val < lower_limit) || (val > upper_limit); });
}

template<typename T, int N>
bool any_vals_less_than(const Array<T, N>& array, const T lower_limit)
{
    return std::any_of(
            array.v().begin(),
            array.v().end(),
            [lower_limit](T val){ return (val < lower_limit); });
}
#endif
