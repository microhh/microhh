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

#ifndef HDF5_INTERFACE_H
#define HDF5_INTERFACE_H


#include <hdf5.h>
#include <hdf5_hl.h>  // Needed for dimension names/scales

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <map>

template<typename T>
hid_t get_hdf5_type()
{
    if constexpr (std::is_same_v<T, float>)
        return H5T_NATIVE_FLOAT;
    else if constexpr (std::is_same_v<T, double>)
        return H5T_NATIVE_DOUBLE;
    else if constexpr (std::is_same_v<T, int>)
        return H5T_NATIVE_INT;
    else
        throw std::runtime_error("Invalid datatype for HDF5!");
}

enum class Hdf5_mode {Read, Write, Read_write};

class Hdf5_file
{
    public:
        Hdf5_file(const std::string& filename, Hdf5_mode m) : mode(m)
        {
            if (mode == Hdf5_mode::Read)
            {
                std::ifstream test(filename);
                if (!test.good())
                    throw std::runtime_error("File does not exist: " + filename);
                test.close();

                file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

                if (file_id < 0)
                    throw std::runtime_error("Failed to open file: " + filename);
            }
            else if (mode == Hdf5_mode::Write)
            {
                std::ifstream test(filename);
                if (test.good())
                {
                    test.close();
                    throw std::runtime_error("File already exists: " + filename);
                }

                file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

                if (file_id < 0)
                    throw std::runtime_error("Failed to create file: " + filename);
            }
            else
            {
                file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

                if (file_id < 0)
                    throw std::runtime_error("Failed to open file for read/write: " + filename);
            }
        }

        ~Hdf5_file()
        {
            close();
        }

        void close()
        {
            if (file_id >= 0)
            {
                H5Fclose(file_id);
                file_id = -1;
            }
        }

        hid_t get_id() const
        {
            return file_id;
        }

        void add_dimension(const std::string& name, int size)
        {
            dimensions[name] = size;
        }

        int get_dimension_size(const std::string& name) const
        {
            auto it = dimensions.find(name);
            if (it == dimensions.end())
                throw std::runtime_error("Dimension not found: " + name);
            return it->second;
        }

        std::vector<std::string> list_datasets() const
        {
            std::vector<std::string> names;
            H5Literate(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, nullptr, collect_names, &names);
            return names;
        }

    private:
        hid_t file_id;
        Hdf5_mode mode;
        std::map<std::string, int> dimensions;

        static herr_t collect_names(
                hid_t loc_id,
                const char *name,
                const H5L_info_t *info,
                void *data)
        {
            auto* names = static_cast<std::vector<std::string>*>(data);

            H5O_info_t obj_info;

            #if H5_VERSION_GE(1, 14, 0)
                H5Oget_info_by_name3(loc_id, name, &obj_info, H5O_INFO_BASIC, H5P_DEFAULT);
            #else
                // Boomer workstation.
                H5Oget_info_by_name(loc_id, name, &obj_info, H5P_DEFAULT);
            #endif

            if (obj_info.type == H5O_TYPE_DATASET)
                names->push_back(name);

            return 0;
        }
};

template<typename TF>
class Hdf5_variable
{
    public:
        // Open existing variable
        Hdf5_variable(const Hdf5_file& file, const std::string& path)
        {
            dataset_id = H5Dopen2(file.get_id(), path.c_str(), H5P_DEFAULT);
            if (dataset_id < 0)
                throw std::runtime_error("Failed to open dataset: " + path);
        }

        Hdf5_variable(
                const Hdf5_file& file,
                const std::string& path,
                const std::vector<std::string>& dim_names)
        {
            // Get dimension sizes
            std::vector<hsize_t> h5_dims;
            std::vector<hsize_t> h5_maxdims;

            for (const auto& dim_name : dim_names)
            {
                int size = file.get_dimension_size(dim_name);
                if (size == -1)
                {
                    h5_dims.push_back(0);
                    h5_maxdims.push_back(H5S_UNLIMITED);
                }
                else
                {
                    h5_dims.push_back(size);
                    h5_maxdims.push_back(size);
                }
            }

            hid_t space = H5Screate_simple(h5_dims.size(), h5_dims.data(), h5_maxdims.data());
            if (space < 0)
                throw std::runtime_error("Failed to create dataspace");

            hid_t plist = H5P_DEFAULT;
            bool has_unlimited = false;
            for (const auto& dim_name : dim_names)
            {
                if (file.get_dimension_size(dim_name) == -1)
                {
                    has_unlimited = true;
                    break;
                }
            }

            if (has_unlimited)
            {
                plist = H5Pcreate(H5P_DATASET_CREATE);
                std::vector<hsize_t> chunk_dims = h5_dims;
                if (chunk_dims[0] == 0) chunk_dims[0] = 1;
                H5Pset_chunk(plist, chunk_dims.size(), chunk_dims.data());
            }

            hid_t h5_type = get_hdf5_type<TF>();
            dataset_id = H5Dcreate2(
                    file.get_id(),
                    path.c_str(),
                    h5_type,
                    space,
                    H5P_DEFAULT,
                    plist,
                    H5P_DEFAULT);

            if (plist != H5P_DEFAULT)
                H5Pclose(plist);
            H5Sclose(space);

            if (dataset_id < 0)
                throw std::runtime_error("Failed to create dataset: " + path);

            for (size_t i = 0; i < dim_names.size(); i++)
            {
                hid_t dim_id = H5Dopen2(file.get_id(), dim_names[i].c_str(), H5P_DEFAULT);
                if (dim_id >= 0)
                {
                    H5DSattach_scale(dataset_id, dim_id, i);
                    H5Dclose(dim_id);
                }
            }
        }

        ~Hdf5_variable()
        {
            if (dataset_id >= 0)
                H5Dclose(dataset_id);
        }

        std::vector<size_t> get_dims() const
        {
            hid_t space = H5Dget_space(dataset_id);
            if (space < 0)
                throw std::runtime_error("Failed to get dataspace");

            int ndims = H5Sget_simple_extent_ndims(space);
            std::vector<hsize_t> dims_tmp(ndims);
            H5Sget_simple_extent_dims(space, dims_tmp.data(), nullptr);
            H5Sclose(space);

            std::vector<size_t> dims(ndims);
            for (int i = 0; i < ndims; i++)
                dims[i] = dims_tmp[i];

            return dims;
        }

        std::vector<TF> read() const
        {
            auto dims = get_dims();

            size_t count = 1;
            for (size_t d : dims)
                count *= d;

            std::vector<TF> data(count);

            hid_t h5_type = get_hdf5_type<TF>();
            herr_t status = H5Dread(
                    dataset_id,
                    h5_type,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    data.data());

            if (status < 0)
                throw std::runtime_error("Failed to read dataset");

            return data;
        }

        // Read hyperslab (partial dataset)
        // start: starting indices for each dimension
        // count: number of elements to read in each dimension
        std::vector<TF> read(const std::vector<size_t>& start,
                             const std::vector<size_t>& count) const
        {
            std::vector<hsize_t> h5_start(start.begin(), start.end());
            std::vector<hsize_t> h5_count(count.begin(), count.end());

            hid_t file_space = H5Dget_space(dataset_id);
            if (file_space < 0)
                throw std::runtime_error("Failed to get dataspace");

            herr_t status = H5Sselect_hyperslab(
                file_space, H5S_SELECT_SET,
                h5_start.data(), nullptr,
                h5_count.data(), nullptr);

            if (status < 0)
            {
                H5Sclose(file_space);
                throw std::runtime_error("Failed to select hyperslab");
            }

            hid_t mem_space = H5Screate_simple(h5_count.size(), h5_count.data(), nullptr);
            if (mem_space < 0)
            {
                H5Sclose(file_space);
                throw std::runtime_error("Failed to create memory dataspace");
            }

            size_t total = 1;
            for (size_t c : count)
                total *= c;

            std::vector<TF> data(total);

            hid_t h5_type = get_hdf5_type<TF>();
            status = H5Dread(
                dataset_id,
                h5_type,
                mem_space,
                file_space,
                H5P_DEFAULT,
                data.data());

            H5Sclose(mem_space);
            H5Sclose(file_space);

            if (status < 0)
                throw std::runtime_error("Failed to read hyperslab");

            return data;
        }

        void insert(const std::vector<TF>& data)
        {
            hid_t h5_type = get_hdf5_type<TF>();
            herr_t status = H5Dwrite(
                    dataset_id,
                    h5_type,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    data.data());

            if (status < 0)
                throw std::runtime_error("Failed to write dataset");
        }

        // Write hyperslab (partial dataset)
        // start: starting indices for each dimension
        // count: number of elements to write in each dimension
        void insert(const std::vector<size_t>& start,
                   const std::vector<size_t>& count,
                   const std::vector<TF>& data)
        {
            // Convert to hsize_t
            std::vector<hsize_t> h5_start(start.begin(), start.end());
            std::vector<hsize_t> h5_count(count.begin(), count.end());

            // Get the dataspace of the dataset
            hid_t file_space = H5Dget_space(dataset_id);
            if (file_space < 0)
                throw std::runtime_error("Failed to get dataspace");

            // Check if we need to extend the dataset (for unlimited dimensions)
            int ndims = H5Sget_simple_extent_ndims(file_space);
            std::vector<hsize_t> current_dims(ndims);
            std::vector<hsize_t> max_dims(ndims);
            H5Sget_simple_extent_dims(file_space, current_dims.data(), max_dims.data());

            bool need_extend = false;
            std::vector<hsize_t> new_dims = current_dims;
            for (int i = 0; i < ndims; i++)
            {
                hsize_t required_size = h5_start[i] + h5_count[i];
                if (required_size > current_dims[i])
                {
                    if (max_dims[i] != H5S_UNLIMITED)
                        throw std::runtime_error("Cannot extend fixed dimension");
                    new_dims[i] = required_size;
                    need_extend = true;
                }
            }

            if (need_extend)
            {
                herr_t status = H5Dset_extent(dataset_id, new_dims.data());
                if (status < 0)
                    throw std::runtime_error("Failed to extend dataset");

                H5Sclose(file_space);
                file_space = H5Dget_space(dataset_id);
            }

            herr_t status = H5Sselect_hyperslab(
                file_space, H5S_SELECT_SET,
                h5_start.data(), nullptr,
                h5_count.data(), nullptr);

            if (status < 0)
            {
                H5Sclose(file_space);
                throw std::runtime_error("Failed to select hyperslab");
            }

            hid_t mem_space = H5Screate_simple(h5_count.size(), h5_count.data(), nullptr);
            if (mem_space < 0)
            {
                H5Sclose(file_space);
                throw std::runtime_error("Failed to create memory dataspace");
            }

            hid_t h5_type = get_hdf5_type<TF>();
            status = H5Dwrite(
                dataset_id,
                h5_type,
                mem_space,
                file_space,
                H5P_DEFAULT,
                data.data());

            H5Sclose(mem_space);
            H5Sclose(file_space);

            if (status < 0)
                throw std::runtime_error("Failed to write hyperslab");
        }

    private:
        hid_t dataset_id;
};
#endif
