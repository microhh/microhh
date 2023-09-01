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

#ifndef NETCDF_INTERFACE_H
#define NETCDF_INTERFACE_H

#include <map>
#include <vector>
#include <netcdf.h>

enum class Netcdf_mode { Create, Read, Write };

class Master;
class Netcdf_handle;
class Netcdf_group;

template<typename T>
class Netcdf_variable
{
    public:
        Netcdf_variable(Master&, Netcdf_handle&, const int, const std::vector<int>&);

        // Do not allow copying of netcdf variable.
        Netcdf_variable(const Netcdf_variable&) = delete;
        Netcdf_variable& operator=(const Netcdf_variable&) = delete;

        // Enable the default move constructor to move initialized variables into containers.
        Netcdf_variable(Netcdf_variable&&) = default;

        void insert(const std::vector<T>&, const std::vector<int>);
        void insert(const std::vector<T>&, const std::vector<int>, const std::vector<int>);
        void insert(const T, const std::vector<int>);

        const std::vector<int> get_dim_sizes() { return dim_sizes; }

        void add_attribute(const std::string&, const std::string&);
        void add_attribute(const std::string&, const double);
        void add_attribute(const std::string&, const float);

    private:
        Master& master;
        Netcdf_handle& nc_handle;
        const int var_id;
        const std::vector<int> dim_sizes;
};

class Netcdf_handle
{
    public:
        Netcdf_handle(Master&);
        virtual ~Netcdf_handle() = default;

        // Do not allow copying or moving of handle.
        Netcdf_handle(const Netcdf_handle&) = delete;
        Netcdf_handle& operator=(const Netcdf_handle&) = delete;

        void add_dimension(const std::string&, const int dim_size = NC_UNLIMITED);

        Netcdf_group& add_group(const std::string&);
        Netcdf_group& get_group(const std::string&);

        int get_dimension_size(const std::string&);

        std::map<std::string, int> get_variable_dimensions(const std::string&);

        bool variable_exists(const std::string&);
        bool group_exists(const std::string&);

        template<typename T>
        Netcdf_variable<T> add_variable(
                const std::string&,
                const std::vector<std::string>&);

        template<typename T>
        T get_variable(
            const std::string&);

        template<typename T>
        std::vector<T> get_variable(
            const std::string&,
            const std::vector<int>&);

        template<typename T>
        void get_variable(
                std::vector<T>&,
                const std::string&,
                const std::vector<int>&,
                const std::vector<int>&,
                const bool required_read = true);

        template<typename T>
        void insert(
                const std::vector<T>&,
                const int var_id,
                const std::vector<int>&,
                const std::vector<int>&);

        template<typename T>
        void insert(
                const T,
                const int var_id,
                const std::vector<int>&,
                const std::vector<int>&);

        void add_attribute(
                const std::string&,
                const std::string&,
                const int);

        void add_attribute(
                const std::string&,
                const double,
                const int);

        void add_attribute(
                const std::string&,
                const float,
                const int);

        virtual int get_dim_id(const std::string&) = 0;

    protected:
        Master& master;
        Netcdf_handle* parent;
        int mpiid_to_write;
        int ncid;
        int root_ncid;
        std::map<std::string, int> dims;
        std::map<std::string, Netcdf_group> groups;
        int record_counter;
};

class Netcdf_file : public Netcdf_handle
{
    public:
        Netcdf_file(Master&, const std::string&, Netcdf_mode, const int mpiid_to_write_int=0);
        virtual ~Netcdf_file();

        // Do not allow copying or moving of file
        Netcdf_file(const Netcdf_file&) = delete;
        Netcdf_file& operator=(const Netcdf_file&) = delete;

        int get_dim_id(const std::string&);

        void sync();
};

class Netcdf_group : public Netcdf_handle
{
    public:
        Netcdf_group(
                Master&, Netcdf_handle*,
                const int, const int, const int);

        // Do not allow copying or moving of groups.
        Netcdf_group(const Netcdf_group&) = delete;
        Netcdf_group& operator=(const Netcdf_group&) = delete;

        int get_dim_id(const std::string&);
};
#endif
