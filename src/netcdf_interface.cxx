#include <iostream>
#include <map>
#include <vector>
#include <netcdf.h>

#include "netcdf_interface.h"
#include "master.h"

namespace
{
    void nc_throw(const int return_value)
    {
        std::string error(nc_strerror(return_value));
        throw std::runtime_error(error);
    }

    void nc_check(const int return_value)
    {
        if (return_value != NC_NOERR)
            nc_throw(return_value);
    }
}

Netcdf_file::Netcdf_file(Master& master, const std::string& name, Netcdf_mode mode) :
    Netcdf_handle(master)
{
    if (mode == Netcdf_mode::Create)
        nc_check( nc_create(name.c_str(), NC_NOCLOBBER | NC_NETCDF4, &ncid) );
    else if (mode == Netcdf_mode::Write)
        nc_check( nc_open(name.c_str(), NC_WRITE | NC_NETCDF4, &ncid) );
    else if (mode == Netcdf_mode::Read)
        nc_check( nc_open(name.c_str(), NC_NOWRITE | NC_NETCDF4, &ncid) );

    root_ncid = ncid;

    if (mode == Netcdf_mode::Create)
        nc_check( nc_enddef(root_ncid) );
}

Netcdf_file::~Netcdf_file()
{
    nc_check( nc_close(ncid) );
}

void Netcdf_handle::add_dimension(const std::string& dim_name, const size_t dim_size)
{
    nc_check( nc_redef(root_ncid) );

    int dim_id;
    int def_out = nc_def_dim(ncid, dim_name.c_str(), dim_size, &dim_id);

    // Dimension is written or already exists.
    if (def_out == NC_NOERR)
        dims.emplace(dim_name, dim_id);
    else if (def_out == NC_ENAMEINUSE)
    {}
    // Error.
    else
        nc_throw(def_out);

    nc_check( nc_enddef(root_ncid) );
}

Netcdf_variable Netcdf_handle::add_variable(
        const std::string& var_name,
        const std::vector<std::string> dim_names)
{
    nc_check ( nc_redef(root_ncid) );

    int ndims = dim_names.size();
    std::vector<int> dim_ids;
    for (const std::string& dim_name : dim_names)
        dim_ids.push_back(dims.at(dim_name));

    int var_id;
    nc_check( nc_def_var(ncid, var_name.c_str(), NC_DOUBLE, ndims, dim_ids.data(), &var_id) );
    nc_check( nc_enddef(root_ncid) );

    std::vector<size_t> dim_sizes;
    for (size_t i=0; i<dim_ids.size(); ++i)
    {
        const int dim_id = dim_ids[i];

        size_t dim_len;
        nc_check( nc_inq_dimlen(ncid, dim_id, &dim_len) );

        if (dim_len == NC_UNLIMITED)
        {
            if (i == 0)
                dim_sizes.push_back(1);
            else
                throw std::runtime_error("Only the outermost dimension is allowed to be an NC_UNLIMITED dimension");
        }
        else
            dim_sizes.push_back(dim_len);
    }

    return Netcdf_variable(*this, var_id, dim_sizes);
}

Netcdf_handle::Netcdf_handle(Master& master) :
    master(master), record_counter(0)
{}

void Netcdf_handle::insert(
        const std::vector<double>& values,
        const int var_id,
        const std::vector<size_t>& i_start,
        const std::vector<size_t>& i_count)
{
    // CvH: Add proper size checking.
    nc_check( nc_put_vara_double(ncid, var_id, i_start.data(), i_count.data(), values.data()) );
}

void Netcdf_handle::insert(
        const double value,
        const int var_id,
        const std::vector<size_t>& i_start,
        const std::vector<size_t>& i_count)
{
    // CvH: Add proper size checking.
    nc_check( nc_put_vara_double(ncid, var_id, i_start.data(), i_count.data(), &value) );
}

Netcdf_group Netcdf_handle::add_group(const std::string& name)
{
    int group_ncid;

    nc_check( nc_redef(root_ncid) );
    nc_check( nc_def_grp(ncid, name.c_str(), &group_ncid) );
    nc_check( nc_enddef(root_ncid) );

    return Netcdf_group(master, group_ncid, root_ncid);
}

Netcdf_group Netcdf_handle::get_group(const std::string& name)
{
    int group_ncid;
    nc_check( nc_inq_ncid(ncid, name.c_str(), &group_ncid) );

    return Netcdf_group(master, group_ncid, root_ncid);
}


Netcdf_group::Netcdf_group(Master& master, int ncid_in, int root_ncid_in) :
    Netcdf_handle(master)
{
    ncid = ncid_in;
    root_ncid = root_ncid_in;
}

void Netcdf_handle::get_variable(
        std::vector<double>& values,
        const std::string& name,
        const std::vector<size_t>& i_start,
        const std::vector<size_t>& i_count)
{
    int var_id;
    nc_check( nc_inq_varid(ncid, name.c_str(), &var_id) );

    // CvH: Add check if the vector is large enough.
    nc_check( nc_get_vara_double(ncid, var_id, i_start.data(), i_count.data(), values.data()) );
}

Netcdf_variable::Netcdf_variable(Netcdf_handle& nc_file, const int var_id, const std::vector<size_t>& dim_sizes) :
    nc_file(nc_file), var_id(var_id), dim_sizes(dim_sizes)
{}

void Netcdf_variable::insert(const std::vector<double>& values, const std::vector<size_t> i_start)
{
    nc_file.insert(values, var_id, i_start, dim_sizes);
}

void Netcdf_variable::insert(
        const std::vector<double>& values,
        const std::vector<size_t> i_start,
        const std::vector<size_t> i_count)
{
    nc_file.insert(values, var_id, i_start, i_count);
}

void Netcdf_variable::insert(const double value, const std::vector<size_t> i_start)
{
    nc_file.insert(value, var_id, i_start, dim_sizes);
}
