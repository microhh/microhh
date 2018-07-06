#include <iostream>
#include <map>
#include <numeric>
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

    void nc_check(Master& master, int return_value)
    {
        master.broadcast(&return_value, 1);
        if (return_value != NC_NOERR)
            nc_throw(return_value);
    }
}

Netcdf_file::Netcdf_file(Master& master, const std::string& name, Netcdf_mode mode) :
    Netcdf_handle(master)
{
    int nc_check_code;

    if (master.get_mpiid() == 0)
    {
        if (mode == Netcdf_mode::Create)
            nc_check_code = nc_create(name.c_str(), NC_NOCLOBBER | NC_NETCDF4, &ncid);
        else if (mode == Netcdf_mode::Write)
            nc_check_code = nc_open(name.c_str(), NC_WRITE | NC_NETCDF4, &ncid);
        else if (mode == Netcdf_mode::Read)
            nc_check_code = nc_open(name.c_str(), NC_NOWRITE | NC_NETCDF4, &ncid);
    }

    nc_check(master, nc_check_code);

    root_ncid = ncid;

    if (master.get_mpiid() == 0)
    {
        if (mode == Netcdf_mode::Create)
            nc_check_code =  nc_enddef(root_ncid);
    }

    nc_check(master, nc_check_code);
}

Netcdf_file::~Netcdf_file()
{
    int nc_check_code;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_close(ncid);
    nc_check(master, nc_check_code);
}

void Netcdf_handle::add_dimension(const std::string& dim_name, const size_t dim_size)
{
    int nc_check_code;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_redef(root_ncid);
    nc_check(master, nc_check_code);

    int dim_id;
    int def_out;

    if (master.get_mpiid() == 0)
        def_out = nc_def_dim(ncid, dim_name.c_str(), dim_size, &dim_id);

    master.broadcast(&def_out, 1);

    // Dimension is written or already exists.
    if (def_out == NC_NOERR)
        dims.emplace(dim_name, dim_id);
    else if (def_out == NC_ENAMEINUSE)
    {}
    // Error.
    else
        nc_throw(def_out);

    if (master.get_mpiid() == 0)
        nc_check_code = nc_enddef(root_ncid);

    nc_check(master, nc_check_code);
}

Netcdf_variable Netcdf_handle::add_variable(
        const std::string& var_name,
        const std::vector<std::string> dim_names)
{
    int nc_check_code;

    int var_id = -1;
    std::vector<size_t> dim_sizes;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_redef(root_ncid);
    nc_check(master, nc_check_code);

    int ndims = dim_names.size();
    std::vector<int> dim_ids;

    for (const std::string& dim_name : dim_names)
        dim_ids.push_back(dims.at(dim_name));

    if (master.get_mpiid() == 0)
        nc_check_code = nc_def_var(ncid, var_name.c_str(), NC_DOUBLE, ndims, dim_ids.data(), &var_id);
    nc_check(master, nc_check_code);

    if (master.get_mpiid() == 0)
        nc_check_code = nc_enddef(root_ncid);
    nc_check(master, nc_check_code);

    // Broadcast the dim_ids size of the main process to run the for loop on all processes.
    int dim_ids_size = dim_ids.size();
    master.broadcast(&dim_ids_size, 1);

    for (int i=0; i<dim_ids_size; ++i)
    {
        const int dim_id = dim_ids.at(i);

        size_t dim_len;

        if (master.get_mpiid() == 0)
            nc_check_code = nc_inq_dimlen(ncid, dim_id, &dim_len);
        nc_check(master, nc_check_code);

        int dim_len_int = static_cast<int>(dim_len);
        master.broadcast(&dim_len_int, 1);

        if (dim_len_int == NC_UNLIMITED)
        {
            if (i == 0)
                dim_sizes.push_back(1);
            else
                throw std::runtime_error("Only the leftmost dimension is allowed to be an NC_UNLIMITED dimension");
        }
        else
            dim_sizes.push_back(dim_len);
    }

    return Netcdf_variable(master, *this, var_id, dim_sizes);
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
    int nc_check_code;

    // CvH: Add proper size checking.
    if (master.get_mpiid() == 0)
        nc_check_code = nc_put_vara_double(ncid, var_id, i_start.data(), i_count.data(), values.data());
    nc_check(master, nc_check_code);
}

void Netcdf_handle::insert(
        const double value,
        const int var_id,
        const std::vector<size_t>& i_start,
        const std::vector<size_t>& i_count)
{
    int nc_check_code;

    // CvH: Add proper size checking.
    if (master.get_mpiid() == 0)
        nc_check_code = nc_put_vara_double(ncid, var_id, i_start.data(), i_count.data(), &value);
    nc_check(master, nc_check_code);
}

Netcdf_group Netcdf_handle::add_group(const std::string& name)
{
    int group_ncid = -1;
    int nc_check_code;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_redef(root_ncid);
    nc_check(master, nc_check_code);

    if (master.get_mpiid() == 0)
        nc_check_code = nc_def_grp(ncid, name.c_str(), &group_ncid);
    nc_check(master, nc_check_code);

    if (master.get_mpiid() == 0)
        nc_check_code = nc_enddef(root_ncid);
    nc_check(master, nc_check_code);

    return Netcdf_group(master, group_ncid, root_ncid);
}

Netcdf_group Netcdf_handle::get_group(const std::string& name)
{
    int group_ncid = -1;
    int nc_check_code;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_inq_ncid(ncid, name.c_str(), &group_ncid);
    nc_check(master, nc_check_code);

    return Netcdf_group(master, group_ncid, root_ncid);
}


Netcdf_group::Netcdf_group(Master& master, int ncid_in, int root_ncid_in) :
    Netcdf_handle(master)
{
    ncid = ncid_in;
    root_ncid = root_ncid_in;
}

template<>
void Netcdf_handle::get_variable(
        std::vector<double>& values,

        const std::string& name,
        const std::vector<size_t>& i_start,
        const std::vector<size_t>& i_count)
{
    int nc_check_code;
    int var_id;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
    nc_check(master, nc_check_code);

    // CvH: Add check if the vector is large enough.
    if (master.get_mpiid() == 0)
        nc_check_code = nc_get_vara_double(ncid, var_id, i_start.data(), i_count.data(), values.data());
    nc_check(master, nc_check_code);

    // If the vector is long enough, it can be copied. We assume that this routine does NOT resize vectors.
    int total_count = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());
    master.broadcast(&total_count, 1);
    master.broadcast(values.data(), total_count);
}

template<>
void Netcdf_handle::get_variable(
        std::vector<float>& values,
        const std::string& name,
        const std::vector<size_t>& i_start,
        const std::vector<size_t>& i_count)
{
    int nc_check_code;

    int var_id;
    if (master.get_mpiid() == 0)
        nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
    nc_check(master, nc_check_code);

    // CvH: This needs to be removed once the stats is properly templated.
    std::vector<double>values_double(values.size());

    // CvH: Add check if the vector is large enough.
    if (master.get_mpiid() == 0)
    {
        std::vector<double>values_double(values.size());
        nc_check_code = nc_get_vara_double(ncid, var_id, i_start.data(), i_count.data(), values_double.data());
    }
    nc_check(master, nc_check_code);

    if (master.get_mpiid() == 0)
        std::copy(values_double.begin(), values_double.end(), values.begin());
}

// Variable does not communicate with NetCDF library directly.
Netcdf_variable::Netcdf_variable(Master& master, Netcdf_handle& nc_file, const int var_id, const std::vector<size_t>& dim_sizes) :
    master(master), nc_file(nc_file), var_id(var_id), dim_sizes(dim_sizes)
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

// template void Netcdf_handle::get_variable<double>(std::vector<double>&, const std::string&, const std::vector<size_t>&, const std::vector<size_t>&);
// template void Netcdf_handle::get_variable<float> (std::vector<float> &, const std::string&, const std::vector<size_t>&, const std::vector<size_t>&);
