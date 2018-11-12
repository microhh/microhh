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
    int nc_check_code = 0;

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
    int nc_check_code = 0;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_close(ncid);
    nc_check(master, nc_check_code);
}

void Netcdf_handle::add_dimension(const std::string& dim_name, const int dim_size)
{
    int nc_check_code = 0;

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
    int nc_check_code = 0;

    int var_id = -1;
    std::vector<int> dim_sizes;

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

        size_t dim_len = 0;

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
        const std::vector<int>& i_start,
        const std::vector<int>& i_count)
{
    const std::vector<size_t> i_start_size_t (i_start.begin(), i_start.end());
    const std::vector<size_t> i_count_size_t (i_count.begin(), i_count.end());

    int nc_check_code = 0;

    // CvH: Add proper size checking.
    if (master.get_mpiid() == 0)
        nc_check_code = nc_put_vara_double(ncid, var_id, i_start_size_t.data(), i_count_size_t.data(), values.data());
    nc_check(master, nc_check_code);
}

void Netcdf_handle::insert(
        const double value,
        const int var_id,
        const std::vector<int>& i_start,
        const std::vector<int>& i_count)
{
    const std::vector<size_t> i_start_size_t (i_start.begin(), i_start.end());
    const std::vector<size_t> i_count_size_t (i_count.begin(), i_count.end());

    int nc_check_code = 0;

    // CvH: Add proper size checking.
    if (master.get_mpiid() == 0)
        nc_check_code = nc_put_vara_double(ncid, var_id, i_start_size_t.data(), i_count_size_t.data(), &value);
    nc_check(master, nc_check_code);
}

Netcdf_group Netcdf_handle::add_group(const std::string& name)
{
    int group_ncid = -1;
    int nc_check_code = 0;

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
    int nc_check_code = 0;

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

std::map<std::string, int> Netcdf_handle::get_variable_dimensions(const std::string& name)
{
    int nc_check_code = 0;
    int var_id;

    if (master.get_mpiid() == 0)
        nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
    nc_check(master, nc_check_code);

    int ndims;
    int dimids[NC_MAX_VAR_DIMS];

    if (master.get_mpiid() == 0)
        nc_check_code = nc_inq_var(ncid, var_id, NULL, NULL, &ndims, dimids, NULL);
    nc_check(master, nc_check_code);

    // Broadcast ndims
    master.broadcast(&ndims, 1);

    std::map<std::string, int> dims;

    for (int n=0; n<ndims; ++n)
    {
        char dim_name[NC_MAX_NAME+1];
        size_t dim_length_size_t;

        nc_check_code = nc_inq_dim(ncid, dimids[n], dim_name, &dim_length_size_t);
        nc_check(master, nc_check_code);

        int dim_length = dim_length_size_t;
        master.broadcast(&dim_length, 1);

        // Broadcast the entire buffer to avoid broadcasting of length.
        master.broadcast(dim_name, NC_MAX_NAME+1);

        dims.emplace(std::string(dim_name), dim_length);
    }

    return dims;
}

namespace
{
    template<typename TF>
    int nc_get_vara_wrapper(
            int, int, const std::vector<size_t>&, const std::vector<size_t>&, std::vector<TF>&);

    template<>
    int nc_get_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, std::vector<double>& values)
    {
        return nc_get_vara_double(ncid, var_id, start.data(), count.data(), values.data());
    }

    template<>
    int nc_get_vara_wrapper(
            int ncid, int var_id, const std::vector<size_t>& start, const std::vector<size_t>& count, std::vector<float>& values)
    {
        return nc_get_vara_float(ncid, var_id, start.data(), count.data(), values.data());
    }
}

template<typename TF>
void Netcdf_handle::get_variable(
        std::vector<TF>& values,
        const std::string& name,
        const std::vector<int>& i_start,
        const std::vector<int>& i_count)
{
    std::string message = "Retrieving from NetCDF: " + name;
    master.print_message(message);

    const std::vector<size_t> i_start_size_t (i_start.begin(), i_start.end());
    const std::vector<size_t> i_count_size_t (i_count.begin(), i_count.end());

    int nc_check_code = 0;
    int var_id;

    bool zero_fill = false;
    try
    {
        if (master.get_mpiid() == 0)
            nc_check_code = nc_inq_varid(ncid, name.c_str(), &var_id);
        nc_check(master, nc_check_code);
    }
    catch (std::runtime_error& e)
    {
        std::string warning = "Netcdf variable " + name + " not found, filling with zeros";
        master.print_warning(warning);
        zero_fill = true;
    }

    // CvH: Add check if the vector is large enough.
    int total_count = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());
    // If the vector is long enough, it can be copied. We assume that this routine does NOT resize vectors.
    master.broadcast(&total_count, 1);

    if (zero_fill)
    {
        std::fill(values.begin(), values.begin() + total_count, 0);
    }
    else
    {
        if (master.get_mpiid() == 0)
            nc_check_code = nc_get_vara_wrapper(ncid, var_id, i_start_size_t, i_count_size_t, values);
        nc_check(master, nc_check_code);
        master.broadcast(values.data(), total_count);
    }
}

// Variable does not communicate with NetCDF library directly.
Netcdf_variable::Netcdf_variable(Master& master, Netcdf_handle& nc_file, const int var_id, const std::vector<int>& dim_sizes) :
    master(master), nc_file(nc_file), var_id(var_id), dim_sizes(dim_sizes)
{}

void Netcdf_variable::insert(const std::vector<double>& values, const std::vector<int> i_start)
{
    nc_file.insert(values, var_id, i_start, dim_sizes);
}

void Netcdf_variable::insert(
        const std::vector<double>& values,
        const std::vector<int> i_start,
        const std::vector<int> i_count)
{
    nc_file.insert(values, var_id, i_start, i_count);
}

void Netcdf_variable::insert(const double value, const std::vector<int> i_start)
{
    nc_file.insert(value, var_id, i_start, dim_sizes);
}

template void Netcdf_handle::get_variable<double>(std::vector<double>&, const std::string&, const std::vector<int>&, const std::vector<int>&);
template void Netcdf_handle::get_variable<float> (std::vector<float> &, const std::string&, const std::vector<int>&, const std::vector<int>&);
