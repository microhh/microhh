#ifndef NETCDF_INTERFACE_H
#define NETCDF_INTERFACE_H

#include <iostream>
#include <map>
#include <vector>
#include <netcdf.h>

enum class Netcdf_mode { Create, Read, Write };

class Master;
class Netcdf_handle;
class Netcdf_group;

class Netcdf_variable
{
    public:
        Netcdf_variable(Master&, Netcdf_handle&, const int, const std::vector<int>&);
        void insert(const std::vector<double>&, const std::vector<int>);
        void insert(const std::vector<double>&, const std::vector<int>, const std::vector<int>);
        void insert(const double, const std::vector<int>);

    private:
        Master& master;
        Netcdf_handle& nc_file;
        const int var_id;
        const std::vector<int> dim_sizes;
};

class Netcdf_handle
{
    public:
        Netcdf_handle(Master&);
        void add_dimension(const std::string&, const int dim_size = NC_UNLIMITED);

        Netcdf_group add_group(const std::string&);
        Netcdf_group get_group(const std::string&);

        Netcdf_variable add_variable(
                const std::string&,
                const std::vector<std::string>);

        template<typename T>
        void get_variable(
                std::vector<T>&,
                const std::string&,
                const std::vector<int>&,
                const std::vector<int>&);

        void insert(
                const std::vector<double>&,
                const int var_id,
                const std::vector<int>&,
                const std::vector<int>&);

        void insert(
                const double,
                const int var_id,
                const std::vector<int>&,
                const std::vector<int>&);

    protected:
        Master& master;
        int ncid;
        int root_ncid;
        std::map<std::string, int> dims;
        int record_counter;
};

class Netcdf_file : public Netcdf_handle
{
    public:
        Netcdf_file(Master&, const std::string&, Netcdf_mode);
        ~Netcdf_file();
};

class Netcdf_group : public Netcdf_handle
{
    public:
        Netcdf_group(Master&, const int, const int);
};
#endif
