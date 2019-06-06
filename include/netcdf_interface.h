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
        Netcdf_variable(const Netcdf_variable&) = default;
        void insert(const std::vector<T>&, const std::vector<int>);
        void insert(const std::vector<T>&, const std::vector<int>, const std::vector<int>);
        void insert(const T, const std::vector<int>);
        const std::vector<int> get_dim_sizes() { return dim_sizes; }
        void add_attribute(const std::string&, const std::string&);
        void add_attribute(const std::string&, const double);
        void add_attribute(const std::string&, const float);

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

        int get_dimension_size(const std::string&);

        std::map<std::string, int> get_variable_dimensions(const std::string&);

        bool variable_exists(const std::string&);
        bool group_exists(const std::string&);

        template<typename T>
        Netcdf_variable<T> add_variable(
                const std::string&,
                const std::vector<std::string>);

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
                const std::vector<int>&);

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

    protected:
        Master& master;
        int mpiid_to_write;
        int ncid;
        int root_ncid;
        std::map<std::string, int> dims;
        int record_counter;
};

class Netcdf_file : public Netcdf_handle
{
    public:
        Netcdf_file(Master&, const std::string&, Netcdf_mode, const int mpiid_to_write_int=0);
        ~Netcdf_file();

        void sync();
};

class Netcdf_group : public Netcdf_handle
{
    public:
        Netcdf_group(Master&, const int, const int, const int);
};
#endif
