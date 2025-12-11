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

#ifndef PARTICLE_LAGRANGIAN_IO_H
#define PARTICLE_LAGRANGIAN_IO_H

#include <hdf5.h>
#include "hdf5_interface.h"

namespace Particle_lagrangian_io
{
    template<typename T>
    void read_coordinate(
        hid_t file_id,
        const char* dset_name,
        T* buffer,
        hsize_t start,
        hsize_t count)
    {
        // Open dataset.
        hid_t dset = H5Dopen(file_id, dset_name, H5P_DEFAULT);
        hid_t filespace = H5Dget_space(dset);
        hid_t memspace;

        if (count > 0)
        {
            // Select hyperslab.
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);

            // Local memory layout is simple; continous 1D array of size `count`.
            memspace = H5Screate_simple(1, &count, NULL);
        }
        else
        {
            // No particles to read - select empty hyperslab.
            H5Sselect_none(filespace);
            memspace = H5Scopy(filespace);
            H5Sselect_none(memspace);
        }

        // Setup collective IO where all tasks participate.
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // Read data (all tasks must call this even if count=0).
        H5Dread(dset, get_hdf5_type<T>(), memspace, filespace, plist_id, buffer);

        // Cleanup!
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Dclose(dset);
    }


    template<typename T>
    void write_coordinate(
        hid_t file_id,
        const char* dset_name,
        const T* buffer,
        hsize_t time_idx,
        hsize_t start,
        hsize_t count,
        hsize_t n_total)
    {
        // Open dataset.
        hid_t dset = H5Dopen(file_id, dset_name, H5P_DEFAULT);
        hid_t filespace = H5Dget_space(dset);
        hid_t memspace;

        if (count > 0)
        {
            // Select hyperslab in file: [time_idx, start:start+count]
            hsize_t file_start[2] = {time_idx, start};
            hsize_t file_count[2] = {1, count};
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

            // Local memory layout is simple; continuous 1D array of size `count`.
            memspace = H5Screate_simple(1, &count, NULL);
        }
        else
        {
            // No particles to write - select empty hyperslab.
            H5Sselect_none(filespace);
            memspace = H5Scopy(filespace);
            H5Sselect_none(memspace);
        }

        // Setup collective IO where all tasks participate.
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // Write data (all tasks must call this even if count=0).
        H5Dwrite(dset, get_hdf5_type<T>(), memspace, filespace, plist_id, buffer);

        // Cleanup!
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Dclose(dset);
    }


    template <typename TF>
    void read_particles_parallel(
        const std::string& filename,
        std::vector<int>& uid,
        std::vector<TF>& x,
        std::vector<TF>& y,
        std::vector<TF>& z,
        int& n_particles,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        // Open file with parallel HDF5.
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, md.commxy, MPI_INFO_NULL);
        H5Pset_all_coll_metadata_ops(plist_id, true);
        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
        H5Pclose(plist_id);

        // Get total number of particles.
        hid_t dset_x = H5Dopen(file_id, "/x", H5P_DEFAULT);
        hid_t dspace = H5Dget_space(dset_x);
        hsize_t n_total;
        H5Sget_simple_extent_dims(dspace, &n_total, NULL);
        H5Sclose(dspace);
        H5Dclose(dset_x);

        n_particles = static_cast<int>(n_total);

        // Each task has `ceil(n_total / nprocs)` particles (except last task, see next code block).
        const int np_per_task = std::ceil(TF(n_total) / md.nprocs);
        hsize_t start = md.mpiid * np_per_task;
        hsize_t count = np_per_task;

        // Handle case where start >= n_total (more tasks than particles per task).
        if (start >= n_total)
            count = 0;
        else if (start + count > n_total)
            count = n_total - start;

        //std::cout << "n_total=" << n_total << ", n_per_task=" << np_per_task << ", mpiid=" << mpiid << ", start=" << start << ", count=" << count << std::endl;

        // Resize local vectors.
        // Use local size without buffer; these particles are only
        // temporarely on this task and shipped elsewhwere soon.
        uid.resize(count);
        x.resize(count);
        y.resize(count);
        z.resize(count);

        // Read each coordinate with hyperslab selection
        read_coordinate(file_id, "/particle_id", uid.data(), start, count);
        read_coordinate(file_id, "/x", x.data(), start, count);
        read_coordinate(file_id, "/y", y.data(), start, count);
        read_coordinate(file_id, "/z", z.data(), start, count);

        // Cleanup!
        H5Fclose(file_id);
    }


    template <typename TF>
    void read_particles_serial(
        const std::string& filename,
        std::vector<int>& uid,
        std::vector<TF>& x,
        std::vector<TF>& y,
        std::vector<TF>& z)
    {
        Hdf5_file h5_file_in(filename, Hdf5_mode::Read);

        Hdf5_variable<int> var_uid(h5_file_in, "particle_id");
        Hdf5_variable<TF> var_x(h5_file_in, "x");
        Hdf5_variable<TF> var_y(h5_file_in, "y");
        Hdf5_variable<TF> var_z(h5_file_in, "z");

        uid = var_uid.read();
        x = var_x.read();
        y = var_y.read();
        z = var_z.read();
    }


    template<typename TF>
    struct Particle_io
    {
        int uid;
        TF x, y, z;     // Location.
    };


    template<typename T>
    MPI_Datatype get_mpi_type()
    {
        if constexpr (std::is_same_v<T, float>)
            return MPI_FLOAT;
        else if constexpr (std::is_same_v<T, double>)
            return MPI_DOUBLE;
        else
            throw std::runtime_error("Invalid float type for MPI!");
    }


    template<typename TF>
    bool particle_uid_compare(const Particle_io<TF>& a, const Particle_io<TF>& b)
    {
        return a.uid < b.uid;
    }


    template<typename TF>
    MPI_Datatype create_particle_type_io()
    {
        // Create MPI type specific for particle I/O, which typically uses less
        // elements per particle than the neighbour-neighbour MPI communication.
        MPI_Datatype particle_type;

        // 1=uid, 3=number of particle properties (currently x,y,z).
        int blocklengths[2] = {1, 3};
        MPI_Aint displacements[2];
        MPI_Datatype types[2] = {MPI_INT, get_mpi_type<TF>()};

        Particle_io<TF> particle;
        MPI_Aint base_address;
        MPI_Get_address(&particle, &base_address);
        MPI_Get_address(&particle.uid, &displacements[0]);
        MPI_Get_address(&particle.x, &displacements[1]);

        displacements[0] = MPI_Aint_diff(displacements[0], base_address);
        displacements[1] = MPI_Aint_diff(displacements[1], base_address);

        MPI_Type_create_struct(2, blocklengths, displacements, types, &particle_type);
        MPI_Type_commit(&particle_type);

        return particle_type;
    }


    template<typename TF>
    void distribute_particles(
        std::vector<int>& uid_out,
        std::vector<TF>& x_out,
        std::vector<TF>& y_out,
        std::vector<TF>& z_out,
        std::vector<int>& uid_in,
        std::vector<TF>& x_in,
        std::vector<TF>& y_in,
        std::vector<TF>& z_in,
        const TF xsize,
        const TF ysize,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        const TF xsize_sub = xsize / md.npx;
        const TF ysize_sub = ysize / md.npy;

        const int np_local = x_in.size();

        // Calculate target `mpiid` for each particle.
        std::vector<int> target_rank(np_local);
        for (int n=0; n<np_local; ++n)
        {
            const int mpicoordx = int(x_in[n] / xsize_sub);
            const int mpicoordy = int(y_in[n] / ysize_sub);

            target_rank[n] = master.get_mpiid(mpicoordx, mpicoordy);
        }

        // Count particles going to each processor.
        std::vector<int> send_counts(md.nprocs, 0);
        for (int n=0; n<np_local; ++n)
            send_counts[target_rank[n]] += 1;

        // Exchange send counts across all tasks.
        std::vector<int> recv_counts(md.nprocs);
        const int size = 1;
        MPI_Alltoall(
            send_counts.data(),
            size,
            MPI_INT,
            recv_counts.data(),
            size,
            MPI_INT,
            md.commxy);

        // Calculate send/receive offsets (cumulative sum send/recv counts).
        // If e.g. send_counts = {3,5,2,4}, then
        //         send_offsets = {0,3,8,10}.
        std::vector<int> send_offsets(md.nprocs, 0);
        std::vector<int> recv_offsets(md.nprocs, 0);
        for (int i=1; i<md.nprocs; ++i)
        {
            send_offsets[i] = send_offsets[i-1] + send_counts[i-1];
            recv_offsets[i] = recv_offsets[i-1] + recv_counts[i-1];
        }

        const int total_send = send_offsets[md.nprocs-1] + send_counts[md.nprocs-1];
        const int total_recv = recv_offsets[md.nprocs-1] + recv_counts[md.nprocs-1];

        // Pack particles into vector of Particle structs. This should make it easier
        // to add other properties like velocity or mass at a later point.
        std::vector<Particle_io<TF>> particles_send(total_send);
        std::vector<int> current_offset = send_offsets;

        for (int n=0; n < np_local; ++n)
        {
            const int rank = target_rank[n];
            const int pos = current_offset[rank];
            current_offset[rank] += 1;

            particles_send[pos].uid = uid_in[n];
            particles_send[pos].x = x_in[n];
            particles_send[pos].y = y_in[n];
            particles_send[pos].z = z_in[n];
        }

        // Allocate receive buffer.
        std::vector<Particle_io<TF>> particles_recv(total_recv);

        // Exchange particles.
        MPI_Datatype particle_type = create_particle_type_io<TF>();

        MPI_Alltoallv(
            particles_send.data(),
            send_counts.data(),
            send_offsets.data(),
            particle_type,
            particles_recv.data(),
            recv_counts.data(),
            recv_offsets.data(),
            particle_type,
            md.commxy);

        MPI_Type_free(&particle_type);

        // Unpack Particle structs in local vectors.
        uid_out.resize(total_recv);
        x_out.resize(total_recv);
        y_out.resize(total_recv);
        z_out.resize(total_recv);

        for (int n=0; n<total_recv; ++n)
        {
            uid_out[n] = particles_recv[n].uid;
            x_out[n] = particles_recv[n].x;
            y_out[n] = particles_recv[n].y;
            z_out[n] = particles_recv[n].z;
        }

        // DEBUG.
        //for (int n=0; n<total_recv; ++n)
        //    std::cout << "mpiidx/y= " << md.mpicoordx << "/" << md.mpicoordy << " has x=" << x_out[n] << ", y=" << y_out[n] << std::endl;
    }


    template<typename TF>
    void create_particle_dump(
        const std::string& filename,
        hid_t& file_id,
        const int n_particles,
        const int szip_compression,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        // Throw error if file exists.
        int file_exists = 0;
        if (md.mpiid == 0)
        {
            std::ifstream test(filename);
            file_exists = test.good() ? 1 : 0;
            test.close();
        }
        MPI_Bcast(&file_exists, 1, MPI_INT, 0, md.commxy);

        if (file_exists)
        {
            std::string error_msg = "ERROR: Particle output file already exists: " + filename;
            throw std::runtime_error(error_msg);
        }

        // Create file with parallel HDF5.
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, md.commxy, MPI_INFO_NULL);
        H5Pset_all_coll_metadata_ops(plist_id, true);
        H5Pset_coll_metadata_write(plist_id, true);
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);

        // Create unlimited time dimension.
        hsize_t time_dims[1] = {0};
        hsize_t time_maxdims[1] = {H5S_UNLIMITED};
        hid_t dspace_time = H5Screate_simple(1, time_dims, time_maxdims);

        hid_t plist_create = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        H5Pset_chunk(plist_create, 1, chunk_dims);

        hid_t dset_time = H5Dcreate(
            file_id, "/time", H5T_NATIVE_DOUBLE, dspace_time,
            H5P_DEFAULT, plist_create, H5P_DEFAULT);

        // Make time a dimension scale
        H5DSset_scale(dset_time, "time");

        H5Dclose(dset_time);
        H5Pclose(plist_create);
        H5Sclose(dspace_time);

        // Create particle_id dimension scale dataset.
        hsize_t particle_dims[1] = {static_cast<hsize_t>(n_particles)};
        hid_t dspace_particle_id = H5Screate_simple(1, particle_dims, particle_dims);
        hid_t dset_particle_id = H5Dcreate(
            file_id, "/particle_id", H5T_NATIVE_INT, dspace_particle_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Make particle_id a dimension scale
        H5DSset_scale(dset_particle_id, "particle_id");

        // Write particle IDs (0, 1, 2, ..., n_particles-1) as coordinate values
        std::vector<int> particle_ids(n_particles);
        for (int i = 0; i < n_particles; ++i)
            particle_ids[i] = i;

        hid_t plist_write = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_write, H5FD_MPIO_COLLECTIVE);
        H5Dwrite(dset_particle_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_write, particle_ids.data());
        H5Pclose(plist_write);

        H5Dclose(dset_particle_id);
        H5Sclose(dspace_particle_id);

        // Create particle datasets with unlimited time dimension.
        hsize_t dims[2] = {0, static_cast<hsize_t>(n_particles)};
        hsize_t maxdims[2] = {H5S_UNLIMITED, static_cast<hsize_t>(n_particles)};
        hid_t dspace = H5Screate_simple(2, dims, maxdims);

        hid_t plist_create_2d = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims_2d[2] = {1, static_cast<hsize_t>(n_particles)};
        H5Pset_chunk(plist_create_2d, 2, chunk_dims_2d);

        // Apply SZIP compression if requested.
        if (szip_compression > 0)
            H5Pset_szip(plist_create_2d, H5_SZIP_NN_OPTION_MASK, szip_compression);

        hid_t h5_type = get_hdf5_type<TF>();

        hid_t dset_x = H5Dcreate(
            file_id, "/x", h5_type, dspace,
            H5P_DEFAULT, plist_create_2d, H5P_DEFAULT);
        hid_t dset_y = H5Dcreate(
            file_id, "/y", h5_type, dspace,
            H5P_DEFAULT, plist_create_2d, H5P_DEFAULT);
        hid_t dset_z = H5Dcreate(
            file_id, "/z", h5_type, dspace,
            H5P_DEFAULT, plist_create_2d, H5P_DEFAULT);

        // Attach dimension scales
        dset_time = H5Dopen(file_id, "/time", H5P_DEFAULT);
        dset_particle_id = H5Dopen(file_id, "/particle_id", H5P_DEFAULT);

        H5DSattach_scale(dset_x, dset_time, 0);
        H5DSattach_scale(dset_x, dset_particle_id, 1);
        H5DSattach_scale(dset_y, dset_time, 0);
        H5DSattach_scale(dset_y, dset_particle_id, 1);
        H5DSattach_scale(dset_z, dset_time, 0);
        H5DSattach_scale(dset_z, dset_particle_id, 1);

        H5Dclose(dset_time);
        H5Dclose(dset_particle_id);
        H5Dclose(dset_x);
        H5Dclose(dset_y);
        H5Dclose(dset_z);
        H5Pclose(plist_create_2d);
        H5Sclose(dspace);

        // Flush to ensure datasets are written to file.
        H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    }


    template<typename TF>
    void create_particle_restart(
        const std::string& filename,
        hid_t& file_id,
        const int n_particles,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        // Throw error if file exists.
        int file_exists = 0;
        if (md.mpiid == 0)
        {
            std::ifstream test(filename);
            file_exists = test.good() ? 1 : 0;
            test.close();
        }
        MPI_Bcast(&file_exists, 1, MPI_INT, 0, md.commxy);

        if (file_exists)
        {
            std::string error_msg = "ERROR: Particle output file already exists: " + filename;
            throw std::runtime_error(error_msg);
        }

        // Create file with parallel HDF5.
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, md.commxy, MPI_INFO_NULL);
        H5Pset_all_coll_metadata_ops(plist_id, true);
        H5Pset_coll_metadata_write(plist_id, true);
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);

        // Create particle_id dimension scale dataset.
        hsize_t particle_dims[1] = {static_cast<hsize_t>(n_particles)};
        hid_t dspace_particle_id = H5Screate_simple(1, particle_dims, particle_dims);
        hid_t dset_particle_id = H5Dcreate(
            file_id, "/particle_id", H5T_NATIVE_INT, dspace_particle_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Make particle_id a dimension scale
        H5DSset_scale(dset_particle_id, "particle_id");

        // Write particle IDs (0, 1, 2, ..., n_particles-1) as coordinate values
        std::vector<int> particle_ids(n_particles);
        for (int i = 0; i < n_particles; ++i)
            particle_ids[i] = i;

        hid_t plist_write = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_write, H5FD_MPIO_COLLECTIVE);
        H5Dwrite(dset_particle_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_write, particle_ids.data());
        H5Pclose(plist_write);

        H5Dclose(dset_particle_id);
        H5Sclose(dspace_particle_id);

        // Create particle datasets without time dimension (single time step).
        hsize_t dims[1] = {static_cast<hsize_t>(n_particles)};
        hid_t dspace = H5Screate_simple(1, dims, dims);

        hid_t plist_create = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {static_cast<hsize_t>(n_particles)};
        H5Pset_chunk(plist_create, 1, chunk_dims);

        hid_t h5_type = get_hdf5_type<TF>();

        hid_t dset_x = H5Dcreate(
            file_id, "/x", h5_type, dspace,
            H5P_DEFAULT, plist_create, H5P_DEFAULT);
        hid_t dset_y = H5Dcreate(
            file_id, "/y", h5_type, dspace,
            H5P_DEFAULT, plist_create, H5P_DEFAULT);
        hid_t dset_z = H5Dcreate(
            file_id, "/z", h5_type, dspace,
            H5P_DEFAULT, plist_create, H5P_DEFAULT);

        // Attach dimension scale
        dset_particle_id = H5Dopen(file_id, "/particle_id", H5P_DEFAULT);

        H5DSattach_scale(dset_x, dset_particle_id, 0);
        H5DSattach_scale(dset_y, dset_particle_id, 0);
        H5DSattach_scale(dset_z, dset_particle_id, 0);

        H5Dclose(dset_particle_id);
        H5Dclose(dset_x);
        H5Dclose(dset_y);
        H5Dclose(dset_z);
        H5Pclose(plist_create);
        H5Sclose(dspace);

        // Flush to ensure datasets are written to file.
        H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    }


    template<typename TF>
    void write_particles_restart(
        hid_t file_id,
        const std::vector<int>& uid,
        const std::vector<TF>& x,
        const std::vector<TF>& y,
        const std::vector<TF>& z,
        const int n_particles,
        const int mpiid,
        const int nprocs)
    {
        // Calculate hyperslab parameters for this MPI task.
        const int np_per_task = std::ceil(TF(n_particles) / nprocs);
        hsize_t start = static_cast<hsize_t>(mpiid * np_per_task);
        hsize_t count = static_cast<hsize_t>(np_per_task);

        // Handle case where start >= n_particles (more tasks than particles per task).
        if (start >= static_cast<hsize_t>(n_particles))
            count = 0;
        else if (start + count > static_cast<hsize_t>(n_particles))
            count = static_cast<hsize_t>(n_particles) - start;

        // Open datasets.
        hid_t dset_x = H5Dopen(file_id, "/x", H5P_DEFAULT);
        hid_t dset_y = H5Dopen(file_id, "/y", H5P_DEFAULT);
        hid_t dset_z = H5Dopen(file_id, "/z", H5P_DEFAULT);

        hid_t filespace_x = H5Dget_space(dset_x);
        hid_t filespace_y = H5Dget_space(dset_y);
        hid_t filespace_z = H5Dget_space(dset_z);

        hid_t memspace;

        if (count > 0)
        {
            // Select hyperslab in file: [start:start+count]
            hsize_t file_start[1] = {start};
            hsize_t file_count[1] = {count};
            H5Sselect_hyperslab(filespace_x, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
            H5Sselect_hyperslab(filespace_y, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
            H5Sselect_hyperslab(filespace_z, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

            // Local memory layout is simple; continuous 1D array of size `count`.
            memspace = H5Screate_simple(1, &count, NULL);
        }
        else
        {
            // No particles to write - select empty hyperslab.
            H5Sselect_none(filespace_x);
            H5Sselect_none(filespace_y);
            H5Sselect_none(filespace_z);
            memspace = H5Scopy(filespace_x);
            H5Sselect_none(memspace);
        }

        // Setup collective IO where all tasks participate.
        hid_t plist_xfer = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_xfer, H5FD_MPIO_COLLECTIVE);

        // Write data (all tasks must call this even if count=0).
        H5Dwrite(dset_x, get_hdf5_type<TF>(), memspace, filespace_x, plist_xfer, x.data());
        H5Dwrite(dset_y, get_hdf5_type<TF>(), memspace, filespace_y, plist_xfer, y.data());
        H5Dwrite(dset_z, get_hdf5_type<TF>(), memspace, filespace_z, plist_xfer, z.data());

        // Cleanup!
        H5Pclose(plist_xfer);
        H5Sclose(memspace);
        H5Sclose(filespace_x);
        H5Sclose(filespace_y);
        H5Sclose(filespace_z);
        H5Dclose(dset_x);
        H5Dclose(dset_y);
        H5Dclose(dset_z);

        // Flush to disk to ensure data is written.
        H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    }


    template<typename TF>
    void write_particles_parallel(
        hid_t file_id,
        const std::vector<int>& uid,
        const std::vector<TF>& x,
        const std::vector<TF>& y,
        const std::vector<TF>& z,
        const double time,
        const int n_particles,
        const int mpiid,
        const int nprocs)
    {
        // Get current time dimension size_dump
        hid_t dset_time = H5Dopen(file_id, "/time", H5P_DEFAULT);
        hid_t dspace_time = H5Dget_space(dset_time);
        hsize_t time_idx;
        H5Sget_simple_extent_dims(dspace_time, &time_idx, NULL);
        H5Sclose(dspace_time);
        H5Dclose(dset_time);

        // Extend datasets to add new time step.
        hsize_t new_dims_time[1] = {time_idx + 1};
        hsize_t new_dims[2] = {time_idx + 1, static_cast<hsize_t>(n_particles)};

        // Extend time dataset.
        dset_time = H5Dopen(file_id, "/time", H5P_DEFAULT);
        H5Dset_extent(dset_time, new_dims_time);

        // Write time value.
        hid_t dspace_time_write = H5Dget_space(dset_time);
        hsize_t time_start[1] = {time_idx};
        hsize_t time_count[1] = {1};
        H5Sselect_hyperslab(dspace_time_write, H5S_SELECT_SET, time_start, NULL, time_count, NULL);
        hid_t memspace_time = H5Screate_simple(1, time_count, NULL);

        hid_t plist_xfer = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_xfer, H5FD_MPIO_COLLECTIVE);

        H5Dwrite(dset_time, H5T_NATIVE_DOUBLE, memspace_time, dspace_time_write, plist_xfer, &time);

        H5Pclose(plist_xfer);
        H5Sclose(memspace_time);
        H5Sclose(dspace_time_write);
        H5Dclose(dset_time);

        // Extend particle datasets.
        hid_t dset_x = H5Dopen(file_id, "/x", H5P_DEFAULT);
        hid_t dset_y = H5Dopen(file_id, "/y", H5P_DEFAULT);
        hid_t dset_z = H5Dopen(file_id, "/z", H5P_DEFAULT);

        H5Dset_extent(dset_x, new_dims);
        H5Dset_extent(dset_y, new_dims);
        H5Dset_extent(dset_z, new_dims);

        H5Dclose(dset_x);
        H5Dclose(dset_y);
        H5Dclose(dset_z);

        // Calculate hyperslab parameters for this MPI task.
        const int np_per_task = std::ceil(TF(n_particles) / nprocs);
        hsize_t start = static_cast<hsize_t>(mpiid * np_per_task);
        hsize_t count = static_cast<hsize_t>(np_per_task);

        // Handle case where start >= n_particles (more tasks than particles per task).
        if (start >= static_cast<hsize_t>(n_particles))
            count = 0;
        else if (start + count > static_cast<hsize_t>(n_particles))
            count = static_cast<hsize_t>(n_particles) - start;

        // Write particle coordinates (write_coordinate opens/closes datasets internally).
        write_coordinate(file_id, "/x", x.data(), time_idx, start, count, n_particles);
        write_coordinate(file_id, "/y", y.data(), time_idx, start, count, n_particles);
        write_coordinate(file_id, "/z", z.data(), time_idx, start, count, n_particles);

        // Flush to disk to ensure data is written.
        H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    }


    template<typename TF>
    void gather_particles(
        std::vector<int>& uid_out,
        std::vector<TF>& x_out,
        std::vector<TF>& y_out,
        std::vector<TF>& z_out,
        const std::vector<int>& uid_in,
        const std::vector<TF>& x_in,
        const std::vector<TF>& y_in,
        const std::vector<TF>& z_in,
        const int n_particles,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        // Each task originally has `ceil(n_particles / nprocs)` particles.
        const int np_per_task = std::ceil(TF(n_particles) / md.nprocs);

        const int np_local = uid_in.size();

        // Calculate which original MPI task each particle belongs to based on uid.
        std::vector<int> target_rank(np_local);
        for (int n=0; n<np_local; ++n)
            target_rank[n] = uid_in[n] / np_per_task;

        // Count particles going to each processor.
        std::vector<int> send_counts(md.nprocs, 0);
        for (int n=0; n<np_local; ++n)
            send_counts[target_rank[n]] += 1;

        // Exchange send counts across all tasks.
        std::vector<int> recv_counts(md.nprocs);
        const int size = 1;
        MPI_Alltoall(
            send_counts.data(),
            size,
            MPI_INT,
            recv_counts.data(),
            size,
            MPI_INT,
            md.commxy);

        // Calculate send/receive offsets.
        std::vector<int> send_offsets(md.nprocs, 0);
        std::vector<int> recv_offsets(md.nprocs, 0);
        for (int i=1; i<md.nprocs; ++i)
        {
            send_offsets[i] = send_offsets[i-1] + send_counts[i-1];
            recv_offsets[i] = recv_offsets[i-1] + recv_counts[i-1];
        }

        const int total_send = send_offsets[md.nprocs-1] + send_counts[md.nprocs-1];
        const int total_recv = recv_offsets[md.nprocs-1] + recv_counts[md.nprocs-1];

        // Pack particles into send buffer.
        std::vector<Particle_io<TF>> particles_send(total_send);
        std::vector<int> current_offset = send_offsets;

        for (int n=0; n < np_local; ++n)
        {
            const int rank = target_rank[n];
            const int pos = current_offset[rank];
            current_offset[rank] += 1;

            particles_send[pos].uid = uid_in[n];
            particles_send[pos].x = x_in[n];
            particles_send[pos].y = y_in[n];
            particles_send[pos].z = z_in[n];
        }

        // Allocate receive buffer.
        std::vector<Particle_io<TF>> particles_recv(total_recv);

        // Exchange particles.
        MPI_Datatype particle_type = create_particle_type_io<TF>();

        MPI_Alltoallv(
            particles_send.data(),
            send_counts.data(),
            send_offsets.data(),
            particle_type,
            particles_recv.data(),
            recv_counts.data(),
            recv_offsets.data(),
            particle_type,
            md.commxy);

        MPI_Type_free(&particle_type);

        // Sort received particles by uid to restore original order.
        std::sort(particles_recv.begin(), particles_recv.end(), particle_uid_compare<TF>);

        // Unpack Particle structs in local vectors.
        uid_out.resize(total_recv);
        x_out.resize(total_recv);
        y_out.resize(total_recv);
        z_out.resize(total_recv);

        for (int n=0; n<total_recv; ++n)
        {
            uid_out[n] = particles_recv[n].uid;
            x_out[n] = particles_recv[n].x;
            y_out[n] = particles_recv[n].y;
            z_out[n] = particles_recv[n].z;
        }
    }
}
#endif
