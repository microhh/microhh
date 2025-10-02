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
    template<typename TF>
    void read_coordinate(
        hid_t file_id, const char* dset_name,
        TF* buffer,
        hsize_t start,
        hsize_t count)
    {
        // Open dataset.
        hid_t dset = H5Dopen(file_id, dset_name, H5P_DEFAULT);
        hid_t filespace = H5Dget_space(dset);

        // Select hyperslab.
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);

        // Local memory layout is simple; continous 1D array of size `count`.
        hid_t memspace = H5Screate_simple(1, &count, NULL);

        // Setup collective IO where all tasks participate.
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        // Read data.
        H5Dread(dset, get_hdf5_type<TF>(), memspace, filespace, plist_id, buffer);

        // Cleanup!
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Dclose(dset);
    }


    template <typename TF>
    void read_particles_parallel(
        const std::string& filename,
        std::vector<TF>& x,
        std::vector<TF>& y,
        std::vector<TF>& z,
        const int mpiid,
        const int nprocs)
    {
        // Open file with parallel HDF5.
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id);
        H5Pclose(plist_id);

        // Get total number of particles.
        hid_t dset_x = H5Dopen(file_id, "/x", H5P_DEFAULT);
        hid_t dspace = H5Dget_space(dset_x);
        hsize_t n_total;
        H5Sget_simple_extent_dims(dspace, &n_total, NULL);
        H5Sclose(dspace);

        // Each task has `ceil(n_total / nprocs)` particles (except last task, see next code block).
        const int np_per_task = std::ceil(TF(n_total) / nprocs);
        hsize_t start = mpiid * np_per_task;
        hsize_t count = np_per_task;

        // Last MPI task might have less if `n_total % nprocs != 0`.
        if (start + count > n_total)
            count = n_total - start;

        //std::cout << "n_total=" << n_total << ", n_per_task=" << np_per_task << ", mpiid=" << mpiid << ", start=" << start << ", count=" << count << std::endl;

        // Resize local vectors.
        // Use local size without buffer; these particles are only
        // temporarely on this task and shipped elsewhwere soon.
        x.resize(count);
        y.resize(count);
        z.resize(count);

        // Read each coordinate with hyperslab selection
        read_coordinate(file_id, "/x", x.data(), start, count);
        read_coordinate(file_id, "/y", y.data(), start, count);
        read_coordinate(file_id, "/z", z.data(), start, count);

        // Cleanup!
        H5Dclose(dset_x);
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
    struct Particle
    {
        TF x, y, z;
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
    MPI_Datatype create_particle_type()
    {
        MPI_Datatype particle_type;
        MPI_Type_contiguous(3, get_mpi_type<TF>(), &particle_type);
        MPI_Type_commit(&particle_type);
        return particle_type;
    }


    template<typename TF>
    void distribute_particles(
        std::vector<TF>& x_out,
        std::vector<TF>& y_out,
        std::vector<TF>& z_out,
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

        // Calculate target `mpiid` for each particle
        std::vector<int> target_rank(np_local);
        for (int n=0; n<np_local; ++n)
        {
            const int mpicoordx = int(x_in[n] / xsize_sub);
            const int mpicoordy = int(y_in[n] / ysize_sub);

            target_rank[n] = master.calc_mpiid(mpicoordx, mpicoordy);
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

        // Calculate send/receive offsets (cumulative sum send/recv counts)/
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
        std::vector<Particle<TF>> particles_send(total_send);
        std::vector<int> current_offset = send_offsets;

        for (int n=0; n < np_local; ++n)
        {
            const int rank = target_rank[n];
            const int pos = current_offset[rank];
            current_offset[rank] += 1;

            particles_send[pos].x = x_in[n];
            particles_send[pos].y = y_in[n];
            particles_send[pos].z = z_in[n];
        }

        // Allocate receive buffer.
        std::vector<Particle<TF>> particles_recv(total_recv);

        // Exchange particles.
        MPI_Datatype particle_type = create_particle_type<TF>();

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
        x_out.resize(total_recv);
        y_out.resize(total_recv);
        z_out.resize(total_recv);

        for (int n=0; n<total_recv; ++n)
        {
            x_out[n] = particles_recv[n].x;
            y_out[n] = particles_recv[n].y;
            z_out[n] = particles_recv[n].z;
        }

        // DEBUG.
        //for (int n=0; n<total_recv; ++n)
        //    std::cout << "mpiidx/y= " << md.mpicoordx << "/" << md.mpicoordy << " has x=" << x_out[n] << ", y=" << y_out[n] << std::endl;
    }
}
#endif
