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

#include <algorithm>
#include <numeric>

#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5_interface.h"

namespace Particle_lagrangian_io
{
    // Negative return values of HDF5 functions (`hid_t` and `herr_t`) indicate failure.
    inline hid_t h5_check(const hid_t status, const std::string& msg)
    {
        if (status < 0)
            throw std::runtime_error("HDF5 error: " + msg);
        return status;
    }


    // Each task has `ceil(n_particles / nprocs)` particles, except the last task(s).
    inline int get_np_per_task(const int n_particles, const int nprocs)
    {
        return (n_particles + nprocs - 1) / nprocs;
    }


    // Contiguous block [start, start+count) of particles that is read/written by this task.
    inline void get_io_block(
        hsize_t& start,
        hsize_t& count,
        const int n_particles,
        const int mpiid,
        const int nprocs)
    {
        const int np_per_task = get_np_per_task(n_particles, nprocs);
        start = static_cast<hsize_t>(mpiid) * np_per_task;
        count = np_per_task;

        // Handle case where start >= n_particles (more tasks than particles per task).
        if (start >= static_cast<hsize_t>(n_particles))
            count = 0;
        else if (start + count > static_cast<hsize_t>(n_particles))
            count = static_cast<hsize_t>(n_particles) - start;
    }


    // File access property list; parallel HDF5 with collective metadata operations for MPI builds.
    inline hid_t create_fapl(Master& master)
    {
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

        #ifdef USEMPI
        auto& md = master.get_MPI_data();
        h5_check(H5Pset_fapl_mpio(plist_id, md.commxy, MPI_INFO_NULL), "setting MPI-IO file driver");
        h5_check(H5Pset_all_coll_metadata_ops(plist_id, true), "setting collective metadata reads");
        h5_check(H5Pset_coll_metadata_write(plist_id, true), "setting collective metadata writes");
        #endif

        return plist_id;
    }


    // Dataset transfer property list; collective IO for MPI builds.
    inline hid_t create_dxpl()
    {
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);

        #ifdef USEMPI
        h5_check(H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE), "setting collective IO");
        #endif

        return plist_id;
    }


    template<typename T>
    void read_coordinate(
        hid_t file_id,
        const char* dset_name,
        T* buffer,
        hsize_t start,
        hsize_t count)
    {
        // Open dataset.
        hid_t dset = h5_check(H5Dopen(file_id, dset_name, H5P_DEFAULT), "opening dataset " + std::string(dset_name));
        hid_t filespace = H5Dget_space(dset);
        hid_t memspace;

        if (count > 0)
        {
            // Select hyperslab.
            h5_check(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL),
                    "selecting hyperslab in " + std::string(dset_name));

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
        hid_t plist_id = create_dxpl();

        // Read data (all tasks must call this even if count=0).
        h5_check(H5Dread(dset, get_hdf5_type<T>(), memspace, filespace, plist_id, buffer),
                "reading dataset " + std::string(dset_name));

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
        hid_t dset = h5_check(H5Dopen(file_id, dset_name, H5P_DEFAULT), "opening dataset " + std::string(dset_name));
        hid_t filespace = H5Dget_space(dset);
        hid_t memspace;

        if (count > 0)
        {
            // Select hyperslab in file: [time_idx, start:start+count]
            hsize_t file_start[2] = {time_idx, start};
            hsize_t file_count[2] = {1, count};
            h5_check(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL),
                    "selecting hyperslab in " + std::string(dset_name));

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
        hid_t plist_id = create_dxpl();

        // Write data (all tasks must call this even if count=0).
        h5_check(H5Dwrite(dset, get_hdf5_type<T>(), memspace, filespace, plist_id, buffer),
                "writing dataset " + std::string(dset_name));

        // Cleanup!
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Dclose(dset);
    }


    template<typename T>
    void write_coordinate_1d(
        hid_t file_id,
        const char* dset_name,
        const T* buffer,
        hsize_t start,
        hsize_t count)
    {
        // Open dataset.
        hid_t dset = h5_check(H5Dopen(file_id, dset_name, H5P_DEFAULT), "opening dataset " + std::string(dset_name));
        hid_t filespace = H5Dget_space(dset);
        hid_t memspace;

        if (count > 0)
        {
            // Select hyperslab in file: [start:start+count]
            h5_check(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL),
                    "selecting hyperslab in " + std::string(dset_name));

            // Local memory layout is simple; continuous 1D array of size `count`.
            memspace = H5Screate_simple(1, &count, NULL);
        }
        else
        {
            // Nothing to write - select empty hyperslab.
            H5Sselect_none(filespace);
            memspace = H5Scopy(filespace);
            H5Sselect_none(memspace);
        }

        // Setup collective IO where all tasks participate.
        hid_t plist_id = create_dxpl();

        // Write data (all tasks must call this even if count=0).
        h5_check(H5Dwrite(dset, get_hdf5_type<T>(), memspace, filespace, plist_id, buffer),
                "writing dataset " + std::string(dset_name));

        // Cleanup!
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Dclose(dset);
    }


    // Particle IDs in the input file must be sorted and equal to 0, 1, ..., n_particles-1,
    // since the particle ID determines the location of each particle in the output files.
    inline void check_particle_uids(
        const std::vector<int>& uid,
        const hsize_t start,
        Master& master)
    {
        int invalid = 0;
        for (size_t n=0; n<uid.size(); ++n)
            if (uid[n] != static_cast<int>(start + n))
            {
                invalid = 1;
                break;
            }

        master.sum(&invalid, 1);

        if (invalid > 0)
            throw std::runtime_error("Particle IDs in input file must be sorted and equal to 0, 1, ..., n_particles-1");
    }


    // Read `particle_id` and properties `names` from file. Each task reads
    // its own block of particles, see `get_io_block()`. Works with and without MPI.
    template <typename TF>
    void read_particles(
        const std::string& filename,
        const std::vector<std::string>& names,
        std::vector<int>& uid,
        const std::vector<std::vector<TF>*>& fields,
        int& n_particles,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        // Open file with (parallel) HDF5.
        hid_t plist_id = create_fapl(master);
        hid_t file_id = h5_check(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_id), "opening " + filename);
        H5Pclose(plist_id);

        // Get total number of particles.
        hid_t dset = h5_check(H5Dopen(file_id, "/particle_id", H5P_DEFAULT), "opening dataset /particle_id in " + filename);
        hid_t dspace = H5Dget_space(dset);
        hsize_t n_total;
        H5Sget_simple_extent_dims(dspace, &n_total, NULL);
        H5Sclose(dspace);
        H5Dclose(dset);

        n_particles = static_cast<int>(n_total);

        hsize_t start, count;
        get_io_block(start, count, n_particles, md.mpiid, md.nprocs);

        // Read each property with hyperslab selection.
        // Use local size without buffer; with MPI, these particles are only
        // temporarely on this task and shipped elsewhwere soon.
        uid.resize(count);
        read_coordinate(file_id, "/particle_id", uid.data(), start, count);

        check_particle_uids(uid, start, master);

        for (size_t i=0; i<names.size(); ++i)
        {
            const std::string path = "/" + names[i];

            if (H5Lexists(file_id, path.c_str(), H5P_DEFAULT) <= 0)
                throw std::runtime_error("Particle property \"" + names[i] + "\" not found in " + filename);

            fields[i]->resize(count);
            read_coordinate(file_id, path.c_str(), fields[i]->data(), start, count);
        }

        // Cleanup!
        h5_check(H5Fclose(file_id), "closing " + filename);
    }


    #ifdef USEMPI
    // Send each particle `n` to MPI task `target_rank[n]`. The received particles are not sorted.
    template<typename TF>
    void alltoallv_particles(
        std::vector<int>& uid_out,
        const std::vector<std::vector<TF>*>& fields_out,
        const std::vector<int>& uid_in,
        const std::vector<std::vector<TF>*>& fields_in,
        const std::vector<int>& target_rank,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        const int np_local = uid_in.size();
        const size_t record_size = get_record_size<TF>(fields_in.size());

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

        // Pack particles into send buffer, ordered by target task.
        std::vector<char> send_buffer(total_send * record_size);
        std::vector<int> current_offset = send_offsets;

        for (int n=0; n<np_local; ++n)
        {
            const int pos = current_offset[target_rank[n]];
            current_offset[target_rank[n]] += 1;

            pack_record(send_buffer.data() + pos * record_size, n, uid_in, fields_in);
        }

        // Exchange particles.
        std::vector<char> recv_buffer(total_recv * record_size);
        MPI_Datatype record_type = create_record_type(record_size);

        MPI_Alltoallv(
            send_buffer.data(),
            send_counts.data(),
            send_offsets.data(),
            record_type,
            recv_buffer.data(),
            recv_counts.data(),
            recv_offsets.data(),
            record_type,
            md.commxy);

        MPI_Type_free(&record_type);

        // Unpack records into output arrays.
        uid_out.resize(total_recv);
        for (std::vector<TF>* field : fields_out)
            field->resize(total_recv);

        for (int n=0; n<total_recv; ++n)
            unpack_record(recv_buffer.data() + n * record_size, n, uid_out, fields_out);
    }


    // Send particles from their "home" task to the task that owns their (x,y) location.
    template<typename TF>
    void distribute_particles(
        std::vector<int>& uid_out,
        const std::vector<std::vector<TF>*>& fields_out,
        const std::vector<int>& uid_in,
        const std::vector<std::vector<TF>*>& fields_in,
        const std::vector<TF>& x_in,
        const std::vector<TF>& y_in,
        const TF xsize,
        const TF ysize,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        const TF xsize_sub = xsize / md.npx;
        const TF ysize_sub = ysize / md.npy;

        // Calculate target `mpiid` for each particle.
        std::vector<int> target_rank(uid_in.size());
        for (size_t n=0; n<uid_in.size(); ++n)
        {
            const int mpicoordx = int(x_in[n] / xsize_sub);
            const int mpicoordy = int(y_in[n] / ysize_sub);

            target_rank[n] = master.get_mpiid(mpicoordx, mpicoordy);
        }

        alltoallv_particles(uid_out, fields_out, uid_in, fields_in, target_rank, master);
    }
    #endif


    // Create new particle file with `particle_id` dimension scale.
    inline hid_t create_particle_file(
        const std::string& filename,
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
        master.broadcast(&file_exists, 1, 0);

        if (file_exists)
        {
            std::string error_msg = "ERROR: Particle output file already exists: " + filename;
            throw std::runtime_error(error_msg);
        }

        // Create file with parallel HDF5.
        hid_t plist_id = create_fapl(master);
        hid_t file_id = h5_check(H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, plist_id), "creating " + filename);
        H5Pclose(plist_id);

        // Create particle_id dimension scale dataset.
        hsize_t particle_dims[1] = {static_cast<hsize_t>(n_particles)};
        hid_t dspace_particle_id = H5Screate_simple(1, particle_dims, particle_dims);
        hid_t dset_particle_id = h5_check(H5Dcreate(
            file_id, "/particle_id", H5T_NATIVE_INT, dspace_particle_id,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), "creating dataset /particle_id");

        // Make particle_id a dimension scale
        h5_check(H5DSset_scale(dset_particle_id, "particle_id"), "setting dimension scale /particle_id");

        H5Dclose(dset_particle_id);
        H5Sclose(dspace_particle_id);

        // Write particle IDs (0, 1, 2, ..., n_particles-1) as coordinate values.
        // Each task only writes its own block of particle IDs.
        hsize_t start, count;
        get_io_block(start, count, n_particles, md.mpiid, md.nprocs);

        std::vector<int> particle_ids(count);
        for (hsize_t i=0; i<count; ++i)
            particle_ids[i] = static_cast<int>(start + i);

        write_coordinate_1d(file_id, "/particle_id", particle_ids.data(), start, count);

        return file_id;
    }


    // Create dataset for particle property `name`, with `long_name` and `units` attributes,
    // and `particle_id` (and optionally `time`) as dimension scales.
    template<typename TF>
    void create_property_dataset(
        hid_t file_id,
        const std::string& name,
        const Particle_property<TF>& prop,
        hid_t dspace,
        hid_t plist_create,
        const bool has_time)
    {
        const std::string path = "/" + name;

        hid_t dset = h5_check(H5Dcreate(
            file_id, path.c_str(), get_hdf5_type<TF>(), dspace,
            H5P_DEFAULT, plist_create, H5P_DEFAULT), "creating dataset " + path);

        h5_check(H5LTset_attribute_string(file_id, path.c_str(), "long_name", prop.long_name.c_str()),
                "setting long_name of " + path);
        h5_check(H5LTset_attribute_string(file_id, path.c_str(), "units", prop.units.c_str()),
                "setting units of " + path);

        // Attach dimension scales.
        hid_t dset_particle_id = h5_check(H5Dopen(file_id, "/particle_id", H5P_DEFAULT), "opening dataset /particle_id");

        if (has_time)
        {
            hid_t dset_time = h5_check(H5Dopen(file_id, "/time", H5P_DEFAULT), "opening dataset /time");
            h5_check(H5DSattach_scale(dset, dset_time, 0), "attaching dimension scale /time to " + path);
            h5_check(H5DSattach_scale(dset, dset_particle_id, 1), "attaching dimension scale /particle_id to " + path);
            H5Dclose(dset_time);
        }
        else
            h5_check(H5DSattach_scale(dset, dset_particle_id, 0), "attaching dimension scale /particle_id to " + path);

        H5Dclose(dset_particle_id);
        H5Dclose(dset);
    }


    template<typename TF>
    void create_particle_dump(
        const std::string& filename,
        hid_t& file_id,
        const int n_particles,
        const int szip_compression,
        const std::vector<std::string>& names,
        Particle_data<TF>& particles,
        Master& master)
    {
        file_id = create_particle_file(filename, n_particles, master);

        // Create unlimited time dimension.
        hsize_t time_dims[1] = {0};
        hsize_t time_maxdims[1] = {H5S_UNLIMITED};
        hid_t dspace_time = H5Screate_simple(1, time_dims, time_maxdims);

        hid_t plist_create = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        h5_check(H5Pset_chunk(plist_create, 1, chunk_dims), "setting chunk size /time");

        hid_t dset_time = h5_check(H5Dcreate(
            file_id, "/time", H5T_NATIVE_DOUBLE, dspace_time,
            H5P_DEFAULT, plist_create, H5P_DEFAULT), "creating dataset /time");

        // Make time a dimension scale
        h5_check(H5DSset_scale(dset_time, "time"), "setting dimension scale /time");

        H5Dclose(dset_time);
        H5Pclose(plist_create);
        H5Sclose(dspace_time);

        // Create particle datasets with unlimited time dimension.
        hsize_t dims[2] = {0, static_cast<hsize_t>(n_particles)};
        hsize_t maxdims[2] = {H5S_UNLIMITED, static_cast<hsize_t>(n_particles)};
        hid_t dspace = H5Screate_simple(2, dims, maxdims);

        hid_t plist_create_2d = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims_2d[2] = {1, static_cast<hsize_t>(n_particles)};
        h5_check(H5Pset_chunk(plist_create_2d, 2, chunk_dims_2d), "setting chunk size particle datasets");

        // Apply SZIP compression if requested.
        if (szip_compression > 0)
        {
            if (H5Zfilter_avail(H5Z_FILTER_SZIP) <= 0)
                throw std::runtime_error("SZIP compression requested, but SZIP filter is not available in HDF5");

            h5_check(H5Pset_szip(plist_create_2d, H5_SZIP_NN_OPTION_MASK, szip_compression), "setting SZIP compression");
        }

        for (const std::string& name : names)
            create_property_dataset(file_id, name, particles.get_property(name), dspace, plist_create_2d, true);

        H5Pclose(plist_create_2d);
        H5Sclose(dspace);

        // Flush to ensure datasets are written to file.
        h5_check(H5Fflush(file_id, H5F_SCOPE_GLOBAL), "flushing file");
    }


    template<typename TF>
    void create_particle_restart(
        const std::string& filename,
        hid_t& file_id,
        const int n_particles,
        const std::vector<std::string>& names,
        Particle_data<TF>& particles,
        Master& master)
    {
        file_id = create_particle_file(filename, n_particles, master);

        // Create particle datasets without time dimension (single time step).
        // Fixed size and no compression, so use default contiguous layout.
        hsize_t dims[1] = {static_cast<hsize_t>(n_particles)};
        hid_t dspace = H5Screate_simple(1, dims, dims);

        for (const std::string& name : names)
            create_property_dataset(file_id, name, particles.get_property(name), dspace, H5P_DEFAULT, false);

        H5Sclose(dspace);

        // Flush to ensure datasets are written to file.
        h5_check(H5Fflush(file_id, H5F_SCOPE_GLOBAL), "flushing file");
    }


    template<typename TF>
    void write_particles_restart(
        hid_t file_id,
        const std::vector<std::string>& names,
        const std::vector<std::vector<TF>*>& fields,
        const int n_particles,
        const int mpiid,
        const int nprocs)
    {
        // Calculate hyperslab parameters for this MPI task.
        hsize_t start, count;
        get_io_block(start, count, n_particles, mpiid, nprocs);

        // Write particle properties (all tasks must call this even if count=0).
        for (size_t i=0; i<names.size(); ++i)
            write_coordinate_1d(file_id, ("/" + names[i]).c_str(), fields[i]->data(), start, count);

        // Flush to disk to ensure data is written.
        h5_check(H5Fflush(file_id, H5F_SCOPE_GLOBAL), "flushing file");
    }


    template<typename TF>
    void write_particles_dump(
        hid_t file_id,
        const std::vector<std::string>& names,
        const std::vector<std::vector<TF>*>& fields,
        const double time,
        const int n_particles,
        const int mpiid,
        const int nprocs)
    {
        // Get current time dimension size_dump
        hid_t dset_time = h5_check(H5Dopen(file_id, "/time", H5P_DEFAULT), "opening dataset /time");
        hid_t dspace_time = H5Dget_space(dset_time);
        hsize_t time_idx;
        H5Sget_simple_extent_dims(dspace_time, &time_idx, NULL);
        H5Sclose(dspace_time);
        H5Dclose(dset_time);

        // Extend datasets to add new time step.
        hsize_t new_dims_time[1] = {time_idx + 1};
        hsize_t new_dims[2] = {time_idx + 1, static_cast<hsize_t>(n_particles)};

        // Extend time dataset.
        dset_time = h5_check(H5Dopen(file_id, "/time", H5P_DEFAULT), "opening dataset /time");
        h5_check(H5Dset_extent(dset_time, new_dims_time), "extending dataset /time");

        // Write time value.
        hid_t dspace_time_write = H5Dget_space(dset_time);
        hsize_t time_start[1] = {time_idx};
        hsize_t time_count[1] = {1};
        h5_check(H5Sselect_hyperslab(dspace_time_write, H5S_SELECT_SET, time_start, NULL, time_count, NULL),
                "selecting hyperslab in /time");
        hid_t memspace_time = H5Screate_simple(1, time_count, NULL);

        hid_t plist_xfer = create_dxpl();

        h5_check(H5Dwrite(dset_time, H5T_NATIVE_DOUBLE, memspace_time, dspace_time_write, plist_xfer, &time),
                "writing dataset /time");

        H5Pclose(plist_xfer);
        H5Sclose(memspace_time);
        H5Sclose(dspace_time_write);
        H5Dclose(dset_time);

        // Calculate hyperslab parameters for this MPI task.
        hsize_t start, count;
        get_io_block(start, count, n_particles, mpiid, nprocs);

        // Extend and write particle properties (all tasks must call this even if count=0).
        for (size_t i=0; i<names.size(); ++i)
        {
            const std::string path = "/" + names[i];

            hid_t dset = h5_check(H5Dopen(file_id, path.c_str(), H5P_DEFAULT), "opening dataset " + path);
            h5_check(H5Dset_extent(dset, new_dims), "extending dataset " + path);
            H5Dclose(dset);

            write_coordinate(file_id, path.c_str(), fields[i]->data(), time_idx, start, count, n_particles);
        }

        // Flush to disk to ensure data is written.
        h5_check(H5Fflush(file_id, H5F_SCOPE_GLOBAL), "flushing file");
    }


    // Gather particles back to their "home" task based on `uid`, sorted by `uid`.
    template<typename TF>
    void gather_particles(
        std::vector<int>& uid_out,
        const std::vector<std::vector<TF>*>& fields_out,
        const std::vector<int>& uid_in,
        const std::vector<std::vector<TF>*>& fields_in,
        const int n_particles,
        Master& master)
    {
        auto& md = master.get_MPI_data();

        const int n_fields = fields_in.size();

        // Unsorted particles on their "home" task.
        std::vector<int> uid_tmp;
        std::vector<std::vector<TF>> fields_tmp(n_fields);
        std::vector<std::vector<TF>*> fields_tmp_ptr(n_fields);
        for (int f=0; f<n_fields; ++f)
            fields_tmp_ptr[f] = &fields_tmp[f];

        #ifdef USEMPI
        // Each task originally has `ceil(n_particles / nprocs)` particles.
        const int np_per_task = get_np_per_task(n_particles, md.nprocs);

        // Calculate which original MPI task each particle belongs to based on uid.
        std::vector<int> target_rank(uid_in.size());
        for (size_t n=0; n<uid_in.size(); ++n)
            target_rank[n] = uid_in[n] / np_per_task;

        alltoallv_particles(uid_tmp, fields_tmp_ptr, uid_in, fields_in, target_rank, master);

        #else
        // Single task: all particles are already local.
        uid_tmp = uid_in;
        for (int f=0; f<n_fields; ++f)
            fields_tmp[f] = *fields_in[f];
        #endif

        const int total_recv = uid_tmp.size();

        // Sort particles by uid to restore original order.
        std::vector<int> order(total_recv);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](const int a, const int b) { return uid_tmp[a] < uid_tmp[b]; });

        uid_out.resize(total_recv);
        for (int n=0; n<total_recv; ++n)
            uid_out[n] = uid_tmp[order[n]];

        for (int f=0; f<n_fields; ++f)
        {
            std::vector<TF>& field_out = *fields_out[f];
            field_out.resize(total_recv);

            for (int n=0; n<total_recv; ++n)
                field_out[n] = fields_tmp[f][order[n]];
        }

        // Each task should now have exactly its own block of particles with uid [start, start+count).
        // If not, particles have been lost or duplicated, and the output would be corrupted.
        hsize_t start, count;
        get_io_block(start, count, n_particles, md.mpiid, md.nprocs);

        int n_invalid = (static_cast<hsize_t>(total_recv) == count) ? 0 : 1;
        for (int n=0; n<total_recv && n_invalid == 0; ++n)
            if (uid_out[n] != static_cast<int>(start + n))
                n_invalid = 1;

        master.sum(&n_invalid, 1);

        if (n_invalid > 0)
            throw std::runtime_error(
                "Lagrangian particles lost or duplicated on " + std::to_string(n_invalid) + " task(s)");
    }
}
#endif
