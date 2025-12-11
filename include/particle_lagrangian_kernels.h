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

#ifndef PARTICLE_LAGRANGIAN_KERNELS_H
#define PARTICLE_LAGRANGIAN_KERNELS_H


namespace Particle_lagrangian_kernels
{
    template<typename TF>
    struct Particle_comm
    {
        int uid;
        TF x, y, z;     // Location.
        TF xt, yt, zt;  // Location tendency.
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
    MPI_Datatype create_particle_type_comm()
    {
        // Create MPI type specific for particle communication, which typically
        // uses more elements per particle than the parallel I/O.
        MPI_Datatype particle_type;

        // 1=uid, 6=number of particle properties (currently x,y,z + x,y,z tendencies).
        int blocklengths[2] = {1, 6};
        MPI_Aint displacements[2];
        MPI_Datatype types[2] = {MPI_INT, get_mpi_type<TF>()};

        Particle_comm<TF> particle;
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


    template<typename T>
    int adaptive_resize(
            std::vector<T>& v,
            const int new_size,
            const double reserve_ratio)
    {
        const int size     = v.size();      // Used elements.
        const int capacity = v.capacity();  // Reserved size.

        // Keep track of increases and decreases in size.
        int status = 0;     // 0 = no change, 1 = increase, 2 = decrease.

        // Don't enlarge and shrink using the same ratio.
        const double shrink_ratio = 2. * reserve_ratio - 1.;

        if (new_size * shrink_ratio < capacity)
        {
            // Too large! Decrease capacity.
            status = 2;
            const int new_capacity = int(new_size * reserve_ratio);
            std::vector<T>(v.begin(), v.begin() + new_size).swap(v);
            v.reserve(new_capacity);

            //std::cout << "Too large! Size=" << size << " -> new size=" << new_size << ", new capacity=" << new_capacity << std::endl;
        }
        else if (new_size > capacity)
        {
            // Too small! Increase capacity.
            status = 1;
            const int new_capacity = int(new_size * reserve_ratio);
            v.reserve(new_capacity);

            //std::cout << "Too small! Size=" << size << " -> new size=" << new_size << ", new capacity=" << new_capacity << std::endl;
        }

        v.resize(new_size);

        return status;
    }


    template<typename TF>
    void compact_vectors(
        const std::vector<int>& leaving,
        std::vector<int>& uid,
        std::vector<TF>& xp,
        std::vector<TF>& yp,
        std::vector<TF>& zp,
        std::vector<TF>& xpt,
        std::vector<TF>& ypt,
        std::vector<TF>& zpt)
    {
        // Create lookup table to check if an index is leaving.
        // Prevents using a find operation, which is more expensive.
        std::vector<char> is_leaving(xp.size(), 0);
        for (int idx : leaving)
            is_leaving[idx] = 1;

        int last_idx = xp.size()-1;

        // Compact the vectors.
        for (int pos : leaving)
        {
            // If idx = [0,1,2,3,4,5] and `1` and `4` leave, the compacted vector is:
            //    idx = [0,5,2,3]
            // So if the leaving index is behind `last_idx`, there is no need to move it.
            if (pos > last_idx)
                continue;

            // Find a valid index at the end of the vector.
            while (last_idx > pos && is_leaving[last_idx])
                last_idx--;

            // Same as above; skip if already in empty block.
            if (pos >= last_idx)
                continue;

            // Swap data from end of vector to empty spot.
            uid[pos] = uid[last_idx];

            xp[pos] = xp[last_idx];
            yp[pos] = yp[last_idx];
            zp[pos] = zp[last_idx];

            xpt[pos] = xpt[last_idx];
            ypt[pos] = ypt[last_idx];
            zpt[pos] = zpt[last_idx];

            last_idx--;
        }
    }


    template<typename TF>
    void calc_interpolation_factors_h(
        int* const restrict index,
        TF* const restrict factor,
        const TF* const restrict xp,
        const TF x0,
        const TF dxi,
        const int n_particles)
    {
        // For each particle, find index left of value, and calculate interpolation factor.
        // Equidistant grid in the horizontal, so index can be found directly.

        for (int i=0; i<n_particles; ++i)
        {
            const TF fi = (xp[i] - x0) * dxi;
            index[i] = static_cast<int>(fi);

            factor[i] = fi - index[i];
        }
    }


    template<typename TF>
    void calc_interpolation_factors_v(
        int* const restrict index,
        TF* const restrict factor,
        const TF* const restrict zp,
        const TF* const restrict z,
        const TF* const restrict dzi,
        const bool is_half_level,
        const int n_particles,
        const int kcells)
    {
        // For each particle, find index left of value, and calculate interpolation factor.
        // Non-equidistant grid in the vertical, so requires search using `std::upper_bound`.

        // This is slightly annoying...
        // For half levels, the spacing from zh[k] to zh[k+1] = dz[k]
        // For full levels, the spacing from z[k] to z[k+1] = dzh[k+1]
        const int dk = is_half_level ? 0 : 1;

        for (int i=0; i<n_particles; ++i)
        {
            // Use `upper_bound`; our `zh[0]` and `zh[1]` are both zero!
            const TF* it = std::upper_bound(z, z+kcells, zp[i]);
            const int k0 = static_cast<int>(it - z) - 1;

            index[i] = k0;
            factor[i] = (zp[i] - z[k0]) * dzi[k0+dk];
        }
    }


    template<typename TF>
    void diagnose_velocity(
        TF* const restrict vel_p,
        const TF* const restrict vel_3d,
        const int* const restrict il,
        const int* const restrict jl,
        const int* const restrict kl,
        const TF* const restrict fx,
        const TF* const restrict fy,
        const TF* const restrict fz,
        const int n_particles,
        const int jstride,
        const int kstride)
    {
        // Diagnose particle velocity by tri-linear interpolation of Eulerian velocity field to particle location.
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int n=0; n<n_particles; ++n)
        {
            const int ijk = il[n] + jl[n]*jstride + kl[n]*kstride;

            const TF fx1 = fx[n];
            const TF fy1 = fy[n];
            const TF fz1 = fz[n];

            const TF fx0 = TF(1) - fx1;
            const TF fy0 = TF(1) - fy1;
            const TF fz0 = TF(1) - fz1;

            vel_p[n] =
                fx0 * fy0 * fz0 * vel_3d[ijk               ] +
                fx1 * fy0 * fz0 * vel_3d[ijk + ii          ] +
                fx0 * fy1 * fz0 * vel_3d[ijk + jj          ] +
                fx0 * fy0 * fz1 * vel_3d[ijk + kk          ] +
                fx1 * fy1 * fz0 * vel_3d[ijk + ii + jj     ] +
                fx1 * fy0 * fz1 * vel_3d[ijk + ii + kk     ] +
                fx0 * fy1 * fz1 * vel_3d[ijk + jj + kk     ] +
                fx1 * fy1 * fz1 * vel_3d[ijk + ii + jj + kk];
        }
    }


    template<typename TF>
    void add_tendency(
        TF* const restrict tend_p,
        const TF* const restrict vel_p,
        const int n_particles)
    {
        // Add velocity to tendency.
        for (int n=0; n<n_particles; ++n)
            tend_p[n] += vel_p[n];
    }


    template<typename TF>
    void particle_exchange_serial(
        std::vector<TF>& xp,
        std::vector<TF>& yp,
        const TF xsize,
        const TF ysize)
    {
        // Periodic BCs without MPI.
        const TF xsize_inv = TF(1) / xsize;
        const TF ysize_inv = TF(1) / ysize;

        for (int n=0; n<xp.size(); ++n)
        {
            xp[n] -= std::floor(xp[n] * xsize_inv) * xsize;
            yp[n] -= std::floor(yp[n] * ysize_inv) * ysize;
        }
    }


    enum class Neighbor
    {
        SW = 0,
        S  = 1,
        SE = 2,
        E  = 3,
        NE = 4,
        N  = 5,
        NW = 6,
        W  = 7
    };


    template<typename TF>
    void particle_exchange_parallel(
        std::vector<int>& uid,
        std::vector<TF>& xp,
        std::vector<TF>& yp,
        std::vector<TF>& zp,
        std::vector<TF>& xpt,
        std::vector<TF>& ypt,
        std::vector<TF>& zpt,
        const TF xsize,
        const TF ysize,
        const TF reserve_ratio,
        int& mem_inc,
        int& mem_dec,
        Master& master)
    {
        // Neighbour-neighbour + periodic boundary exchange with MPI.

        // ------------------------------------------
        // 1. Setup MPI communication with neighbors.
        // ------------------------------------------
        auto& md = master.get_MPI_data();

        const TF xsize_sub = xsize / md.npx;
        const TF ysize_sub = ysize / md.npy;

        // Particles are bound by `x0 >= xp > x1` and same for `y`.
        const TF x0 =  md.mpicoordx    * xsize_sub;
        const TF x1 = (md.mpicoordx+1) * xsize_sub;

        const TF y0 =  md.mpicoordy    * ysize_sub;
        const TF y1 = (md.mpicoordy+1) * ysize_sub;

        // Define offsets of neighbors.
        const int n_neighbors = 8;
        const int neighbor_coords[8][2] = {
            {-1, -1},  // 0: SW
            { 0, -1},  // 1: S
            { 1, -1},  // 2: SE
            { 1,  0},  // 3: E
            { 1,  1},  // 4: NE
            { 0,  1},  // 5: N
            {-1,  1},  // 6: NW
            {-1,  0},  // 7: W
        };

        // Find `mpiid`'s of neighbors.
        std::vector<int> mpiid_neighbors(n_neighbors);
        for (int i=0; i<n_neighbors; ++i)
        {
            mpiid_neighbors[i] = master.get_mpiid(
                md.mpicoordx + neighbor_coords[i][0],
                md.mpicoordy + neighbor_coords[i][1]);
        }


        // -------------------------------------
        // 2. Setup outgoing particles:
        //    Which leave, and where do they go?
        // -------------------------------------
        std::vector<std::vector<int>> leaving_indices(n_neighbors);
        const int old_size = xp.size();

        for (int n=0; n<old_size; ++n)
        {
            // First check if particle stays; probability is (probably..) much higher.
            if (xp[n] >= x0 && xp[n] < x1 && yp[n] >= y0 && yp[n] < y1)
                continue;

            // Particle leaves; determine to which neighbor.
            Neighbor neighbor_idx;

            if (xp[n] < x0)
            {
                if (yp[n] < y0)
                    neighbor_idx = Neighbor::SW;
                else if (yp[n] >= y1)
                    neighbor_idx = Neighbor::NW;
                else
                    neighbor_idx = Neighbor::W;
            }
            else if (xp[n] >= x1)
            {
                if (yp[n] < y0)
                    neighbor_idx = Neighbor::SE;
                else if (yp[n] >= y1)
                    neighbor_idx = Neighbor::NE;
                else
                    neighbor_idx = Neighbor::E;
            }
            else
            {
                if (yp[n] < y0)
                    neighbor_idx = Neighbor::S;
                else if (yp[n] >= y1)
                    neighbor_idx = Neighbor::N;
            }

            leaving_indices[int(neighbor_idx)].push_back(n);
        }

        // Pack leaving particles into send buffers.
        std::vector<int> send_counts(n_neighbors);
        std::vector<std::vector<Particle_comm<TF>>> particles_to_send(n_neighbors);

        for (int i=0; i<n_neighbors; ++i)
        {
            send_counts[i] = leaving_indices[i].size();
            particles_to_send[i].resize(send_counts[i]);

            for (int j=0; j<send_counts[i]; ++j)
            {
                const int idx = leaving_indices[i][j];
                particles_to_send[i][j].uid = uid[idx];
                particles_to_send[i][j].x   = xp[idx];
                particles_to_send[i][j].y   = yp[idx];
                particles_to_send[i][j].z   = zp[idx];
                particles_to_send[i][j].xt  = xpt[idx];
                particles_to_send[i][j].yt  = ypt[idx];
                particles_to_send[i][j].zt  = zpt[idx];
            }
        }

        // Collect all indexes that are leaving.
        std::vector<int> all_leaving;
        for (int i = 0; i < n_neighbors; ++i)
        {
            all_leaving.insert(
                all_leaving.end(),
                leaving_indices[i].begin(),
                leaving_indices[i].end());
        }

        const int total_leaving = all_leaving.size();
        std::sort(all_leaving.begin(), all_leaving.end());


        // -------------------------------------------------
        // 3. Compact vectors to make them continuous again.
        // -------------------------------------------------
        // TODO: perhaps there is a smarter method. In theory, the amount of incoming and
        //       outgoing particles should +/- balance. Is the compaction then needed? We could
        //       write the incoming particles directly at the positions of the leaving particles.
        //       This is a bit tricky if the balance is uneven, so I’ve kept it simple for now...

        compact_vectors(all_leaving, uid, xp, yp, zp, xpt, ypt, zpt);


        // ----------------------------
        // 4. Setup incoming particles.
        // ----------------------------
        std::vector<int> recv_counts(n_neighbors);
        std::vector<MPI_Request> requests(2 * n_neighbors);

        const int count = 1;
        const int tag = 0;

        for (int i=0; i<n_neighbors; ++i)
        {
            MPI_Isend(
                &send_counts[i],
                count,
                MPI_INT,
                mpiid_neighbors[i],
                tag,
                md.commxy,
                &requests[2*i]);

            MPI_Irecv(
                &recv_counts[i],
                count,
                MPI_INT,
                mpiid_neighbors[i],
                tag,
                md.commxy,
                &requests[2*i+1]);
        }

        MPI_Waitall(2*n_neighbors, requests.data(), MPI_STATUSES_IGNORE);

        int total_incoming = 0;
        for (int i=0; i<n_neighbors; ++i)
            total_incoming += recv_counts[i];

        // Debug..
        //std::cout << md.mpicoordx << ", " << md.mpicoordy << ", leaving=" << total_leaving << ", incoming=" << total_incoming << std::endl;

        // ------------------
        // 5. Resize vectors.
        // ------------------
        // Resize vectors with a `reserve_ratio` margin. This avoids needed
        // to do resize() operations every model iteration, which is costly.

        const int new_size = old_size + total_incoming - total_leaving;

        const int status = adaptive_resize(uid, new_size, reserve_ratio);

        if (status == 1)
            mem_inc += 1;
        else if (status == 2)
            mem_dec += 1;

        adaptive_resize(xp, new_size, reserve_ratio);
        adaptive_resize(yp, new_size, reserve_ratio);
        adaptive_resize(zp, new_size, reserve_ratio);

        adaptive_resize(xpt, new_size, reserve_ratio);
        adaptive_resize(ypt, new_size, reserve_ratio);
        adaptive_resize(zpt, new_size, reserve_ratio);


        // -------------------------------------
        // 6. Exchange particles with neighbors.
        // -------------------------------------
        MPI_Datatype particle_type = create_particle_type_comm<TF>();

        std::vector<Particle_comm<TF>> recv_buffer(total_incoming);
        int recv_offset = 0;

        requests.clear();
        requests.resize(2 * n_neighbors);

        for (int i=0; i<n_neighbors; ++i)
        {
            if (send_counts[i] > 0)
                MPI_Isend(
                    particles_to_send[i].data(),
                    send_counts[i],
                    particle_type,
                    mpiid_neighbors[i],
                    count,
                    md.commxy,
                    &requests[2*i]);
            else
                requests[2*i] = MPI_REQUEST_NULL;

            if (recv_counts[i] > 0)
            {
                MPI_Irecv(
                    &recv_buffer[recv_offset],
                    recv_counts[i],
                    particle_type,
                    mpiid_neighbors[i],
                    count,
                    md.commxy,
                    &requests[2*i + 1]);
                recv_offset += recv_counts[i];
            }
            else
                requests[2*i + 1] = MPI_REQUEST_NULL;
        }

        MPI_Waitall(2*n_neighbors, requests.data(), MPI_STATUSES_IGNORE);
        MPI_Type_free(&particle_type);


        // -----------------------------------------------
        // 7. Write data from receive buffer into vectors.
        // -----------------------------------------------
        const int istart = old_size - total_leaving;
        const int iend = istart + total_incoming;

        for (int i=istart; i<iend; ++i)
        {
            const int idx = i - istart;
            uid[i] = recv_buffer[idx].uid;

            xp[i] = recv_buffer[idx].x;
            yp[i] = recv_buffer[idx].y;
            zp[i] = recv_buffer[idx].z;

            xpt[i] = recv_buffer[idx].xt;
            ypt[i] = recv_buffer[idx].yt;
            zpt[i] = recv_buffer[idx].zt;
        }


        // --------------------------------------
        // 8. Apply periodic boundary conditions.
        // --------------------------------------
        const TF xsize_inv = TF(1) / xsize;
        const TF ysize_inv = TF(1) / ysize;

        for (int n=0; n<xp.size(); ++n)
        {
            xp[n] -= std::floor(xp[n] * xsize_inv) * xsize;
            yp[n] -= std::floor(yp[n] * ysize_inv) * ysize;
        }

    }
}
#endif
