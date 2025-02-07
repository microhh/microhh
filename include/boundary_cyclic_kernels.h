#pragma once


namespace Boundary_cyclic_kernels
{
    #ifdef USEMPI
    template<typename TF> MPI_Datatype mpi_data_type();
    template<> MPI_Datatype mpi_data_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_data_type<float>() { return MPI_FLOAT; }
    template<> MPI_Datatype mpi_data_type<unsigned int>() { return MPI_UNSIGNED; }


    template<typename TF, bool use_gpu>
    inline void cyclic_kernel(
            TF* __restrict__ const data,
            TF* __restrict__ const buffer_send,
            TF* __restrict__ const buffer_recv,
            const Edge edge,
            const int istart, const int iend, const int jstart, const int jend,
            const int icells, const int jcells, const int kcells,
            const int jj, const int kk,
            const MPI_data& md)
    {
        const int igc = istart;
        const int jgc = jstart;

        if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
        {
            const int jj_buf = igc;
            const int kk_buf = igc*jcells;

            // Pack the buffer.
            #pragma acc parallel loop deviceptr(data, buffer_send) collapse(3) if (use_gpu)
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma GCC ivdep
                    for (int i=0; i<igc; ++i)
                    {
                        const int ijk_buf = i + j*jj_buf + k*kk_buf;
                        const int ijk = iend-igc+i + j*jj + k*kk;
                        buffer_send[ijk_buf] = data[ijk];
                    }

            MPI_Request reqs[2];
            MPI_Isend(buffer_send, igc*jcells*kcells, mpi_data_type<TF>(), md.neast, 1, md.commxy, &reqs[0]);
            MPI_Irecv(buffer_recv, igc*jcells*kcells, mpi_data_type<TF>(), md.nwest, 1, md.commxy, &reqs[1]);
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

            // Unpack the buffer.
            #pragma acc parallel loop deviceptr(data, buffer_recv) collapse(3) if (use_gpu)
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma GCC ivdep
                    for (int i=0; i<igc; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijk_buf = i + j*jj_buf + k*kk_buf;
                        data[ijk] = buffer_recv[ijk_buf];
                    }

            // Pack the buffer.
            #pragma acc parallel loop deviceptr(data, buffer_send) collapse(3) if (use_gpu)
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma GCC ivdep
                    for (int i=0; i<igc; ++i)
                    {
                        const int ijk_buf = i + j*jj_buf + k*kk_buf;
                        const int ijk = i+istart + j*jj + k*kk;
                        buffer_send[ijk_buf] = data[ijk];
                    }

            MPI_Isend(buffer_send, igc*jcells*kcells, mpi_data_type<TF>(), md.nwest, 2, md.commxy, &reqs[0]);
            MPI_Irecv(buffer_recv, igc*jcells*kcells, mpi_data_type<TF>(), md.neast, 2, md.commxy, &reqs[1]);
            MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

            // Unpack the buffer.
            #pragma acc parallel loop deviceptr(data, buffer_recv) collapse(3) if (use_gpu)
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma GCC ivdep
                    for (int i=0; i<igc; ++i)
                    {
                        const int ijk = i+iend + j*jj + k*kk;
                        const int ijk_buf = i + j*jj_buf + k*kk_buf;
                        data[ijk] = buffer_recv[ijk_buf];
                    }
        }

        if (edge == Edge::North_south_edge || edge == Edge::Both_edges)
        {
            // if the run is 3D, apply the BCs
            if ((jend - jstart) > 1)
            {
                const int jj_buf = icells;
                const int kk_buf = icells*jgc;

                // Pack the buffer.
                #pragma acc parallel loop deviceptr(data, buffer_send) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk_buf = i + j*jj_buf + k*kk_buf;
                            const int ijk = i + (jend-jgc+j)*jj + k*kk;
                            buffer_send[ijk_buf] = data[ijk];
                        }

                MPI_Request reqs[2];
                MPI_Isend(buffer_send, icells*jgc*kcells, mpi_data_type<TF>(), md.nnorth, 1, md.commxy, &reqs[0]);
                MPI_Irecv(buffer_recv, icells*jgc*kcells, mpi_data_type<TF>(), md.nsouth, 1, md.commxy, &reqs[1]);
                MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

                // Unpack the buffer.
                #pragma acc parallel loop deviceptr(data, buffer_recv) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk = i + j*jj + k*kk;
                            const int ijk_buf = i + j*jj_buf + k*kk_buf;
                            data[ijk] = buffer_recv[ijk_buf];
                        }

                // Pack the buffer.
                #pragma acc parallel loop deviceptr(data, buffer_send) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk_buf = i + j*jj_buf + k*kk_buf;
                            const int ijk = i + (j+jstart)*jj + k*kk;
                            buffer_send[ijk_buf] = data[ijk];
                        }

                MPI_Isend(buffer_send, icells*jgc*kcells, mpi_data_type<TF>(), md.nsouth, 2, md.commxy, &reqs[0]);
                MPI_Irecv(buffer_recv, icells*jgc*kcells, mpi_data_type<TF>(), md.nnorth, 2, md.commxy, &reqs[1]);
                MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

                // Unpack the buffer.
                #pragma acc parallel loop deviceptr(data, buffer_recv) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk = i + (j+jend)*jj + k*kk;
                            const int ijk_buf = i + j*jj_buf + k*kk_buf;
                            data[ijk] = buffer_recv[ijk_buf];
                        }
            }
            // in case of 2D, fill all the ghost cells with the current value
            else
            {
                #pragma acc parallel loop deviceptr(data) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijkref   = i + jstart*jj   + k*kk;
                            const int ijknorth = i + j*jj        + k*kk;
                            const int ijksouth = i + (jend+j)*jj + k*kk;
                            data[ijknorth] = data[ijkref];
                            data[ijksouth] = data[ijkref];
                        }
            }
        }
    }
    #else
    template<typename TF, bool use_gpu>
    inline void cyclic_kernel(
            TF* __restrict__ const data,
            TF* __restrict__ const buffer_send,
            TF* __restrict__ const buffer_recv,
            const Edge edge,
            const int istart, const int iend, const int jstart, const int jend,
            const int icells, const int jcells, const int kcells,
            const int jj, const int kk,
            const MPI_data& md)
    {
        const int igc = istart;
        const int jgc = jstart;

        if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
        {
            #pragma acc parallel loop deviceptr(data) collapse(3) if (use_gpu)
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma GCC ivdep
                    for (int i=0; i<igc; ++i)
                    {
                        const int ijk0 = i          + j*jj + k*kk;
                        const int ijk1 = iend-igc+i + j*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }

            #pragma acc parallel loop deviceptr(data) collapse(3) if (use_gpu)
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    #pragma GCC ivdep
                    for (int i=0; i<igc; ++i)
                    {
                        const int ijk0 = i+iend   + j*jj + k*kk;
                        const int ijk1 = i+istart + j*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }
        }

        if (edge == Edge::North_south_edge || edge == Edge::Both_edges)
        {
            // if the run is 3D, apply the BCs
            if ((jend - jstart) > 1)
            {
                // second, send and receive the ghost cells in the north-south direction
                #pragma acc parallel loop deviceptr(data) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk0 = i + j           *jj + k*kk;
                            const int ijk1 = i + (jend-jgc+j)*jj + k*kk;
                            data[ijk0] = data[ijk1];
                        }

                #pragma acc parallel loop deviceptr(data) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijk0 = i + (j+jend  )*jj + k*kk;
                            const int ijk1 = i + (j+jstart)*jj + k*kk;
                            data[ijk0] = data[ijk1];
                        }
            }
            // in case of 2D, fill all the ghost cells with the current value
            else
            {
                #pragma acc parallel loop deviceptr(data) collapse(3) if (use_gpu)
                for (int k=0; k<kcells; ++k)
                    for (int j=0; j<jgc; ++j)
                        #pragma GCC ivdep
                        for (int i=0; i<icells; ++i)
                        {
                            const int ijkref   = i + jstart*jj   + k*kk;
                            const int ijknorth = i + j*jj        + k*kk;
                            const int ijksouth = i + (jend+j)*jj + k*kk;
                            data[ijknorth] = data[ijkref];
                            data[ijksouth] = data[ijkref];
                        }
            }
        }
    }
    #endif
}
