/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
#include <stdexcept>

#include "master.h"
#include "grid.h"
#include "boundary_outflow.h"

namespace
{
    template<typename TF, Edge_location location>
    void set_neumann(
            TF* const restrict a,
            const int istart, const int iend, const int igc,
            const int jstart, const int jend, const int jgc,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii = 1;
        int ijk;
        int ijk_gc;

        // Set the ghost cells using extrapolation.
        if (location == Edge_location::West || location == Edge_location::East)
        {
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    for (int i=0; i<igc; ++i)
                    {
                        if (location == Edge_location::West)
                        {
                            ijk    = (istart+i  ) + j*icells + k*ijcells;
                            ijk_gc = (istart-1-i) + j*icells + k*ijcells;
                        }
                        else if (location == Edge_location::East)
                        {
                            ijk    = (iend-1-i) + j*icells + k*ijcells;
                            ijk_gc = (iend+i  ) + j*icells + k*ijcells;
                        }

                        a[ijk_gc] = a[ijk];
                    }
        }
        else if (location == Edge_location::North || location == Edge_location::South)
        {
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jgc; ++j)
                    for (int i=0; i<icells; ++i)
                    {
                        if (location == Edge_location::South)
                        {
                            ijk    = i + (jstart+j  )*icells + k*ijcells;
                            ijk_gc = i + (jstart-1-j)*icells + k*ijcells;
                        }
                        else if (location == Edge_location::North)
                        {
                            ijk    = i + (jend-1-j)*icells + k*ijcells;
                            ijk_gc = i + (jend+j  )*icells + k*ijcells;
                        }

                        a[ijk_gc] = a[ijk];
                    }
        }
    }

    template<typename TF, Edge_location location, Flow_direction direction>
    void compute_inoutflow_2nd(
            TF* const restrict a, const TF* const restrict inflow_prof,
            const int istart, const int iend, const int igc,
            const int jstart, const int jend, const int jgc,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii = 1;
        int ijk;
        int ijk_gc;
        int ijk_d;

        // Set the ghost cells using extrapolation.
        if (location == Edge_location::West || location == Edge_location::East)
        {
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jcells; ++j)
                    for (int i=0; i<igc; ++i)
                    {
                        if (location == Edge_location::West)
                        {
                            ijk_d  = (istart    ) + j*icells + k*ijcells;
                            ijk    = (istart+i  ) + j*icells + k*ijcells;
                            ijk_gc = (istart-1-i) + j*icells + k*ijcells;
                        }
                        else if (location == Edge_location::East)
                        {
                            ijk_d  = (iend-1  ) + j*icells + k*ijcells;
                            ijk    = (iend-1-i) + j*icells + k*ijcells;
                            ijk_gc = (iend+i  ) + j*icells + k*ijcells;
                        }

                        if (direction == Flow_direction::Inflow)
                            a[ijk_gc] = a[ijk_d] - (i+1)*TF(2)*(a[ijk_d] - inflow_prof[k]);
                        else
                            a[ijk_gc] = a[ijk];
                    }
        }
        else if (location == Edge_location::North || location == Edge_location::South)
        {
            for (int k=0; k<kcells; ++k)
                for (int j=0; j<jgc; ++j)
                    for (int i=0; i<icells; ++i)
                    {
                        if (location == Edge_location::South)
                        {
                            ijk_d  = i + (jstart    )*icells + k*ijcells;
                            ijk    = i + (jstart+j  )*icells + k*ijcells;
                            ijk_gc = i + (jstart-1-j)*icells + k*ijcells;
                        }
                        else if (location == Edge_location::North)
                        {
                            ijk_d  = i + (jend-1  )*icells + k*ijcells;
                            ijk    = i + (jend-1-j)*icells + k*ijcells;
                            ijk_gc = i + (jend+j  )*icells + k*ijcells;
                        }

                        if (direction == Flow_direction::Inflow)
                            a[ijk_gc] = a[ijk_d] - (j+1)*TF(2)*(a[ijk_d] - inflow_prof[k]);
                        else
                            a[ijk_gc] = a[ijk];
                    }
        }
    }

    template<typename TF>
    void compute_outflow_4th(
            TF* const restrict a,
            const int iend,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        // Set the ghost cells using extrapolation.
        for (int k=0; k<kcells; ++k)
            for (int j=0; j<jcells; ++j)
            {
                const int ijk = (iend-1) + j*icells + k*ijcells;
                a[ijk+ii1] = TF(2.)*a[ijk] - TF( 3./2.)*a[ijk-ii1] + TF(1./2.)*a[ijk-ii2];
                a[ijk+ii2] = TF(3.)*a[ijk] - TF( 7./2.)*a[ijk-ii1] + TF(3./2.)*a[ijk-ii2];
                a[ijk+ii3] = TF(5.)*a[ijk] - TF(15./2.)*a[ijk-ii1] + TF(7./2.)*a[ijk-ii2];
            }
    }

    template<typename TF>
    void compute_inflow_4th(
            TF* const restrict a, const TF value,
            const int istart,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        // Set the ghost cells using extrapolation.
        for (int k=0; k<kcells; ++k)
            for (int j=0; j<jcells; ++j)
            {
                const int ijk = istart + j*icells + k*ijcells;
                a[ijk-ii1] = value + TF( 9./8.)*a[ijk] - TF( 14./8.)*a[ijk+ii1] + TF( 5./8.)*a[ijk+ii2];
                a[ijk-ii2] = value + TF(33./8.)*a[ijk] - TF( 54./8.)*a[ijk+ii1] + TF(21./8.)*a[ijk+ii2];
                a[ijk-ii3] = value + TF(65./8.)*a[ijk] - TF(110./8.)*a[ijk+ii1] + TF(45./8.)*a[ijk+ii2];
            }
    }
}

template<typename TF>
Boundary_outflow<TF>::Boundary_outflow(Master& masterin, Grid<TF>& gridin, Input& inputin) :
    master(masterin),
    grid(gridin)
{
    std::vector<std::string> outflow_list = inputin.get_list<std::string>(
            "boundary", "scalar_outflow", "", std::vector<std::string>());

    if (outflow_list.size() > 0)
    {
        // Process the in/outflow locations:
        auto process_lbc = [&](std::string edge_name, Edge_location edge)
        {
            std::string direction = inputin.get_item<std::string>("boundary", "flow_direction", edge_name);

            if (direction == "inflow")
                flow_direction[edge] = Flow_direction::Inflow;
            else if (direction == "outflow")
                flow_direction[edge] = Flow_direction::Outflow;
            else
            {
                std::string error = "Inflow direction \"" + direction + "\" is invalid\n";
                throw std::runtime_error(error);
            }
        };

        process_lbc("north", Edge_location::North);
        process_lbc("east",  Edge_location::East);
        process_lbc("south", Edge_location::South);
        process_lbc("west",  Edge_location::West);
    }
}

template<typename TF>
Boundary_outflow<TF>::~Boundary_outflow()
{
}

#ifndef USECUDA
template<typename TF>
void Boundary_outflow<TF>::exec(
        TF* const restrict data,
        const TF* const restrict inflow_prof)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        // Dirichlet BCs on west boundary, Neumann on east boundary,
        // cyclic BCs on north-south boundaries
        if (md.mpicoordx == 0)
            compute_inflow_4th(
                    data, TF(0.),
                    gd.istart,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);

        if (md.mpicoordx == md.npx-1)
            compute_outflow_4th(
                    data,
                    gd.iend,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);
    }
    else if (grid.get_spatial_order() == Grid_order::Second)
    {
        if (md.mpicoordx == 0)
        {
            const Edge_location edge = Edge_location::West;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd<TF, edge, Flow_direction::Inflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd<TF, edge, Flow_direction::Outflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }

        if (md.mpicoordx == md.npx-1)
        {
            const Edge_location edge = Edge_location::East;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd<TF, edge, Flow_direction::Inflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd<TF, edge, Flow_direction::Outflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }

        if (md.mpicoordy == 0)
        {
            const Edge_location edge = Edge_location::South;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd<TF, edge, Flow_direction::Inflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd<TF, edge, Flow_direction::Outflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }

        if (md.mpicoordy == md.npy-1)
        {
            const Edge_location edge = Edge_location::North;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd<TF, edge, Flow_direction::Inflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd<TF, edge, Flow_direction::Outflow>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }
    }
}
#endif

template class Boundary_outflow<double>;
template class Boundary_outflow<float>;
