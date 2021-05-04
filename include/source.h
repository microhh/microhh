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

#ifndef SOURCE_H
#define SOURCE_H

class Master;
template<typename> class Grid;
template<typename> class Fields;
class Input;

template<typename TF>
class Source
{
    public:
        Source(Master&, Grid<TF>&, Fields<TF>&, Input&); // Constructor of the buffer class.
        ~Source(); // Destructor of the buffer class.

        void init();
        void create(Input&);

        void exec(); // Add the tendencies created by the damping.

    private:
        Master& master;     // Reference to master class.
        Grid<TF>& grid;     // Reference to grid class.
        Fields<TF>& fields; // Reference to fields class.

        struct Shape
        {
            std::vector<int> range_x;
            std::vector<int> range_y;
            std::vector<int> range_z;
        };

        std::vector<Shape> shape;

        std::string swsource;

        TF x0;
        TF y0;
        TF z0;

        std::vector<std::string> sourcelist;
        std::vector<TF> source_x0;
        std::vector<TF> source_y0;
        std::vector<TF> source_z0;
        std::vector<TF> sigma_x;
        std::vector<TF> sigma_y;
        std::vector<TF> sigma_z;
        std::vector<TF> strength;
        std::vector<TF> line_x;
        std::vector<TF> line_y;
        std::vector<TF> line_z;
        std::vector<TF> norm;

        TF calc_norm(
                const TF* const, const TF, const TF, const TF,
                const TF* const, const TF, const TF, const TF,
                const TF* const, const TF, const TF, const TF,
                std::vector<int>, std::vector<int>, std::vector<int>);

        void calc_source(
                TF* const, const TF* const, const TF, const TF, const TF,
                const TF* const, const TF, const TF, const TF,
                const TF* const, const TF, const TF, const TF,
                std::vector<int>, std::vector<int>, std::vector<int>,
                const TF, TF);

        void add_source(TF* const, const TF* const, std::vector<int>, std::vector<int>, std::vector<int>);
};
#endif

