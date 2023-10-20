/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef DUMP_H
#define DUMP_H

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;

template<typename TF>
class Dump
{
    public:
        Dump(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Dump();

        void init(double);
        void create();

        unsigned long get_time_limit(unsigned long);
        bool get_switch() { return swdump; }
        std::vector<std::string>& get_dumplist();

        bool do_dump(unsigned long, unsigned long);
        void save_dump(TF*, const std::string&, int);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_io<TF> field3d_io;

        std::vector<std::string> dumplist; // List with all dumps from the ini file.
        bool swdump;                       // Statistics on/off switch
        bool swdoubledump;                 // On/off switch for two consecutive dumps in time
        double sampletime;
        unsigned long isampletime;
};
#endif
