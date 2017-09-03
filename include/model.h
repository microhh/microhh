/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#ifndef MODEL
#define MODEL

#include <string>
#include <memory>

class Master;
class Input;
class Data_block;

template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timeloop;
template<typename> class Boundary;
template<typename> class Pres;

template<typename TF>
class Model
{
    public:
        Model(Master&, int, char**);
        ~Model();

        void init();
        void load_or_save();
        void exec();

    private:
        Master& master;

        std::shared_ptr<Input> input;
        std::shared_ptr<Data_block> profs;

        std::shared_ptr<Grid<TF>> grid;
        std::shared_ptr<Fields<TF>> fields;
        std::shared_ptr<Timeloop<TF>> timeloop;
        std::shared_ptr<Boundary<TF>> boundary;
        std::shared_ptr<Pres<TF>> pres;

        std::string simmode;
        std::string simname;

        void load();
        void save();

        void delete_objects();
        void print_status();
        void calc_stats(std::string);
        void set_time_step();
};
#endif
