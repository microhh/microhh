/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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

#ifndef THERMO_DISABLED
#define THERMO_DISABLED

class Master;
class Input;
class Grid;
class Fields;
struct Mask;

class Thermo_disabled : public Thermo
{
  public:
    Thermo_disabled(Model*, Input*);
    virtual ~Thermo_disabled();

    void init();
    void create(Input*);
    void exec();
    void exec_stats(Mask*);

    void exec_cross();
    void exec_dump();

    void get_mask(Field3d*, Field3d*, Mask*);

    // interfacing functions to get buoyancy properties from other classes
    bool check_thermo_field(std::string name);
    void get_thermo_field(Field3d*, Field3d*, std::string name);
    void get_buoyancy_surf(Field3d*);
    void get_buoyancy_fluxbot(Field3d*);
    void get_prog_vars(std::vector<std::string>*);
};
#endif
