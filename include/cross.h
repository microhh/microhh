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

#ifndef CROSS
#define CROSS

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;



template<typename TF>
class Cross
{
    public:
        Cross(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Cross();

        void init(double);
        void create();
        bool get_switch() { return swcross; }

        std::vector<std::string>* get_crosslist();

        unsigned long get_time_limit(unsigned long);
        //int exec(double, unsigned long, int);


        bool do_cross(unsigned long);

        enum Direction {Top_to_bottom, Bottom_to_top};

        int cross_simple(TF*, std::string, int);
        int cross_lngrad(TF*, std::string, int);
        int cross_plane (TF*, std::string, int);
        int cross_path  (TF*, std::string, int);
        int cross_height_threshold(TF*, TF, Direction, std::string, int);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_io<TF> field3d_io;

        bool swcross;
        TF sampletime;
        unsigned long isampletime;

        std::vector<std::string> crosslist; ///< List with all crosses from the ini file.

        std::vector<int> jxz;   ///< Index of nearest full y position of xz input
        std::vector<int> ixz;   ///< Index of nearest full x position of yz input
        std::vector<int> kxy;   ///< Index of nearest full height level of xy input
        std::vector<int> jxzh;  ///< Index of nearest half y position of xz input
        std::vector<int> ixzh;  ///< Index of nearest half x position of yz input
        std::vector<int> kxyh;  ///< Index of nearest half height level of xy input
        std::vector<TF> xz; ///< Y-position [m] xz cross from ini file
        std::vector<TF> yz; ///< X-position [m] yz cross from ini file
        std::vector<TF> xy; ///< Z-position [m] xy cross from ini file

        std::vector<std::string> simple;
        std::vector<std::string> bot;
        std::vector<std::string> fluxbot;
        std::vector<std::string> lngrad;
        std::vector<std::string> path;

        //int check_list(std::vector<std::string> *, FieldMap *, std::string crossname);
        int check_save(int, char *);
};
#endif
