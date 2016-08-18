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

#ifndef IMMERSED_BOUNDARY
#define IMMERSED_BOUNDARY

#include <vector>
#include <bitset>

class Model;
class Grid;
class Fields;
class Input;
class Stats;
struct Mask;

struct Neighbour
{
    int i;
    int j;
    int k;
    double distance;
};

struct Ghost_cell
{
    int i;  ///< x-index of ghost cell
    int j;  ///< y-index of ghost cell
    int k;  ///< z-index of ghost cell

    double xB;  ///< Nearest x location on immersed boundary
    double yB;  ///< Nearest y location on immersed boundary
    double zB;  ///< Nearest z location on immersed boundary

    double xI;  ///< x-Location of image point of ghost cell
    double yI;  ///< x-Location of image point of ghost cell
    double zI;  ///< x-Location of image point of ghost cell

    std::vector< std::vector<double> > B; ///< Inversed matrix with the locations of all interpolation points
    std::vector<Neighbour> neighbours;    ///< Neighbouring fluid points used in interpolation

    std::vector<double> c_idw; ///< Weights for inverse distance weighted interpolation
    double c_idw_sum;          ///< Sum inverse distance weights
};

class Immersed_boundary
{
    public:
        enum IB_type{None_type, Sine_type, Gaus_type, Agnesi_type, Block_type, User_type};
        enum Interpol_type{Linear_type, Distance_type};

        Immersed_boundary(Model*, Input*); ///< Constructor of the class.
        ~Immersed_boundary();              ///< Destructor of the class.

        void init();
        void create();
        void exec();            ///< Set the immersed boundary ghost cells
        void exec_stats(Mask*); ///< Execute statistics of immersed boundaries
        void get_mask(Field3d*, Field3d*);

    private:
        Model*  model;  ///< Pointer to model class.
        Fields* fields; ///< Pointer to fields class.
        Grid*   grid;   ///< Pointer to grid class.
        Stats*  stats;  ///< Pointer to grid class.
        
        ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_u;  
        std::vector<Ghost_cell> ghost_cells_v;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_w;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_s;  ///< Vector holding info on all the ghost cells within the boundary

        template<IB_type, int> 
        void calc_mask(double*, double*, double*, int*, int*, int*, const double*, const double*, const double*, const double*); 
        template<IB_type, int> 
        void find_ghost_cells(std::vector<Ghost_cell>*, const double*, const double*, const double*); ///< Function which determines the ghost cells
        void read_ghost_cells(std::vector<Ghost_cell>*, std::string, const double*, const double*, const double*); ///< Function to read user input IB
        template<IB_type, int> 
        double boundary_function(double, double); ///< Function describing boundary
        template<IB_type, int> 
        bool is_ghost_cell(const double*, const double*, const double*, const int, const int, const int); ///< Function which checks if a cell is a ghost cell
        template<IB_type, int> 
        void find_nearest_location_wall(double&, double&, double&, double&, 
                                        const double, const double, const double, const int, const int, const int); ///< Function which checks if a cell is a ghost cell
        template<IB_type, int> 
        void find_interpolation_points(Ghost_cell&, const double*, const double*, const double*, const int, const int, const int, const int); ///< Function which checks if a cell is a ghost cell

//        void define_distance_matrix(Ghost_cell&, const double*, const double*, const double*);

        // General settings IB
        std::string sw_ib; ///< Namelist IB switch
        IB_type ib_type;   ///< Internal IB switch

        // Interpolation settings
        std::string sw_interpol;     ///< Namelist interpolation switch
        Interpol_type interpol_type; ///< Internal interpolation switch
        int idw_points;              ///< Number of points used for inverse distance weighted interpolation

        double amplitude; ///< Height of IB object (Gaussian, sine or blocks)
        double z_offset;  ///< Vertical offset of IB objects
        int xy_dims;      ///< Hill dimension (1=x, 2=xy)

        // Sine type of boundary
        double wavelength_x; ///< Wave length sine in x-direction
        double wavelength_y; ///< Wave length sine in y-direction

        // Gaussian hill
        double x0_hill;      ///< Center of hill in x-direction 
        double y0_hill;      ///< Center of hill in y-direction 
        double sigma_x_hill; ///< Std.dev hill width in x-direction
        double sigma_y_hill; ///< Std.dev hill width in y-direction
};
#endif
