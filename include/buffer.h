/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

#ifndef BUFFER
#define BUFFER

// forward declarations to reduce compilation time
class cmaster;
class cmodel;
class cgrid;
class cfields;

/**
 * Class for the buffer layer in the top of the domain.
 * This class performs the gravity wave damping in the top of the domain to
 * prevent reflection at the top boundary.
 */
class cbuffer
{
  public:
    cbuffer(cmodel *, cinput *); ///< Constructor of the buffer class.
    ~cbuffer();                  ///< Destructor of the buffer class.

    void init();           ///< Initialize the arrays that contain the profiles.
    void create(cinput *); ///< Read the profiles of the forces from the input.
    int exec();            ///< Add the tendencies created by the damping.

    // GPU functions and variables
    int prepareGPU(); ///< Allocate and copy buffer profiles at/to GPU                             


  private:
    cmaster *master; ///< Pointer to master class.
    cmodel  *model;  ///< Pointer to model class.
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.

    double zstart; ///< Height above which the buffer layer is starting.
    double sigma;  ///< Damping frequency.
    double beta;   ///< Exponent for damping increase with height.

    int bufferkstart;  ///< Grid point at cell center at which damping starts.
    int bufferkstarth; ///< Grid point at cell face at which damping starts.

    std::map<std::string, double*> bufferprofs;   ///< Map containing the buffer profiles.

    std::string swbuffer; ///< Switch for buffer.

    int buffer(double * const, const double * const, 
               const double * const, const double * const); ///< Calculate the tendency 

    // GPU functions and variables
    std::map<std::string, double*> bufferprofs_g; ///< Map containing the buffer profiles at GPU.

};
#endif
