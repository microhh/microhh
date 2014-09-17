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

namespace fd
{
  namespace o4
  {
    // 4th order interpolation
    const double ci0  = -1./16.;
    const double ci1  =  9./16.;
    const double ci2  =  9./16.;
    const double ci3  = -1./16.;
                             
    const double bi0  =  5./16.;
    const double bi1  = 15./16.;
    const double bi2  = -5./16.;
    const double bi3  =  1./16.;
                               ;
    const double ti0  =  1./16.;
    const double ti1  = -5./16.;
    const double ti2  = 15./16.;
    const double ti3  =  5./16.;

    // 4th order gradient
    const double cg0  =   1.;
    const double cg1  = -27.;
    const double cg2  =  27.;
    const double cg3  =  -1.;
    const double cgi  =   1./24.;

    const double bg0  = -23.;
    const double bg1  =  21.;
    const double bg2  =   3.;
    const double bg3  =  -1.;

    const double tg0  =   1.;
    const double tg1  =  -3.;
    const double tg2  = -21.;
    const double tg3  =  23.;

    //// 4th order divgrad
    const double cdg0 = -1460./576.;
    const double cdg1 =   783./576.;
    const double cdg2 =   -54./576.;
    const double cdg3 =     1./576.;
  }
}
