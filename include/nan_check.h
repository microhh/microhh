/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#ifndef NAN_CHECK_H
#define NAN_CHECK_H

#include <cmath>
#include <cstdint>
#include <cstring>

namespace Nan_check
{
    // Bit-level finiteness test; immune to fast-math / finite-math-only assumptions.
    inline bool is_finite_bits(const double x)
    {
        std::uint64_t u;
        std::memcpy(&u, &x, sizeof(u));
        return (u & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
    }

    inline bool is_finite_bits(const float x)
    {
        std::uint32_t u;
        std::memcpy(&u, &x, sizeof(u));
        return (u & 0x7f800000u) != 0x7f800000u;
    }

    template<typename TF>
    inline bool is_finite(const TF x)
    {
        #if defined(_CRAYC) || defined(__cray__)
        // Cray `-Ofast` assumes finite math and optimizes `std::isfinite` away: test the exponent bits.
        return is_finite_bits(x);
        #else
        return std::isfinite(x);
        #endif
    }
}
#endif
