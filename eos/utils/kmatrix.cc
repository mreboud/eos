/*
 * Copyright (c) 2019 Stephan Kürten
 * Copyright (c) 2019 Danny van Dyk
 * Copyright (c) 2021 Méril Reboud
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <eos/utils/kmatrix.hh>
#include <eos/maths/power-of.hh>

namespace eos
{
    namespace Kmatrix_utils
    {
        double BWfactor(const unsigned & l, const double & q, const double & qR)
        {
            if (qR == 0)
            {
                throw InternalError("Blatt-Weisskopf factors are not defined for a resonance of zero size.");
            }

            if (q < 0)
            {
                throw InternalError("Blatt-Weisskopf factors are not defined for negative momenta.");
            }

            const double z = q / qR;

            switch (l)
            {
                case 0:
                    return 1;
                case 1:
                    return std::sqrt(2 * z / (z + 1));
                case 2:
                    return std::sqrt(13 * power_of<2>(z) / (power_of<2>(z - 3) + 9 * z));
                case 3:
                    return std::sqrt(277 * power_of<3>(z) / (z * power_of<2>(z - 15) + 9 * power_of<2>(2 * z - 5))); // This one is wrong in CBHKSS:1995
                case 4:
                    return std::sqrt(12746 * power_of<4>(z) / (power_of<2>(power_of<2>(z) - 45 * z + 105) + 25 * z * power_of<2>(2 * z - 21)));
                default:
                    throw InternalError("Blatt-Weisskopf factors are not implemented for l > 4.");
            }
        }
    }
}