/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef CONSTEXPR_MATH_H
#define CONSTEXPR_MATH_H
#include "ac/f2c.h"
#include <limits>

namespace constexpr_math
{
     namespace implementation
     {
          inline doublereal constexpr sqrtNewtonRaphson(doublereal x, doublereal curr, doublereal prev)
          {
               return curr == prev
                    ? curr
                    : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
          }
     }

/*
 * Constexpr version of the square root
 * Return value:
 *   - For a finite and non-negative value of "x", returns an approximation for the square root of "x"
 *   - Otherwise, returns NaN
 */
     inline doublereal constexpr sqrt(doublereal x)
     {
          return x >= 0 && x < std::numeric_limits<doublereal>::infinity()
               ? implementation::sqrtNewtonRaphson(x, x, 0)
               : std::numeric_limits<doublereal>::quiet_NaN();
     }

     // Perform a few tests with sqrtNewtonRaphson
     static_assert(sqrt(1234.567890123456 * 1234.567890123456) == 1234.567890123456, "unit test for compile time square root failed");
     static_assert(sqrt(87654.32123456789 * 87654.32123456789) == 87654.32123456789, "unit test for compile time square root failed");
}

#endif
