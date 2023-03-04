/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2023) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "solidinteg.h"

constexpr sp_grad::index_type Gauss2::iGaussOrder;
constexpr doublereal Gauss2::ri[];
constexpr doublereal Gauss2::alphai[];
constexpr doublereal Gauss2::ri_lumped[];
constexpr doublereal Gauss2::alphai_lumped[];

constexpr sp_grad::index_type Gauss2x2::iNumEvalPoints;
constexpr sp_grad::index_type Gauss2x2::ridx[];
constexpr sp_grad::index_type Gauss2x2::sidx[];

constexpr sp_grad::index_type Gauss3x3::iNumEvalPoints;
constexpr sp_grad::index_type Gauss3x3::ridx[];
constexpr sp_grad::index_type Gauss3x3::sidx[];

constexpr sp_grad::index_type CollocTria6h::iNumEvalPoints;
constexpr doublereal CollocTria6h::A;
constexpr doublereal CollocTria6h::B;
constexpr doublereal CollocTria6h::P1;
constexpr doublereal CollocTria6h::P2;
constexpr doublereal CollocTria6h::zeta[7];
constexpr doublereal CollocTria6h::eta[7];
constexpr doublereal CollocTria6h::w[7];

constexpr sp_grad::index_type Gauss2x2x2::iNumEvalPointsStiffness;
constexpr sp_grad::index_type Gauss2x2x2::iNumEvalPointsMass;
constexpr sp_grad::index_type Gauss2x2x2::iNumEvalPointsMassLumped;
constexpr sp_grad::index_type Gauss2x2x2::ridx[];
constexpr sp_grad::index_type Gauss2x2x2::sidx[];
constexpr sp_grad::index_type Gauss2x2x2::tidx[];

constexpr sp_grad::index_type Gauss3::iGaussOrder;
constexpr sp_grad::index_type Gauss3x3x3::iNumEvalPointsStiffness;
constexpr sp_grad::index_type Gauss3x3x3::iNumEvalPointsMass;
constexpr sp_grad::index_type Gauss3x3x3::iNumEvalPointsMassLumped;
constexpr sp_grad::index_type Gauss3x3x3::ridx[];
constexpr sp_grad::index_type Gauss3x3x3::sidx[];
constexpr sp_grad::index_type Gauss3x3x3::tidx[];

constexpr doublereal Gauss3::ri[];
constexpr doublereal Gauss3::alphai[];
constexpr doublereal Gauss3::ri_lumped[];
constexpr doublereal Gauss3::alphai_lumped[];

constexpr sp_grad::index_type CollocPenta15::M;
constexpr sp_grad::index_type CollocPenta15::N;
constexpr sp_grad::index_type CollocPenta15::M_lumped;
constexpr sp_grad::index_type CollocPenta15::iNumEvalPointsStiffness;
constexpr sp_grad::index_type CollocPenta15::iNumEvalPointsMass;
constexpr sp_grad::index_type CollocPenta15::iNumEvalPointsMassLumped;

constexpr doublereal CollocPenta15::ri[];
constexpr doublereal CollocPenta15::si[];
constexpr doublereal CollocPenta15::wi[];
constexpr doublereal CollocPenta15::ti[];
constexpr doublereal CollocPenta15::alphai[];

constexpr doublereal CollocPenta15::ri_lumped[];
constexpr doublereal CollocPenta15::si_lumped[];
constexpr doublereal CollocPenta15::wi_lumped[];
constexpr doublereal CollocPenta15::ti_lumped[];
constexpr doublereal CollocPenta15::alphai_lumped[];

constexpr sp_grad::index_type CollocTet10h::iNumEvalPointsStiffness;
constexpr sp_grad::index_type CollocTet10h::iNumEvalPointsMass;
constexpr sp_grad::index_type CollocTet10h::iNumEvalPointsMassLumped;

constexpr doublereal CollocTet10h::r1[];
constexpr doublereal CollocTet10h::s1[];
constexpr doublereal CollocTet10h::t1[];
constexpr doublereal CollocTet10h::w1[];

constexpr doublereal CollocTet10h::r2[];
constexpr doublereal CollocTet10h::s2[];
constexpr doublereal CollocTet10h::t2[];
constexpr doublereal CollocTet10h::w2[];

constexpr doublereal CollocTet10h::r3[];
constexpr doublereal CollocTet10h::s3[];
constexpr doublereal CollocTet10h::t3[];
constexpr doublereal CollocTet10h::w3[];
