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

#include "solidshape.h"

constexpr sp_grad::index_type Quadrangle4::iNumNodes;
constexpr doublereal Quadrangle4::ri[];
constexpr doublereal Quadrangle4::si[];

constexpr sp_grad::index_type Quadrangle8::iNumNodes;
constexpr doublereal Quadrangle8::ri[];
constexpr doublereal Quadrangle8::si[];

constexpr sp_grad::index_type Quadrangle9::iNumNodes;
constexpr doublereal Quadrangle9::ri[];
constexpr doublereal Quadrangle9::si[];

constexpr sp_grad::index_type Quadrangle8r::iNumNodes;
constexpr doublereal Quadrangle8r::ri[];
constexpr doublereal Quadrangle8r::si[];

constexpr sp_grad::index_type Triangle6h::iNumNodes;
constexpr doublereal Triangle6h::ri[];
constexpr doublereal Triangle6h::si[];

constexpr sp_grad::index_type Hexahedron1p::iNumNodes;
constexpr sp_grad::index_type Hexahedron1p::iNumNodesExtrap;
constexpr doublereal Hexahedron1p::ri[];
constexpr doublereal Hexahedron1p::si[];
constexpr doublereal Hexahedron1p::ti[];

constexpr sp_grad::index_type Hexahedron8::iNumNodes;
constexpr sp_grad::index_type Hexahedron8::iNumNodesExtrap;
constexpr doublereal Hexahedron8::ri[];
constexpr doublereal Hexahedron8::si[];
constexpr doublereal Hexahedron8::ti[];

constexpr sp_grad::index_type Hexahedron8p::iNumNodes;
constexpr sp_grad::index_type Hexahedron8p::iNumNodesExtrap;
constexpr doublereal Hexahedron8p::ri[];
constexpr doublereal Hexahedron8p::si[];
constexpr doublereal Hexahedron8p::ti[];

constexpr sp_grad::index_type Hexahedron20::iNumNodes;
constexpr sp_grad::index_type Hexahedron20::iNumNodesExtrap;
constexpr doublereal Hexahedron20::ri[];
constexpr doublereal Hexahedron20::si[];
constexpr doublereal Hexahedron20::ti[];

constexpr sp_grad::index_type Hexahedron27::iNumNodes;
constexpr sp_grad::index_type Hexahedron27::iNumNodesExtrap;
constexpr doublereal Hexahedron27::ri[];
constexpr doublereal Hexahedron27::si[];
constexpr doublereal Hexahedron27::ti[];

constexpr sp_grad::index_type Hexahedron20r::iNumNodes;
constexpr sp_grad::index_type Hexahedron20r::iNumNodesExtrap;
constexpr doublereal Hexahedron20r::ri[];
constexpr doublereal Hexahedron20r::si[];
constexpr doublereal Hexahedron20r::ti[];

constexpr sp_grad::index_type Pentahedron15::iNumNodes;
constexpr sp_grad::index_type Pentahedron15::iNumNodesExtrap;
constexpr doublereal Pentahedron15::ri[];
constexpr doublereal Pentahedron15::si[];
constexpr doublereal Pentahedron15::ti[];

constexpr sp_grad::index_type Tetrahedron10h::iNumNodes;
constexpr sp_grad::index_type Tetrahedron10h::iNumNodesExtrap;
constexpr doublereal Tetrahedron10h::ri[];
constexpr doublereal Tetrahedron10h::si[];
constexpr doublereal Tetrahedron10h::ti[];
