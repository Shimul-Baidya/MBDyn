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

#ifndef ___SOLID_PRESS_H__INCLUDED___
#define ___SOLID_PRESS_H__INCLUDED___

#include "dataman.h"

// 2D elements
class Quadrangle4;
class Quadrangle8;
class Quadrangle9;
class Quadrangle8r;
class Triangle6h;

// 2D collocation rules
class Gauss2x2;
class Gauss2x2Lumped;
class Gauss3x3;
class Gauss3x3Lumped;
class CollocTria6h;

// 2D base class
class SurfaceLoadElem: virtual public Elem, public InitialAssemblyElem {
public:
     SurfaceLoadElem(unsigned uLabel,
                      flag fOut);
     virtual ~SurfaceLoadElem();

     virtual Elem::Type GetElemType() const override;

     virtual void
     SetValue(DataManager *pDM,
              VectorHandler& X, VectorHandler& XP,
              SimulationEntity::Hints *ph) override;

     virtual std::ostream& Restart(std::ostream& out) const override;

     virtual unsigned int iGetInitialNumDof() const override;

     virtual bool bIsDeformable() const override;

     virtual void Output(OutputHandler& OH) const override;

protected:
     Vec3 Ftot;
};

template <typename ElementType, typename CollocationType>
SurfaceLoadElem*
ReadPressureLoad(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

template <typename ElementType, typename CollocationType>
SurfaceLoadElem*
ReadTractionLoad(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

template<typename ElementType, typename CollocationType>
SurfaceLoadElem*
ReadUnilateralInPlaneContact(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
#endif
