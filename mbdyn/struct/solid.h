/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef ___SOLID_H__INCLUDED___
#define ___SOLID_H__INCLUDED___

#include <array>
#include <memory>

#include "strnodead.h"
#include "sp_matvecass.h"
#include "elem.h"
#include "constltp.h"
#include "gravity.h"

class Hexahedron8 {
protected:
     static constexpr sp_grad::index_type NumberOfNodes = 8;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, NumberOfNodes, 3>& h0d1);
};

class Gauss2 {
public:
     static inline constexpr sp_grad::index_type
     iGetNumEvalPoints();

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);

private:
     static constexpr sp_grad::index_type NumberOfEvalPoints = 2;
     static constexpr sp_grad::index_type NumberOfCollocationPoints = NumberOfEvalPoints * NumberOfEvalPoints * NumberOfEvalPoints;
     static constexpr sp_grad::index_type ridx[NumberOfCollocationPoints] = {0, 0, 0, 0, 1, 1, 1, 1};
     static constexpr sp_grad::index_type sidx[NumberOfCollocationPoints] = {0, 0, 1, 1, 0, 0, 1, 1};
     static constexpr sp_grad::index_type tidx[NumberOfCollocationPoints] = {0, 1, 0, 1, 0, 1, 0, 1};
     static constexpr doublereal ri[NumberOfEvalPoints] = {0.577350269189626, -0.577350269189626};
     static constexpr doublereal alphai[NumberOfEvalPoints] = {1.0, 1.0};
};

template <typename ElementType, typename CollocationType>
class SolidElemStatic: virtual public Elem, public ElemGravityOwner, protected ElementType, protected CollocationType {
public:
     static constexpr sp_grad::index_type iGetNumNodes() { return ElementType::NumberOfNodes; }
     static constexpr sp_grad::index_type iGetNumEvalPoints() { return CollocationType::iGetNumEvalPoints(); }

     SolidElemStatic(unsigned uLabel,
                     const std::array<const StructDispNodeAd*, ElementType::NumberOfNodes>& rgNodes,
                     std::array<std::unique_ptr<ConstitutiveLaw6D>, CollocationType::iGetNumEvalPoints()>&& rgCSL,
                     flag fOut);
     virtual ~SolidElemStatic();

     virtual Elem::Type GetElemType(void) const override;

     virtual void
     SetValue(DataManager *pDM,
              VectorHandler& X, VectorHandler& XP,
              SimulationEntity::Hints *ph) override;

     virtual std::ostream& Restart(std::ostream& out) const override;

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            enum sp_grad::SpFunctionCall func);

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     // Used by Newton Krylov solvers only
     // JacY += Jac * Y
     virtual void
     AssJac(VectorHandler& JacY,
            const VectorHandler& Y,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr,
            VariableSubMatrixHandler& WorkMat) override;

protected:
     inline void
     Jacobian(const sp_grad::SpMatrix<doublereal, ElementType::NumberOfNodes, 3>& hd,
              sp_grad::SpMatrix<doublereal, 3, 3>& J);

     sp_grad::SpMatrixA<doublereal, 3, ElementType::NumberOfNodes> x0;
     const std::array<const StructDispNodeAd*, ElementType::NumberOfNodes> rgNodes;
     std::array<std::unique_ptr<ConstitutiveLaw6D>, CollocationType::iGetNumEvalPoints()> rgConstLaw;
};

template <typename SolidElemType>
Elem*
ReadSolid(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif
