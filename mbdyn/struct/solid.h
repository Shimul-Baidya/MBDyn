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
     static constexpr sp_grad::index_type iNumNodes = 8;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);
};

class Gauss2 {
private:
     static constexpr sp_grad::index_type iGaussOrder = 2;
public:
     static constexpr sp_grad::index_type iNumEvalPoints = std::pow(iGaussOrder, 3);

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);

private:
     static constexpr sp_grad::index_type ridx[iNumEvalPoints] = {0, 1, 1, 0, 0, 1, 1, 0};
     static constexpr sp_grad::index_type sidx[iNumEvalPoints] = {0, 0, 1, 1, 0, 0, 1, 1};
     static constexpr sp_grad::index_type tidx[iNumEvalPoints] = {0, 0, 0, 0, 1, 1, 1, 1};
     static constexpr doublereal ri[iGaussOrder] = {0.577350269189626, -0.577350269189626};
     static constexpr doublereal alphai[iGaussOrder] = {1.0, 1.0};
};

template <typename ElementType, typename CollocationType>
class SolidElemStatic: virtual public Elem, public ElemGravityOwner, protected ElementType, protected CollocationType {
public:
     static constexpr sp_grad::index_type iNumNodes = ElementType::iNumNodes;
     static constexpr sp_grad::index_type iNumEvalPoints = CollocationType::iNumEvalPoints;
     static constexpr sp_grad::index_type iNumDof = iNumNodes * 3;

     SolidElemStatic(unsigned uLabel,
                     const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                     std::array<std::unique_ptr<ConstitutiveLaw6D>, iNumEvalPoints>&& rgCSL,
                     flag fOut);
     virtual ~SolidElemStatic();

     virtual Elem::Type GetElemType(void) const override;

     virtual void Output(OutputHandler& OH) const override;
     
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

     virtual void
     AssJac(VectorHandler& JacY,
            const VectorHandler& Y,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr,
            VariableSubMatrixHandler& WorkMat) override;

protected:
     template <typename T>
     inline void
     GetDeformations(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                     doublereal dCoef,
                     sp_grad::SpFunctionCall func) const;
     
     inline void
     ComputeStressStrain(sp_grad::SpMatrixA<doublereal, 6, iNumEvalPoints>& epsilon,
                         sp_grad::SpMatrixA<doublereal, 6, iNumEvalPoints>& tau) const;

     template <typename T>
     void
     AssStiffnessVec(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                     sp_grad::SpColVector<T, iNumDof>& R,
                     doublereal dCoef,
                     sp_grad::SpFunctionCall func);
     
     struct CollocData {
          static constexpr sp_grad::index_type iNumNodes = SolidElemStatic::iNumNodes;
          static constexpr sp_grad::index_type iNumDof = SolidElemStatic::iNumDof;
          
          inline void
          Init(sp_grad::index_type iColloc,
               const sp_grad::SpMatrixA<doublereal, 3, iNumNodes>& x0,
               std::unique_ptr<ConstitutiveLaw6D>&& pCSL);

          template <typename T>
          inline void
          StrainMatrix1(const sp_grad::SpMatrix<T, 3, 3>& dF,
                        sp_grad::SpMatrix<T, 6, iNumDof>& BL1) const;

          template <typename T>
          inline sp_grad::SpColVector<T, 6>
          ComputeStress(const sp_grad::SpMatrix<T, 3, 3>& G,
                        const sp_grad::SpMatrix<T, 3, 3>& dF);

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& G_tmp,
                             const sp_grad::SpColVector<doublereal, 6>& sigma_tmp,
                             const sp_grad::SpMatrix<doublereal, 3, 3>& dF_tmp);

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& G,
                             const sp_grad::SpColVector<sp_grad::SpGradient, 6>& sigma,
                             const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& dF) {
          }

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& G,
                             const sp_grad::SpColVector<sp_grad::GpGradProd, 6>& sigma,
                             const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& dF) {
          }          
          
          std::unique_ptr<ConstitutiveLaw6D> pConstLaw;
          sp_grad::SpMatrixA<doublereal, iNumNodes, 3> h0d;
          sp_grad::SpMatrixA<doublereal, 3, 3> J, invJ;
          doublereal detJ;
          sp_grad::SpMatrixA<doublereal, 6, iNumDof> BL0;
          sp_grad::SpMatrixA<doublereal, 3, 3> G;
          sp_grad::SpMatrixA<doublereal, 3, 3> F;
          sp_grad::SpColVectorA<doublereal, 6> sigma;
     };

     sp_grad::SpMatrixA<doublereal, 3, iNumNodes> x0;
     const std::array<const StructDispNodeAd*, iNumNodes> rgNodes;
     std::array<CollocData, iNumEvalPoints> rgCollocData;
};

template <typename SolidElemType>
Elem*
ReadSolid(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif
