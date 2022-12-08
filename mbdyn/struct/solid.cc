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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "solid.h"
#include "dataman.h"

constexpr sp_grad::index_type Hexahedron8::NumberOfNodes;

void
Hexahedron8::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                sp_grad::SpMatrix<doublereal, NumberOfNodes, 3>& h0d1)
{
     h0d1(1,1) = ((r(2)+1)*(r(3)+1))/8.0E+0;
     h0d1(1,2) = ((r(1)+1)*(r(3)+1))/8.0E+0;
     h0d1(1,3) = ((r(1)+1)*(r(2)+1))/8.0E+0;
     h0d1(2,1) = -((r(2)+1)*(r(3)+1))/8.0E+0;
     h0d1(2,2) = ((1-r(1))*(r(3)+1))/8.0E+0;
     h0d1(2,3) = ((1-r(1))*(r(2)+1))/8.0E+0;
     h0d1(3,1) = -((1-r(2))*(r(3)+1))/8.0E+0;
     h0d1(3,2) = -((1-r(1))*(r(3)+1))/8.0E+0;
     h0d1(3,3) = ((1-r(1))*(1-r(2)))/8.0E+0;
     h0d1(4,1) = ((1-r(2))*(r(3)+1))/8.0E+0;
     h0d1(4,2) = -((r(1)+1)*(r(3)+1))/8.0E+0;
     h0d1(4,3) = ((r(1)+1)*(1-r(2)))/8.0E+0;
     h0d1(5,1) = ((r(2)+1)*(1-r(3)))/8.0E+0;
     h0d1(5,2) = ((r(1)+1)*(1-r(3)))/8.0E+0;
     h0d1(5,3) = -((r(1)+1)*(r(2)+1))/8.0E+0;
     h0d1(6,1) = -((r(2)+1)*(1-r(3)))/8.0E+0;
     h0d1(6,2) = ((1-r(1))*(1-r(3)))/8.0E+0;
     h0d1(6,3) = -((1-r(1))*(r(2)+1))/8.0E+0;
     h0d1(7,1) = -((1-r(2))*(1-r(3)))/8.0E+0;
     h0d1(7,2) = -((1-r(1))*(1-r(3)))/8.0E+0;
     h0d1(7,3) = -((1-r(1))*(1-r(2)))/8.0E+0;
     h0d1(8,1) = ((1-r(2))*(1-r(3)))/8.0E+0;
     h0d1(8,2) = -((r(1)+1)*(1-r(3)))/8.0E+0;
     h0d1(8,3) = -((r(1)+1)*(1-r(2)))/8.0E+0;
}

constexpr sp_grad::index_type Gauss2::NumberOfEvalPoints;
constexpr sp_grad::index_type Gauss2::NumberOfCollocationPoints;
constexpr sp_grad::index_type Gauss2::ridx[NumberOfCollocationPoints];
constexpr sp_grad::index_type Gauss2::sidx[NumberOfCollocationPoints];
constexpr sp_grad::index_type Gauss2::tidx[NumberOfCollocationPoints];
constexpr doublereal Gauss2::ri[NumberOfEvalPoints];
constexpr doublereal Gauss2::alphai[NumberOfEvalPoints];

constexpr sp_grad::index_type
Gauss2::iGetNumEvalPoints()
{
     return NumberOfCollocationPoints;
}

void
Gauss2::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < NumberOfCollocationPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < NumberOfEvalPoints);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < NumberOfEvalPoints);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < NumberOfEvalPoints);

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
     r(3) = ri[tidx[i]];
}

doublereal
Gauss2::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < NumberOfCollocationPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < NumberOfEvalPoints);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < NumberOfEvalPoints);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < NumberOfEvalPoints);

     return alphai[ridx[i]] * alphai[sidx[i]] * alphai[tidx[i]];
}

template <typename ElementType, typename CollocationType>
SolidElemStatic<ElementType, CollocationType>::SolidElemStatic(unsigned uLabel,
                                                               const std::array<const StructDispNodeAd*, ElementType::NumberOfNodes>& rgNodes,
                                                               std::array<std::unique_ptr<ConstitutiveLaw6D>, CollocationType::iGetNumEvalPoints()>&& rgCSL,
                                                               flag fOut)
     :Elem{uLabel, fOut}, ElemGravityOwner{uLabel, fOut}, rgNodes{rgNodes}, rgConstLaw(std::move(rgCSL))
{
     using namespace sp_grad;

     for (index_type i = 1; i <= ElementType::NumberOfNodes; ++i) {
          const Vec3& x0i = rgNodes[i - 1]->GetXCurr();
          for (index_type j = 1; j <= 3; ++j) {
               x0(j, i) = x0i(j);
          }
     }
}

template <typename ElementType, typename CollocationType>
SolidElemStatic<ElementType, CollocationType>::~SolidElemStatic()
{
}

template <typename ElementType, typename CollocationType>
Elem::Type SolidElemStatic<ElementType, CollocationType>::GetElemType(void) const
{
     return SOLID;
}

template <typename ElementType, typename CollocationType>
void
SolidElemStatic<ElementType, CollocationType>::SetValue(DataManager *pDM,
                                                        VectorHandler& X, VectorHandler& XP,
                                                        SimulationEntity::Hints *ph)
{
}

template <typename ElementType, typename CollocationType>
std::ostream& SolidElemStatic<ElementType, CollocationType>::Restart(std::ostream& out) const
{
     out << "## solid element: Restart not implemented yet\n";

     return out;
}

template <typename ElementType, typename CollocationType>
void SolidElemStatic<ElementType, CollocationType>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = ElementType::NumberOfNodes * 3;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType>::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                      doublereal dCoef,
                                                      const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                      const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                                      enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, ElementType::NumberOfNodes> u(3, ElementType::NumberOfNodes, 1);
     SpColVector<T, 3> Xj(3, 1);

     for (index_type j = 1; j <= ElementType::NumberOfNodes; ++j) {
          rgNodes[j - 1]->GetXCurr(Xj, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               u(i, j) = std::move(Xj(i));
               u(i, j) -= x0(i, j);
          }
     }

     SpColVectorA<doublereal, 3> r;
     SpMatrixA<doublereal, ElementType::NumberOfNodes, 3> hd;
     SpMatrixA<doublereal, 3, 3> J, invJ;
     doublereal detJ;
     SpMatrixA<doublereal, 6, ElementType::NumberOfNodes * 3> BL0;
     SpMatrixA<T, 6, ElementType::NumberOfNodes * 3, ElementType::NumberOfNodes * 3> BL1;
     SpColVectorA<T, 3 * ElementType::NumberOfNodes, 2 * 3 * ElementType::NumberOfNodes * CollocationType::iGetNumEvalPoints()> R;

     for (index_type iColloc = 0; iColloc < CollocationType::iGetNumEvalPoints(); ++iColloc) {
          CollocationType::GetPosition(iColloc, r);
          ElementType::ShapeFunctionDeriv(r, hd);
          Jacobian(hd, J);
          Inv(J, invJ, detJ);
          SpMatrix<doublereal, ElementType::NumberOfNodes, 3> h0d = hd * invJ;

          for (index_type k = 1; k <= ElementType::NumberOfNodes; ++k) {
               for (index_type i = 1; i <= 3; ++i) {
                    BL0(i, (k - 1) * 3 + i) = h0d(k, i);
               }

               BL0(4, (k - 1) * 3 + 1) = h0d(k, 2);
               BL0(4, (k - 1) * 3 + 2) = h0d(k, 1);
               BL0(5, (k - 1) * 3 + 2) = h0d(k, 3);
               BL0(5, (k - 1) * 3 + 3) = h0d(k, 2);
               BL0(6, (k - 1) * 3 + 1) = h0d(k, 3);
               BL0(6, (k - 1) * 3 + 3) = h0d(k, 1);
          }

          SpMatrix<T, 3, 3> dF = u * h0d;

          for (index_type k=1; k <= ElementType::NumberOfNodes; ++k) {
               for (index_type i = 1; i <= 3; ++i) {
                    for (index_type j = 1; j <= 3; ++j) {
                         BL1(i, (k - 1) * 3 + j) = dF(j , i) * h0d(k, i);
                    }
               }

               BL1(4, (k - 1) * 3 + 1) = dF(1, 1) * h0d(k, 2) + dF(1, 2) * h0d(k, 1);
               BL1(4, (k - 1) * 3 + 2) = dF(2, 1) * h0d(k, 2) + dF(2, 2) * h0d(k, 1);
               BL1(4, (k - 1) * 3 + 3) = dF(3, 1) * h0d(k, 2) + dF(3, 2) * h0d(k, 1);
               BL1(5, (k - 1) * 3 + 1) = dF(1, 2) * h0d(k, 3) + dF(1, 3) * h0d(k, 2);
               BL1(5, (k - 1) * 3 + 2) = dF(2, 2) * h0d(k, 3) + dF(2, 3) * h0d(k, 2);
               BL1(5, (k - 1) * 3 + 3) = dF(3, 2) * h0d(k, 3) + dF(3, 3) * h0d(k, 2);
               BL1(6, (k - 1) * 3 + 1) = dF(1, 1) * h0d(k, 3) + dF(1, 3) * h0d(k, 1);
               BL1(6, (k - 1) * 3 + 2) = dF(2, 1) * h0d(k, 3) + dF(2, 3) * h0d(k, 1);
               BL1(6, (k - 1) * 3 + 3) = dF(3, 1) * h0d(k, 3) + dF(3, 3) * h0d(k, 1);
          }

          const SpMatrix<T, 3, 3> G = 0.5 * (dF + Transpose(dF) + Transpose(dF) * dF);

          const SpColVector<T, 6> eps{G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1)};

          const SpColVector<T, 6> S = rgConstLaw[iColloc]->Update(eps);

#ifdef DEBUG
          if (std::is_same<T, doublereal>::value) {
               DEBUGCERR("eps=" << eps << "\n");
               DEBUGCERR("S=" << S << "\n");
          }
#endif
          const doublereal alpha = CollocationType::dGetWeight(iColloc);

          R -= (Transpose(BL0) * S + Transpose(BL1) * S) * alpha * detJ;
     }

     for (index_type i = 1; i <= ElementType::NumberOfNodes; ++i) {
          const index_type iEqIndex = rgNodes[i - 1]->iGetFirstMomentumIndex();

          for (index_type j = 1; j <= 3; ++j) {
               WorkVec.AddItem(iEqIndex + j, R((i - 1) * 3 + j));
          }
     }
}

template <typename ElementType, typename CollocationType>
SubVectorHandler&
SolidElemStatic<ElementType, CollocationType>::AssRes(SubVectorHandler& WorkVec,
                                                      doublereal dCoef,
                                                      const VectorHandler& XCurr,
                                                      const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemStatic::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType>
VariableSubMatrixHandler&
SolidElemStatic<ElementType, CollocationType>::AssJac(VariableSubMatrixHandler& WorkMat,
                                                      doublereal dCoef,
                                                      const VectorHandler& XCurr,
                                                      const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemStatic::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

template <typename ElementType, typename CollocationType>
void
SolidElemStatic<ElementType, CollocationType>::AssJac(VectorHandler& JacY,
                                                      const VectorHandler& Y,
                                                      doublereal dCoef,
                                                      const VectorHandler& XCurr,
                                                      const VectorHandler& XPrimeCurr,
                                                      VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

template <typename ElementType, typename CollocationType>
void
SolidElemStatic<ElementType, CollocationType>::Jacobian(const sp_grad::SpMatrix<doublereal, ElementType::NumberOfNodes, 3>& hd,
                                                        sp_grad::SpMatrix<doublereal, 3, 3>& J)
{
     using namespace sp_grad;

     J = Transpose(x0 * hd);
}

template <typename SolidElemType>
Elem*
ReadSolid(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
        DEBUGCOUTFNAME("ReadSolid");

        using namespace sp_grad;

        constexpr index_type iNumNodes = SolidElemType::iGetNumNodes();
        constexpr index_type iNumEvalPoints = SolidElemType::iGetNumEvalPoints();

        std::array<const StructDispNodeAd*, iNumNodes> rgNodes;

        for (index_type i = 0; i < iNumNodes; ++i) {
             rgNodes[i] = pDM->ReadNode<const StructDispNodeAd, Node::STRUCTURAL>(HP);
        }

        std::array<std::unique_ptr<ConstitutiveLaw6D>, iNumEvalPoints> rgConstLaw;

        for (index_type i = 0; i < iNumEvalPoints; ++i) {
             ConstLawType::Type CLType_I = ConstLawType::UNKNOWN;

             rgConstLaw[i].reset(HP.GetConstLaw6D(CLType_I));

             if (rgConstLaw[i]->iGetNumDof() != 0) {
                  silent_cerr("line " << HP.GetLineData()
                              << ": solid(" << uLabel << ") "
                              "does not support dynamic constitutive laws yet"
                              << std::endl);
                  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
             }

             if (CLType_I != ConstLawType::ELASTIC) {
                  silent_cerr("solid(" << uLabel << ") does not support this constitutive law type "
                              "at line " << HP.GetLineData() << "\n");
                  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
             }
        }

        const flag fOut = pDM->fReadOutput(HP, Elem::SOLID);

        std::ostream& out = pDM->GetLogFile();

        out << "hexahedron8: " << uLabel;

        for (sp_grad::index_type i = 0; i < iNumNodes; ++i) {
             out << ',' << rgNodes[i]->GetLabel();
        }

        out << "\n";

        Elem* pEl = nullptr;

        SAFENEWWITHCONSTRUCTOR(pEl,
                               SolidElemType,
                               SolidElemType(uLabel,
                                             rgNodes,
                                             std::move(rgConstLaw),
                                             fOut));

        if (HP.IsArg()) {
                silent_cerr("semicolon expected "
                        "at line " << HP.GetLineData() << std::endl);
                throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return pEl;
}

template Elem* ReadSolid<SolidElemStatic<Hexahedron8, Gauss2>>(DataManager*, MBDynParser&, unsigned int);
