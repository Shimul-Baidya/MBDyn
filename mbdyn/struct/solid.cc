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

constexpr sp_grad::index_type Hexahedron8::iNumNodes;

void
Hexahedron8::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
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

constexpr sp_grad::index_type Gauss2::iGaussOrder;
constexpr sp_grad::index_type Gauss2::iNumEvalPoints;
constexpr sp_grad::index_type Gauss2::ridx[iNumEvalPoints];
constexpr sp_grad::index_type Gauss2::sidx[iNumEvalPoints];
constexpr sp_grad::index_type Gauss2::tidx[iNumEvalPoints];
constexpr doublereal Gauss2::ri[iGaussOrder];
constexpr doublereal Gauss2::alphai[iGaussOrder];

void
Gauss2::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
     r(3) = ri[tidx[i]];
}

doublereal
Gauss2::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     return alphai[ridx[i]] * alphai[sidx[i]] * alphai[tidx[i]];
}

template <typename ElementType, typename CollocationType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType>::iNumNodes;

template <typename ElementType, typename CollocationType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType>::iNumEvalPoints;

template <typename ElementType, typename CollocationType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType>::iNumDof;

template <typename ElementType, typename CollocationType>
SolidElemStatic<ElementType, CollocationType>::SolidElemStatic(unsigned uLabel,
                                                               const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                                                               std::array<std::unique_ptr<ConstitutiveLaw6D>, iNumEvalPoints>&& rgCSL,
                                                               flag fOut)
     :Elem{uLabel, fOut}, ElemGravityOwner{uLabel, fOut}, rgNodes{rgNodes}
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumNodes; ++i) {
          const Vec3& x0i = rgNodes[i - 1]->GetXCurr();
          for (index_type j = 1; j <= 3; ++j) {
               x0(j, i) = x0i(j);
          }
     }

     for (index_type iColloc = 0; iColloc < iNumEvalPoints; ++iColloc) {
          rgCollocData[iColloc].Init(iColloc, x0, std::move(rgCSL[iColloc]));
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
void SolidElemStatic<ElementType, CollocationType>::Output(OutputHandler& OH) const
{
     using namespace sp_grad;

     if (bToBeOutput()) {
          sp_grad::SpMatrixA<doublereal, 6, iNumEvalPoints> epsilon;
          sp_grad::SpMatrixA<doublereal, 6, iNumEvalPoints> tau;
          
          ComputeStressStrain(epsilon, tau);

          if (OH.UseText(OutputHandler::SOLIDS)) {
               std::ostream& of = OH.Solids();

               of << std::setw(8) << GetLabel();

               for (index_type j = 1; j <= iNumEvalPoints; ++j) {
                    for (index_type i = 1; i <= 6; ++i) {
                         of << ' ' << epsilon(i, j);
                    }

                    for (index_type i = 1; i <= 6; ++i) {
                         of << ' ' << tau(i, j);
                    }
               }

               of << '\n';
          }
     }
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
     *piNumRows = iNumDof;
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

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1);

     GetDeformations(u, dCoef, func);

     constexpr index_type iNumColsR = 2 * iNumDof * iNumEvalPoints;
     
     SpColVectorA<T, iNumDof, iNumColsR> R;

     AssStiffnessVec(u, R, dCoef, func);

     for (index_type i = 1; i <= iNumNodes; ++i) {
          const index_type iEqIndex = rgNodes[i - 1]->iGetFirstMomentumIndex();

          for (index_type j = 1; j <= 3; ++j) {
               WorkVec.AddItem(iEqIndex + j, R((i - 1) * 3 + j));
          }
     }
}

template <typename ElementType, typename CollocationType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType>::AssStiffnessVec(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                               sp_grad::SpColVector<T, iNumDof>& R,
                                                               const doublereal dCoef,
                                                               const sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;
     
     SpMatrixA<T, 6, iNumDof, iNumDof> BL1;

     for (index_type iColloc = 0; iColloc < iNumEvalPoints; ++iColloc) {
          const auto& h0d = rgCollocData[iColloc].h0d;

          const SpMatrix<T, 3, 3> dF = u * h0d;

          const SpMatrix<T, 3, 3> G = 0.5 * (dF + Transpose(dF) + Transpose(dF) * dF);

          const SpColVector<T, 6> sigma = rgCollocData[iColloc].ComputeStress(G, dF);

          const doublereal alpha = CollocationType::dGetWeight(iColloc);

          rgCollocData[iColloc].StrainMatrix1(dF, BL1);
          
          const auto& BL0 = rgCollocData[iColloc].BL0;
          
          R -= (Transpose(BL0) * sigma + Transpose(BL1) * sigma) * (alpha * rgCollocData[iColloc].detJ);
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
template <typename T>
void
SolidElemStatic<ElementType, CollocationType>::GetDeformations(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                               const doublereal dCoef,
                                                               const sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     SpColVectorA<T, 3, 1> Xj;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          rgNodes[j - 1]->GetXCurr(Xj, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               u(i, j) = std::move(Xj(i));
               u(i, j) -= x0(i, j);
          }
     }
}

template <typename ElementType, typename CollocationType>
void
SolidElemStatic<ElementType, CollocationType>::ComputeStressStrain(sp_grad::SpMatrixA<doublereal, 6, iNumEvalPoints>& epsilon,
                                                                   sp_grad::SpMatrixA<doublereal, 6, iNumEvalPoints>& tau) const
{
     using namespace sp_grad;

     static constexpr index_type idxr[6] = {1, 2, 3, 1, 2, 3};
     static constexpr index_type idxc[6] = {1, 2, 3, 2, 3, 1};

     SpMatrixA<doublereal, 3, 3> S;

     for (index_type iColloc = 0; iColloc < iNumEvalPoints; ++iColloc) {
          const auto& sigma = rgCollocData[iColloc].sigma;

          for (index_type i = 1; i <= 6; ++i) {
               S(idxr[i - 1], idxc[i - 1]) = sigma(i);
          }

          for (index_type i = 4; i <= 6; ++i) {
               S(idxc[i - 1], idxr[i - 1]) = sigma(i);
          }

          const auto& F = rgCollocData[iColloc].F;

          const SpMatrixA<doublereal, 3, 3> taui = F * S * Transpose(F) / Det(F);

          const auto& G = rgCollocData[iColloc].G;

          for (index_type i = 1; i <= 3; ++i) {
               epsilon(i, iColloc + 1) = sqrt(1. + 2. * G(i, i)) - 1.;
          }

          for (index_type i = 4; i <= 6; ++i) {
               epsilon(i, iColloc + 1) = 2. * G(idxr[i - 1], idxc[i - 1]) / ((1. + epsilon(idxr[i - 1], iColloc + 1)) * (1. + epsilon(idxc[i - 1], iColloc + 1)));
          }

          for (index_type i = 1; i <= 6; ++i) {
               tau(i, iColloc + 1) = taui(idxr[i - 1], idxc[i - 1]);
          }
     }
}

template <typename ElementType, typename CollocationType>
void
SolidElemStatic<ElementType, CollocationType>::CollocData::Init(const sp_grad::index_type iColloc, const sp_grad::SpMatrixA<doublereal, 3, iNumNodes>& x0, std::unique_ptr<ConstitutiveLaw6D>&& pCSL)
{
     using namespace sp_grad;

     pConstLaw = std::move(pCSL);

     SpColVectorA<doublereal, 3> r;
     sp_grad::SpMatrixA<doublereal, ElementType::iNumNodes, 3> hd;

     CollocationType::GetPosition(iColloc, r);
     ElementType::ShapeFunctionDeriv(r, hd);

     J = Transpose(x0 * hd);
     Inv(J, invJ, detJ);

     h0d = hd * invJ;

     for (index_type k = 1; k <= iNumNodes; ++k) {
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
}

template <typename ElementType, typename CollocationType>
template <typename T>
sp_grad::SpColVector<T, 6>
SolidElemStatic<ElementType, CollocationType>::CollocData::ComputeStress(const sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpMatrix<T, 3, 3>& dF)
{
     using namespace sp_grad;
     const SpColVector<T, 6> eps{G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1)};
     const SpColVector<T, 6> sigma = pConstLaw->Update(eps);

     UpdateStressStrain(G, sigma, dF);

     return sigma;
}

template <typename ElementType, typename CollocationType>
inline void
SolidElemStatic<ElementType, CollocationType>::CollocData::UpdateStressStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& G_tmp,
                                                                              const sp_grad::SpColVector<doublereal, 6>& sigma_tmp,
                                                                              const sp_grad::SpMatrix<doublereal, 3, 3>& dF_tmp)
{
     using namespace sp_grad;
               
     G = G_tmp;
     sigma = sigma_tmp;
     F = dF_tmp + Eye3;
}

template <typename ElementType, typename CollocationType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType>::CollocData::StrainMatrix1(const sp_grad::SpMatrix<T, 3, 3>& dF,
                                                                         sp_grad::SpMatrix<T, 6, iNumDof>& BL1) const
{
     using namespace sp_grad;

     for (index_type k = 1; k <= iNumNodes; ++k) {
          for (index_type i = 1; i <= 3; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    BL1(i, (k - 1) * 3 + j) = dF(j , i) * h0d(k, i);
               }
          }

          static constexpr index_type idx1[] = {2, 3, 3};
          static constexpr index_type idx2[] = {1, 2, 1};

          for (index_type i = 1; i <= 3; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    BL1(i + 3, (k - 1) * 3 + j) = dF(j, idx2[i - 1]) * h0d(k, idx1[i - 1]) + dF(j, idx1[i - 1]) * h0d(k, idx2[i - 1]);
               }
          }
     }
}

template <typename SolidElemType>
Elem*
ReadSolid(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
        DEBUGCOUTFNAME("ReadSolid");

        using namespace sp_grad;

        constexpr index_type iNumNodes = SolidElemType::iNumNodes;
        constexpr index_type iNumEvalPoints = SolidElemType::iNumEvalPoints;

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
