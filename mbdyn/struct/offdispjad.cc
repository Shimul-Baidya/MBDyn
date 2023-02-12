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
  Copyright (C) 2023(-2023) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <cstring>
#include "output.h"
#include "offdispjad.h"

OffsetDispJointAd::OffsetDispJointAd(unsigned int uL,
                                     const DofOwner* pD,
                                     const StructNodeAd* pNode1,
                                     const Vec3& o1,
                                     const StructDispNodeAd* pNode2,
                                     flag fOut)
     :Elem(uL, fOut),
      Joint(uL, pD, fOut),
      pNode1(pNode1),
      pNode2(pNode2),
      o1(o1),
      F1Tmp(::Zero3),
      M1Tmp(::Zero3)
{
}

OffsetDispJointAd::~OffsetDispJointAd()
{
}

void OffsetDispJointAd::Output(OutputHandler& OH) const
{
     using namespace sp_grad;

     if (bToBeOutput() && OH.UseText(OutputHandler::JOINTS)) {
          const Mat3x3& R1 = pNode1->GetRCurr();
          Joint::Output(OH.Joints(), "OffsetDispJoint", GetLabel(), -R1.MulTV(F1Tmp), -R1.MulTV(M1Tmp), -F1Tmp, -M1Tmp) << '\n';
     }
}

void OffsetDispJointAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12;
     *piNumCols = 0;
}

unsigned int OffsetDispJointAd::iGetNumDof(void) const
{
     return 3u;
}

DofOrder::Order OffsetDispJointAd::GetDofType(unsigned int i) const
{
     return DofOrder::ALGEBRAIC;
}

DofOrder::Order OffsetDispJointAd::GetEqType(unsigned int i) const
{
     return DofOrder::ALGEBRAIC;
}

std::ostream& OffsetDispJointAd::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 3 << ": reaction forces [Fx, Fy, Fz]\n";

     if (bInitial) {
          out << prefix << iFirstIndex + 4 << "->" << iFirstIndex + 6 << ": reaction force derivatives [FPx, FPy, FPz]\n";
     }

     return out;
}

std::ostream& OffsetDispJointAd::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 3 << ": position constraints [Px, Py, Pz]\n";

     if (bInitial) {
          out << prefix << iFirstIndex + 4 << "->" << iFirstIndex + 6 << ": velocity constraints [vx, vy, vz]\n";
     }

     return out;
}

VariableSubMatrixHandler&
OffsetDispJointAd::AssJac(VariableSubMatrixHandler& WorkMat,
                          doublereal dCoef,
                          const VectorHandler& XCurr,
                          const VectorHandler& XPrimeCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                   WorkMat.SetSparseGradient(),
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
OffsetDispJointAd::AssJac(VectorHandler& JacY,
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

SubVectorHandler&
OffsetDispJointAd::AssRes(SubVectorHandler& WorkVec,
                          doublereal dCoef,
                          const VectorHandler& XCurr,
                          const VectorHandler& XPrimeCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

unsigned int OffsetDispJointAd::iGetNumPrivData(void) const
{
     return 3u;
}

unsigned int OffsetDispJointAd::iGetPrivDataIdx(const char *s) const
{
     static constexpr char rgPrivDataName[][3] = {"Fx", "Fy", "Fz", "Mx", "My", "Mz"};

     for (integer i = 0; i < 3; ++i) {
          if (0 == strcmp(rgPrivDataName[i], s)) {
               return i + 1;
          }
     }

     return 0u;
}

doublereal OffsetDispJointAd::dGetPrivData(unsigned int i) const
{
     switch (i) {
     case 1:
     case 2:
     case 3:
          return -F1Tmp(i);
     case 4:
     case 5:
     case 6:
          return -M1Tmp(i - 3);
     default:
          ASSERT(0);
          return 0.;
     }
}

int OffsetDispJointAd::GetNumConnectedNodes(void) const
{
     return 2;
}

void OffsetDispJointAd::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
}

void OffsetDispJointAd::SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                                 SimulationEntity::Hints *ph)
{
}

std::ostream& OffsetDispJointAd::Restart(std::ostream& out) const
{
     return out;
}

unsigned int OffsetDispJointAd::iGetInitialNumDof(void) const
{
     return 6u;
}

void
OffsetDispJointAd::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 24;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
OffsetDispJointAd::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                 const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::InitialAssJac(this,
                                                                   WorkMat.SetSparseGradient(),
                                                                   XCurr,
                                                                   sp_grad::SpFunctionCall::INITIAL_ASS_JAC);

     return WorkMat;
}

SubVectorHandler&
OffsetDispJointAd::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                          WorkVec,
                                                          XCurr,
                                                          sp_grad::SpFunctionCall::INITIAL_ASS_RES);

     return WorkVec;
}

template <typename T>
void OffsetDispJointAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                               doublereal dCoef,
                               const sp_grad::SpGradientVectorHandler<T>& XCurr,
                               const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                               sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const integer iFirstIndexNode1 = pNode1->iGetFirstMomentumIndex();
     const integer iFirstIndexNode2 = pNode2->iGetFirstMomentumIndex();
     const integer iFirstIndexLambda = iGetFirstIndex();

     SpColVector<T, 3> X1(3, 1), X2(3, 1), F2(3, 1);
     SpMatrix<T, 3, 3> R1(3, 3, 3);

     pNode1->GetXCurr(X1, dCoef, func);
     pNode1->GetRCurr(R1, dCoef, func);
     pNode2->GetXCurr(X2, dCoef, func);

     XCurr.GetVec(iFirstIndexLambda + 1, F2, 1.);

     const SpColVector<T, 3> R1_o1 = R1 * o1;
     const SpColVector<T, 3> Phi = (X1 + R1_o1 - X2) / dCoef;

     const SpColVector<T, 3> F1 = -F2;
     const SpColVector<T, 3> M1 = Cross(R1_o1, F1);

     SaveReactionForce(F1, M1);

     WorkVec.AddItem(iFirstIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstIndexLambda + 1, Phi);
}

template <typename T>
void OffsetDispJointAd::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                      const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                      sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const integer iFirstIndexNode1 = pNode1->iGetFirstPositionIndex();
     const integer iFirstIndexNode2 = pNode2->iGetFirstPositionIndex();
     const integer iFirstIndexLambda = iGetFirstIndex();

     SpColVector<T, 3> X1(3, 1), X1P(3, 1), X2(3, 1), X2P(3, 1), F2(3, 1), F2P(3, 1), W1(3, 3);
     SpMatrix<T, 3, 3> R1(3, 3, 3);

     pNode1->GetXCurr(X1, 1., func);
     pNode1->GetVCurr(X1P, 1., func);
     pNode1->GetRCurr(R1, 1., func);
     pNode1->GetWCurr(W1, 1., func);
     pNode2->GetXCurr(X2, 1., func);
     pNode2->GetVCurr(X2P, 1., func);

     for (index_type i = 1; i <= 3; ++i) {
          XCurr.dGetCoef(iFirstIndexLambda + i, F2(i), 1.);
          XCurr.dGetCoef(iFirstIndexLambda + i + 3, F2P(i), 1.);
     }

     const SpColVector<T, 3> R1_o1 = R1 * o1;
     const SpColVector<T, 3> Phi = X1 + R1_o1 - X2;
     const SpColVector<T, 3> PhiP = X1P + Cross(W1, R1_o1) - X2P;
     const SpColVector<T, 3> F1 = -F2;
     const SpColVector<T, 3> M1 = Cross(R1_o1, F1);
     const SpColVector<T, 3> F1P = -F2P;
     const SpColVector<T, 3> M1P = Cross(R1_o1, F1P) + Cross(Cross(W1, R1_o1), F1);

     SaveReactionForce(F1, M1);

     WorkVec.AddItem(iFirstIndexNode1 + 1, F1);
     WorkVec.AddItem(iFirstIndexNode1 + 4, M1);
     WorkVec.AddItem(iFirstIndexNode1 + 7, F1P);
     WorkVec.AddItem(iFirstIndexNode1 + 10, M1P);
     WorkVec.AddItem(iFirstIndexNode2 + 1, F2);
     WorkVec.AddItem(iFirstIndexNode2 + 4, F2P);
     WorkVec.AddItem(iFirstIndexLambda + 1, Phi);
     WorkVec.AddItem(iFirstIndexLambda + 4, PhiP);
}

void OffsetDispJointAd::SaveReactionForce(const sp_grad::SpColVector<doublereal, 3>& F1, const sp_grad::SpColVector<doublereal, 3>& M1)
{
     F1Tmp = F1;
     M1Tmp = M1;
}

const OutputHandler::Dimensions
OffsetDispJointAd::GetEquationDimension(integer index) const {
     switch (index) {
     case 1:
     case 2:
     case 3:
          return OutputHandler::Dimensions::Length;
     case 4:
     case 5:
     case 6:
          return OutputHandler::Dimensions::Velocity;
     default:
          ASSERT(0);
          return OutputHandler::Dimensions::UnknownDimension;
     }
}

Joint::Type OffsetDispJointAd::GetJointType() const
{
     return OFFSETDISPLACEMENTJOINT;
}
