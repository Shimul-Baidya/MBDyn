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
#include "rbdispjad.h"

RigidBodyDispJointAd::RigidBodyDispJointAd(unsigned int uL,
                                           const DofOwner* pD,
                                           const StructNodeAd* pNodeMasterTmp,
                                           std::vector<SlaveNodeData>&& rgNodesSlaveTmp,
                                           flag fOut)
     :Elem(uL, fOut),
      Joint(uL, pD, fOut),
      pNodeMaster{pNodeMasterTmp},
      rgNodesSlave{std::move(rgNodesSlaveTmp)},
      FmTmp(::Zero3),
      MmTmp(::Zero3)
{
}

RigidBodyDispJointAd::~RigidBodyDispJointAd()
{
}

void RigidBodyDispJointAd::Output(OutputHandler& OH) const
{
     using namespace sp_grad;

     if (OH.UseText(OutputHandler::JOINTS)) {
          const Mat3x3& Rm = pNodeMaster->GetRCurr();
          Joint::Output(OH.Joints(), "OffsetDispJoint", GetLabel(), -Rm.MulTV(FmTmp), -Rm.MulTV(MmTmp), -FmTmp, -MmTmp) << '\n';
     }
}

void RigidBodyDispJointAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12 + rgNodesSlave.size() * 3;
     *piNumCols = 0;
}

unsigned int RigidBodyDispJointAd::iGetNumDof(void) const
{
     return 6u;
}

DofOrder::Order RigidBodyDispJointAd::GetDofType(unsigned int i) const
{
     return DofOrder::ALGEBRAIC;
}

DofOrder::Order RigidBodyDispJointAd::GetEqType(unsigned int i) const
{
     return DofOrder::ALGEBRAIC;
}

std::ostream& RigidBodyDispJointAd::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 3 << ": position constraint [lambda1, lambda2, lambda3]\n";
     out << prefix << iFirstIndex + 4 << "->" << iFirstIndex + 6 << ": orientation constraint [lambda4, lambda5, lambda6]\n";

     if (bInitial) {
          out << prefix << iFirstIndex + 7 << "->" << iFirstIndex + 9 << ": position constraint derivatives [lambdaP1, lambdaP2, lambdaP3]\n";
          out << prefix << iFirstIndex + 10 << "->" << iFirstIndex + 12 << ": orientation constraint derivatives [lambdaP4, lambdaP5, lambdaP6]\n";
     }

     return out;
}

std::ostream& RigidBodyDispJointAd::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
     const integer iFirstIndex = iGetFirstIndex();

     out << prefix << iFirstIndex + 1 << "->" << iFirstIndex + 3 << ": position constraints [Phi1, Phi2, Phi3]\n";
     out << prefix << iFirstIndex + 4 << "->" << iFirstIndex + 6 << ": orientation constraints [Phi4, Phi5, Phi6]\n";

     if (bInitial) {
          out << prefix << iFirstIndex + 7 << "->" << iFirstIndex + 9 << ": position constraints derivatives [PhiP1, PhiP2, PhiP3]\n";
          out << prefix << iFirstIndex + 10 << "->" << iFirstIndex + 12 << ": orientation constraints [PhiP4, PhiP5, PhiP6]\n";
     }

     return out;
}

VariableSubMatrixHandler&
RigidBodyDispJointAd::AssJac(VariableSubMatrixHandler& WorkMat,
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
RigidBodyDispJointAd::AssJac(VectorHandler& JacY,
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
RigidBodyDispJointAd::AssRes(SubVectorHandler& WorkVec,
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

unsigned int RigidBodyDispJointAd::iGetNumPrivData(void) const
{
     return 6u;
}

unsigned int RigidBodyDispJointAd::iGetPrivDataIdx(const char *s) const
{
     static constexpr char rgPrivDataName[][3] = {"Fx", "Fy", "Fz", "Mx", "My", "Mz"};

     for (integer i = 0; i < 3; ++i) {
          if (0 == strcmp(rgPrivDataName[i], s)) {
               return i + 1;
          }
     }

     return 0u;
}

doublereal RigidBodyDispJointAd::dGetPrivData(unsigned int i) const
{
     switch (i) {
     case 1:
     case 2:
     case 3:
          return -FmTmp(i);

     case 4:
     case 5:
     case 6:
          return -MmTmp(i - 3);

     default:
          ASSERT(0);
          return 0.;
     }
}

int RigidBodyDispJointAd::GetNumConnectedNodes() const
{
     return 1 + rgNodesSlave.size();
}

void RigidBodyDispJointAd::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
     connectedNodes.resize(0);
     connectedNodes.reserve(GetNumConnectedNodes());

     connectedNodes.push_back(pNodeMaster);

     for (const auto& oNDS: rgNodesSlave) {
          connectedNodes.push_back(oNDS.pNode);
     }
}


void RigidBodyDispJointAd::SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                                    SimulationEntity::Hints *ph)
{
}

std::ostream& RigidBodyDispJointAd::Restart(std::ostream& out) const
{
     return out;
}

unsigned int RigidBodyDispJointAd::iGetInitialNumDof(void) const
{
     return 12;
}

void
RigidBodyDispJointAd::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 6 + rgNodesSlave.size() * 3;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
RigidBodyDispJointAd::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                    const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::InitialAssJac(this,
                                                                   WorkMat.SetSparseGradient(),
                                                                   XCurr,
                                                                   sp_grad::SpFunctionCall::INITIAL_ASS_JAC);

     return WorkMat;
}

SubVectorHandler&
RigidBodyDispJointAd::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
     sp_grad::SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                          WorkVec,
                                                          XCurr,
                                                          sp_grad::SpFunctionCall::INITIAL_ASS_RES);

     return WorkVec;
}

template <typename T>
void RigidBodyDispJointAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                  doublereal dCoef,
                                  const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                  const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                  sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const integer iFirstIndexNodeMaster = pNodeMaster->iGetFirstMomentumIndex();
     const integer iFirstIndexLambda = iGetFirstIndex();

     SpColVector<T, 3> Xm(3, 1), Xj(3, 1);
     SpColVector<T, 3> Phit(3, rgNodesSlave.size() * 5);
     SpColVector<T, 3> Phir(3, rgNodesSlave.size() * 10);
     SpColVector<T, 3> Fm(3, rgNodesSlave.size() * 9);
     SpColVector<T, 3> Mm(3, rgNodesSlave.size() * 30);
     SpMatrix<T, 3, 3> Rm(3, 3, 3);
     SpColVector<T, 3> lambdat(3, 1), lambdar(3, 1);

     pNodeMaster->GetXCurr(Xm, dCoef, func);
     pNodeMaster->GetRCurr(Rm, dCoef, func);
     XCurr.GetVec(iFirstIndexLambda + 1, lambdat, 1.);
     XCurr.GetVec(iFirstIndexLambda + 4, lambdar, 1.);

     for (const auto& oNDS: rgNodesSlave) {
          const integer iFirstIndexSlave = oNDS.pNode->iGetFirstMomentumIndex();
          const Vec3& oj = oNDS.offset;
          const doublereal wj = oNDS.weight;

          oNDS.pNode->GetXCurr(Xj, dCoef, func);

          const SpColVector<T, 3> Rm_oj = Rm * oj;
          const SpColVector<T, 3> DeltaXj = Xm - Xj;
          const SpColVector<T, 3> Fj = (Cross(Rm_oj, lambdar) - lambdat) * wj;

          ASSERT(Fj.iGetMaxSize() <= 9);

          WorkVec.AddItem(iFirstIndexSlave + 1, Fj);

          Phit += (DeltaXj + Rm_oj) * (wj / dCoef);
          Phir += Cross(Rm_oj, DeltaXj) * (wj / dCoef);

          Fm -= Fj;
          Mm += Cross(Rm_oj, lambdat + Cross(DeltaXj, lambdar)) * wj;
     }

     ASSERT(Phit.iGetMaxSize() <= rgNodesSlave.size() * 5);
     ASSERT(Phir.iGetMaxSize() <= rgNodesSlave.size() * 10);
     ASSERT(Fm.iGetMaxSize() <= rgNodesSlave.size() * 9);
     ASSERT(Mm.iGetMaxSize() <= rgNodesSlave.size() * 30);

     DEBUGCERR("Phit.iGetMaxSize() / N = " << Phit.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("Phir.iGetMaxSize() / N = " << Phir.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("Fm.iGetMaxSize() / N = " << Fm.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("Mm.iGetMaxSize() / N = " << Mm.iGetMaxSize() / rgNodesSlave.size() << "\n");

     WorkVec.AddItem(iFirstIndexNodeMaster + 1, Fm);
     WorkVec.AddItem(iFirstIndexNodeMaster + 4, Mm);
     WorkVec.AddItem(iFirstIndexLambda + 1, Phit);
     WorkVec.AddItem(iFirstIndexLambda + 4, Phir);

     SaveReactionForce(Fm, Mm);
}

template <typename T>
void RigidBodyDispJointAd::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                         const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                         sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const integer iFirstIndexNodeMaster = pNodeMaster->iGetFirstPositionIndex();
     const integer iFirstIndexLambda = iGetFirstIndex();

     SpColVector<T, 3> Xm(3, 1), Xj(3, 1);
     SpColVector<T, 3> XPm(3, 1), XPj(3, 1);
     SpColVector<T, 3> Wm(3, 3);
     SpColVector<T, 3> Phit(3, rgNodesSlave.size() * 5);
     SpColVector<T, 3> Phir(3, rgNodesSlave.size() * 10);
     SpColVector<T, 3> PhiPt(3, rgNodesSlave.size() * 10);
     SpColVector<T, 3> PhiPr(3, rgNodesSlave.size() * 30);
     SpColVector<T, 3> Fm(3, rgNodesSlave.size() * 9);
     SpColVector<T, 3> Mm(3, rgNodesSlave.size() * 30);

     SpColVector<T, 3> FPm(3, rgNodesSlave.size() * 27);
     SpColVector<T, 3> MPm(3, rgNodesSlave.size() * 88);
     SpMatrix<T, 3, 3> Rm(3, 3, 3);
     SpColVector<T, 3> lambdat(3, 1), lambdar(3, 1);
     SpColVector<T, 3> lambdaPt(3, 1), lambdaPr(3, 1);

     pNodeMaster->GetXCurr(Xm, 1., func);
     pNodeMaster->GetVCurr(XPm, 1., func);
     pNodeMaster->GetRCurr(Rm, 1., func);
     pNodeMaster->GetWCurr(Wm, 1., func);

     XCurr.GetVec(iFirstIndexLambda + 1, lambdat, 1.);
     XCurr.GetVec(iFirstIndexLambda + 4, lambdar, 1.);
     XCurr.GetVec(iFirstIndexLambda + 7, lambdaPt, 1.);
     XCurr.GetVec(iFirstIndexLambda + 10, lambdaPr, 1.);

     for (const auto& oNDS: rgNodesSlave) {
          const integer iFirstIndexSlave = oNDS.pNode->iGetFirstPositionIndex();
          const Vec3& oj = oNDS.offset;
          const doublereal wj = oNDS.weight;

          oNDS.pNode->GetXCurr(Xj, 1., func);
          oNDS.pNode->GetVCurr(XPj, 1., func);

          const SpColVector<T, 3> Rm_oj = Rm * oj;
          const SpColVector<T, 3> DeltaXj = Xm - Xj;
          const SpColVector<T, 3> DeltaXPj = XPm - XPj;
          const SpColVector<T, 3> WmRm_oj = Cross(Wm, Rm_oj);

          const SpColVector<T, 3> Fj = (Cross(Rm_oj, lambdar) - lambdat) * wj;
          const SpColVector<T, 3> FPj = (Cross(Rm_oj, lambdaPr) - Cross(lambdar, WmRm_oj) - lambdaPt) * wj;

          ASSERT(Fj.iGetMaxSize() <= 9);

          WorkVec.AddItem(iFirstIndexSlave + 1, Fj);
          WorkVec.AddItem(iFirstIndexSlave + 4, FPj);

          Phit += (DeltaXj + Rm_oj) * wj;
          Phir += Cross(Rm_oj, DeltaXj) * wj;
          PhiPt += (XPm + WmRm_oj - XPj) * wj;
          PhiPr += (Cross(Rm_oj, DeltaXPj) - Cross(DeltaXj, WmRm_oj)) * wj;
          Fm -= Fj;
          Mm += Cross(Rm_oj, lambdat + Cross(DeltaXj, lambdar)) * wj;
          FPm -= FPj;
          MPm += (Cross(Rm_oj, lambdaPt + Cross(DeltaXPj, lambdar) + Cross(DeltaXj, lambdaPr))
                  - Cross(lambdat + Cross(DeltaXj, lambdar), WmRm_oj)) * wj;
     }

     ASSERT(Phit.iGetMaxSize() <= rgNodesSlave.size() * 5);
     ASSERT(Phir.iGetMaxSize() <= rgNodesSlave.size() * 10);
     ASSERT(PhiPt.iGetMaxSize() <= rgNodesSlave.size() * 10);
     ASSERT(PhiPr.iGetMaxSize() <= rgNodesSlave.size() * 30);
     ASSERT(Fm.iGetMaxSize() <= rgNodesSlave.size() * 9);
     ASSERT(Mm.iGetMaxSize() <= rgNodesSlave.size() * 30);
     ASSERT(FPm.iGetMaxSize() <= rgNodesSlave.size() * 27);
     ASSERT(MPm.iGetMaxSize() <= rgNodesSlave.size() * 88);

     DEBUGCERR("Phit.iGetMaxSize() / N = " << Phit.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("Phir.iGetMaxSize() / N = " << Phir.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("PhiPt.iGetMaxSize() / N = " << PhiPt.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("PhiPr.iGetMaxSize() / N = " << PhiPr.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("Fm.iGetMaxSize() / N = " << Fm.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("Mm.iGetMaxSize() / N = " << Mm.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("FPm.iGetMaxSize() / N = " << FPm.iGetMaxSize() / rgNodesSlave.size() << "\n");
     DEBUGCERR("MPm.iGetMaxSize() / N = " << MPm.iGetMaxSize() / rgNodesSlave.size() << "\n");

     WorkVec.AddItem(iFirstIndexNodeMaster + 1, Fm);
     WorkVec.AddItem(iFirstIndexNodeMaster + 4, Mm);
     WorkVec.AddItem(iFirstIndexNodeMaster + 7, FPm);
     WorkVec.AddItem(iFirstIndexNodeMaster + 10, MPm);
     WorkVec.AddItem(iFirstIndexLambda + 1, Phit);
     WorkVec.AddItem(iFirstIndexLambda + 4, Phir);
     WorkVec.AddItem(iFirstIndexLambda + 7, PhiPt);
     WorkVec.AddItem(iFirstIndexLambda + 10, PhiPr);

     SaveReactionForce(Fm, Mm);
}

void RigidBodyDispJointAd::SaveReactionForce(const sp_grad::SpColVector<doublereal, 3>& Fm, const sp_grad::SpColVector<doublereal, 3>& Mm)
{
     FmTmp = Fm;
     MmTmp = Mm;
}

const OutputHandler::Dimensions
RigidBodyDispJointAd::GetEquationDimension(integer index) const {
     switch (index) {
     case 1:
     case 2:
     case 3:
          return OutputHandler::Dimensions::Length;
     case 4:
     case 5:
     case 6:
          return OutputHandler::Dimensions::UnknownDimension;
     default:
          ASSERT(0);
          return OutputHandler::Dimensions::UnknownDimension;
     }
}

Joint::Type RigidBodyDispJointAd::GetJointType() const
{
     return RIGIDBODYDISPLACEMENTJOINT;
}
