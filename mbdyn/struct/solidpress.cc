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

#include <array>
#include <memory>

#include "presnodead.h"
#include "strnodead.h"
#include "sp_matvecass.h"
#include "solidpress.h"
#include "solidinteg.h"
#include "solidshape.h"

template <sp_grad::index_type iNumNodes>
class PressureFromNodes {
public:
     PressureFromNodes()
          :rgNodes{nullptr} {
     }

     template <typename T>
     inline void
     GetNodalPressure(sp_grad::SpColVector<T, iNumNodes>& p,
                      doublereal dCoef,
                      sp_grad::SpFunctionCall func) const;


     void SetNode(sp_grad::index_type i, const ScalarNodeAd* pNode) {
          ASSERT(i >= 0);
          ASSERT(i < iNumNodes);

          rgNodes[i] = pNode;
     }

     const ScalarNodeAd* pGetNode(sp_grad::index_type i) const {
          ASSERT(i >= 0);
          ASSERT(i < iNumNodes);

          return rgNodes[i];
     }

     static inline constexpr int
     GetNumConnectedNodes() { return iNumNodes; }

     inline void
     GetConnectedNodes(std::vector<const Node*>& connectedNodes) const {
          for (const ScalarNodeAd* pNode: rgNodes) {
               connectedNodes.push_back(pNode);
          }
     }

     void PrintLogFile(std::ostream& of) const {
          for (const ScalarNodeAd* pNode: rgNodes) {
               of << ' ' << pNode->GetLabel();
          }
     }
private:
     std::array<const ScalarNodeAd*, iNumNodes> rgNodes;
};

template <sp_grad::index_type iNumDrives>
class PressureFromDrives {
public:
     PressureFromDrives() {}
     PressureFromDrives(PressureFromDrives&& oPressureTmp)
          :rgDrives(std::move(oPressureTmp.rgDrives)) {
     }

     template <typename T>
     inline void
     GetNodalPressure(sp_grad::SpColVector<T, iNumDrives>& p,
                      doublereal dCoef,
                      sp_grad::SpFunctionCall func) const;

     void SetDrive(sp_grad::index_type i, std::unique_ptr<DriveCaller>&& pDrive) {
          ASSERT(i >= 0);
          ASSERT(i < iNumDrives);

          rgDrives[i] = std::move(pDrive);
     }

     static inline constexpr int
     GetNumConnectedNodes() { return 0; }

     static inline void
     GetConnectedNodes(std::vector<const Node*>&) {
     }

     void PrintLogFile(std::ostream& of) const {
          for (const auto& pDrive: rgDrives) {
               of << ' ' << pDrive->GetLabel();
          }
     }
private:
     std::array<std::unique_ptr<DriveCaller>, iNumDrives> rgDrives;
};

template <typename ElementType, typename CollocationType, typename PressureSource>
class PressureLoad: public PressureLoadElem {
public:
     static constexpr sp_grad::index_type iNumNodes = ElementType::iNumNodes;
     static constexpr sp_grad::index_type iNumEvalPoints = CollocationType::iNumEvalPoints;
     static constexpr sp_grad::index_type iNumDof = iNumNodes * 3;

     PressureLoad(unsigned uLabel,
                  const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                  PressureSource&& oPressureTmp,
                  flag fOut);
     virtual ~PressureLoad();

     virtual void Output(OutputHandler& OH) const override;

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

     template <typename T>
     inline void
     InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                   const sp_grad::SpGradientVectorHandler<T>& XCurr,
                   enum sp_grad::SpFunctionCall func);

     virtual VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
                   const VectorHandler& XCurr) override;

     virtual SubVectorHandler&
     InitialAssRes(SubVectorHandler& WorkVec,
                   const VectorHandler& XCurr) override;

     virtual void
     InitialWorkSpaceDim(integer* piNumRows,
                         integer* piNumCols) const override;

     virtual int
     GetNumConnectedNodes() const override;

     virtual void
     GetConnectedNodes(std::vector<const Node*>& connectedNodes) const override;

protected:
     template <typename T>
     inline void
     AssPressureLoad(sp_grad::SpColVector<T, iNumNodes * 3>& f,
                     doublereal dCoef,
                     enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssVector(sp_grad::SpGradientAssVec<T>& WorkVec,
               sp_grad::SpColVector<T, iNumDof>& R,
               integer (StructDispNode::*pfnGetFirstIndex)(void) const);

     template <typename T>
     inline void
     GetNodalPosition(sp_grad::SpColVector<T, 3 * iNumNodes>& x,
                      doublereal dCoef,
                      sp_grad::SpFunctionCall func) const;

     inline void
     UpdateTotalForce(const sp_grad::SpColVector<doublereal, iNumDof>& R);

     inline void
     UpdateTotalForce(const sp_grad::SpColVector<sp_grad::SpGradient, iNumDof>& R) {}

     inline void
     UpdateTotalForce(const sp_grad::SpColVector<sp_grad::GpGradProd, iNumDof>& R) {}

     struct CollocData {
          static constexpr sp_grad::index_type iNumNodes = PressureLoad::iNumNodes;

          sp_grad::SpColVectorA<doublereal, iNumNodes> HA;
          sp_grad::SpMatrixA<doublereal, 3, 3 * iNumNodes> Hf, dHf_dr, dHf_ds;
     };

     std::array<const StructDispNodeAd*, iNumNodes> rgNodes;
     PressureSource oPressure;
     std::array<CollocData, iNumEvalPoints> rgCollocData;
     Vec3 Ftot;
};

template <sp_grad::index_type iNumDrives>
template <typename T>
inline void
PressureFromDrives<iNumDrives>::GetNodalPressure(sp_grad::SpColVector<T, iNumDrives>& p,
                                                 doublereal dCoef,
                                                 sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     for (index_type j = 1; j <= iNumDrives; ++j) {
          SpGradientTraits<T>::ResizeReset(p(j), rgDrives[j - 1]->dGet(), 0);
     }
}

template <sp_grad::index_type iNumNodes>
template <typename T>
inline void
PressureFromNodes<iNumNodes>::GetNodalPressure(sp_grad::SpColVector<T, iNumNodes>& p,
                                               doublereal dCoef,
                                               sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          rgNodes[j - 1]->GetX(p(j), dCoef, func);
     }
}

PressureLoadElem::PressureLoadElem(unsigned uLabel,
                                   flag fOut)
     :Elem(uLabel, fOut), InitialAssemblyElem(uLabel, fOut)
{
}

PressureLoadElem::~PressureLoadElem()
{
}

Elem::Type PressureLoadElem::GetElemType() const
{
     return Elem::PRESSURE_LOAD;
}

void
PressureLoadElem::SetValue(DataManager *pDM,
                           VectorHandler& X, VectorHandler& XP,
                           SimulationEntity::Hints *ph)
{
}

std::ostream& PressureLoadElem::Restart(std::ostream& out) const
{
     out << "## pressure load element: Restart not implemented yet\n";

     return out;
}

unsigned int PressureLoadElem::iGetInitialNumDof() const
{
     return 0;
}

bool PressureLoadElem::bIsDeformable() const
{
     return false;
}

template <typename ElementType, typename CollocationType, typename PressureSource>
PressureLoad<ElementType, CollocationType, PressureSource>::PressureLoad(unsigned uLabel,
                                                                         const std::array<const StructDispNodeAd*, iNumNodes>& rgNodesTmp,
                                                                         PressureSource&& oPressureTmp,
                                                                         flag fOut)
     :PressureLoadElem(uLabel, fOut),
      Elem(uLabel, fOut),
      rgNodes(rgNodesTmp),
      oPressure(std::move(oPressureTmp)),
      Ftot(::Zero3)
{
     using namespace sp_grad;

     SpColVector<doublereal, 2> r(2, 0);
     SpMatrix<doublereal, iNumNodes, 2> hd(iNumNodes, 2, 0);

     for (index_type i = 0; i < iNumEvalPoints; ++i) {
          CollocationType::GetPosition(i, r);
          ElementType::ShapeFunction(r, rgCollocData[i].HA);
          ElementType::ShapeFunctionDeriv(r, hd);

          for (index_type k = 1; k <= iNumNodes; ++k) {
               for (index_type j = 1; j <= 3; ++j) {
                    rgCollocData[i].Hf(j, (k - 1) * 3 + j) = rgCollocData[i].HA(k);
                    rgCollocData[i].dHf_dr(j, (k - 1) * 3 + j) = hd(k, 1);
                    rgCollocData[i].dHf_ds(j, (k - 1) * 3 + j) = hd(k, 2);
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource>
PressureLoad<ElementType, CollocationType, PressureSource>::~PressureLoad()
{
}

template <typename ElementType, typename CollocationType, typename PressureSource>
void PressureLoad<ElementType, CollocationType, PressureSource>::Output(OutputHandler& OH) const
{
     using namespace sp_grad;

     if (bToBeOutput() && OH.UseText(OutputHandler::PRESSURE_LOADS)) {
          if (OH.UseText(OutputHandler::PRESSURE_LOADS)) {
               std::ostream& of = OH.PressureLoads();

               of << std::setw(8) << GetLabel() << ' ' << Ftot << '\n';
          }
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource>
void PressureLoad<ElementType, CollocationType, PressureSource>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = iNumDof;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename PressureSource>
template <typename T>
inline void
PressureLoad<ElementType, CollocationType, PressureSource>::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                   doublereal dCoef,
                                                                   const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                   const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                                                   enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpColVector<T, iNumNodes * 3> R(iNumDof, (iNumDof + iNumNodes) * iNumEvalPoints);

     AssPressureLoad(R, dCoef, func);

     AssVector(WorkVec, R, &StructDispNodeAd::iGetFirstMomentumIndex);
}

template <typename ElementType, typename CollocationType, typename PressureSource>
template <typename T>
inline void
PressureLoad<ElementType, CollocationType, PressureSource>::AssPressureLoad(sp_grad::SpColVector<T, iNumNodes * 3>& R,
                                                                            doublereal dCoef,
                                                                            enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     T p_i;

     SpColVector<T, 3 * iNumNodes> x(3 * iNumNodes, 1);
     SpColVector<T, iNumNodes> p(iNumNodes, 1);
     SpColVector<T, 3> n1(3, iNumDof), n2(3, iNumDof), n(3, iNumDof), F_i(3, iNumDof + iNumNodes);
     SpColVector<T, iNumDof> HfT_F_i(iNumDof, iNumDof + iNumNodes);

     GetNodalPosition(x, dCoef, func);
     oPressure.GetNodalPressure(p, dCoef, func);

     sp_grad::SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(x);
     oDofMap.GetDofStat(p);
     oDofMap.Reset();
     oDofMap.InsertDof(x);
     oDofMap.InsertDof(p);
     oDofMap.InsertDone();

     for (index_type i = 0; i < iNumEvalPoints; ++i) {
          const doublereal alpha = CollocationType::dGetWeight(i);

          p_i = Dot(rgCollocData[i].HA, p);
          n1.MapAssign(rgCollocData[i].dHf_dr * x, oDofMap);
          n2.MapAssign(rgCollocData[i].dHf_ds * x, oDofMap);
          n.MapAssign(Cross(n1, n2), oDofMap);
          F_i.MapAssign(n * (-alpha * p_i), oDofMap);
          HfT_F_i.MapAssign(Transpose(rgCollocData[i].Hf) * F_i, oDofMap);
          R += HfT_F_i;
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource>
template <typename T>
inline void
PressureLoad<ElementType, CollocationType, PressureSource>::AssVector(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                      sp_grad::SpColVector<T, iNumDof>& R,
                                                                      integer (StructDispNode::*pfnGetFirstIndex)(void) const)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumNodes; ++i) {
          const index_type iEqIndex = (rgNodes[i - 1]->*pfnGetFirstIndex)();

          for (index_type j = 1; j <= 3; ++j) {
               WorkVec.AddItem(iEqIndex + j, R((i - 1) * 3 + j));
          }
     }

     UpdateTotalForce(R);
}

template <typename ElementType, typename CollocationType, typename PressureSource>
SubVectorHandler&
PressureLoad<ElementType, CollocationType, PressureSource>::AssRes(SubVectorHandler& WorkVec,
                                                                   doublereal dCoef,
                                                                   const VectorHandler& XCurr,
                                                                   const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("PressureLoad::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename PressureSource>
VariableSubMatrixHandler&
PressureLoad<ElementType, CollocationType, PressureSource>::AssJac(VariableSubMatrixHandler& WorkMat,
                                                                   doublereal dCoef,
                                                                   const VectorHandler& XCurr,
                                                                   const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("PressureLoad::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename PressureSource>
void
PressureLoad<ElementType, CollocationType, PressureSource>::AssJac(VectorHandler& JacY,
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

template <typename ElementType, typename CollocationType, typename PressureSource>
template <typename T>
void
PressureLoad<ElementType, CollocationType, PressureSource>::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                          const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                          enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpColVector<T, iNumNodes * 3> R(iNumDof, (iNumDof + iNumNodes) * iNumEvalPoints);

     AssPressureLoad(R, 1., func);

     AssVector(WorkVec, R, &StructDispNodeAd::iGetFirstPositionIndex);
}

template <typename ElementType, typename CollocationType, typename PressureSource>
VariableSubMatrixHandler&
PressureLoad<ElementType, CollocationType, PressureSource>::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                                                          const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::InitialAssJac(this,
                                                 WorkMat.SetSparseGradient(),
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_JAC);

     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename PressureSource>
SubVectorHandler&
PressureLoad<ElementType, CollocationType, PressureSource>::InitialAssRes(SubVectorHandler& WorkVec,
                                                                          const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                 WorkVec,
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename PressureSource>
void
PressureLoad<ElementType, CollocationType, PressureSource>::InitialWorkSpaceDim(integer* piNumRows,
                                                                                integer* piNumCols) const
{
     *piNumRows = iNumDof;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename PressureSource>
int
PressureLoad<ElementType, CollocationType, PressureSource>::GetNumConnectedNodes() const
{
     return iNumNodes + oPressure.GetNumConnectedNodes();
}

template <typename ElementType, typename CollocationType, typename PressureSource>
void
PressureLoad<ElementType, CollocationType, PressureSource>::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
     connectedNodes.reserve(GetNumConnectedNodes());
     connectedNodes.clear();

     for (const Node* pNode:rgNodes) {
          connectedNodes.push_back(pNode);
     }

     oPressure.GetConnectedNodes(connectedNodes);
}

template <typename ElementType, typename CollocationType, typename PressureSource>
template <typename T>
inline void
PressureLoad<ElementType, CollocationType, PressureSource>::GetNodalPosition(sp_grad::SpColVector<T, 3 * iNumNodes>& x,
                                                                             doublereal dCoef,
                                                                             sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     SpColVector<T, 3> Xj(3, 1);

     for (index_type j = 1; j <= iNumNodes; ++j) {
          rgNodes[j - 1]->GetXCurr(Xj, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               x(i + (j - 1) * 3) = std::move(Xj(i));
          }
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource>
inline void
PressureLoad<ElementType, CollocationType, PressureSource>::UpdateTotalForce(const sp_grad::SpColVector<doublereal, iNumDof>& R)
{
     using namespace sp_grad;

     Ftot = ::Zero3;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          for (index_type i = 1; i <= 3; ++i) {
               Ftot(i) += R((j - 1) * 3 + i);
          }
     }
}

template <typename ElementType, typename CollocationType>
PressureLoadElem*
ReadPressureLoad(DataManager* const pDM, MBDynParser& HP, const unsigned int uLabel)
{
     DEBUGCOUTFNAME("ReadPressureLoad");

     using namespace sp_grad;

     constexpr index_type iNumNodes = ElementType::iNumNodes;

     typedef PressureLoad<ElementType, CollocationType, PressureFromNodes<iNumNodes>> PressureLoadFromNodes;
     typedef PressureLoad<ElementType, CollocationType, PressureFromDrives<iNumNodes>> PressureLoadFromDrives;

     std::array<const StructDispNodeAd*, iNumNodes> rgNodes;

     enum { PRESSURE_FROM_NODES,
            PRESSURE_FROM_DRIVES,
            PRESSURE_UNKNOWN
     } ePressureSource = PRESSURE_UNKNOWN;

     PressureFromNodes<iNumNodes> oPressureFromNodes;
     PressureFromDrives<iNumNodes> oPressureFromDrives;

     for (index_type i = 0; i < iNumNodes; ++i) {
          rgNodes[i] = pDM->ReadNode<const StructDispNodeAd, Node::STRUCTURAL>(HP);
     }

     if (HP.IsKeyWord("from" "nodes")) {
          ePressureSource = PRESSURE_FROM_NODES;

          for (index_type i = 0; i < iNumNodes; ++i) {
               oPressureFromNodes.SetNode(i, pDM->ReadNode<const ScalarNodeAd, Node::ABSTRACT>(HP));
          }
     } else if (HP.IsKeyWord("from" "drives")) {
          ePressureSource = PRESSURE_FROM_DRIVES;

          std::unique_ptr<DriveCaller> pDrive;

          for (index_type i = 0; i < iNumNodes; ++i) {
               pDrive.reset(HP.GetDriveCaller());

               oPressureFromDrives.SetDrive(i, std::move(pDrive));
          }
     } else {
          silent_cerr("keyword \"from nodes\" or \"from drives\" expected at line " << HP.GetLineData() << "\n");
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     if (HP.IsArg()) {
          silent_cerr("semicolon expected "
                      "at line " << HP.GetLineData() << std::endl);
          throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     const flag fOut = pDM->fReadOutput(HP, Elem::PRESSURE_LOAD);

     std::ostream& out = pDM->GetLogFile();

     out << ElementType::ElementName() << ": " << uLabel;

     for (sp_grad::index_type i = 0; i < iNumNodes; ++i) {
          out << ' ' << rgNodes[i]->GetLabel();
     }

     PressureLoadElem* pElem = nullptr;

     switch (ePressureSource) {
     case PRESSURE_FROM_NODES:
          oPressureFromNodes.PrintLogFile(out);

          SAFENEWWITHCONSTRUCTOR(pElem,
                                 PressureLoadFromNodes,
                                 PressureLoadFromNodes(uLabel, rgNodes, std::move(oPressureFromNodes), fOut));
          break;
     case PRESSURE_FROM_DRIVES:
          oPressureFromDrives.PrintLogFile(out);

          SAFENEWWITHCONSTRUCTOR(pElem,
                                 PressureLoadFromDrives,
                                 PressureLoadFromDrives(uLabel, rgNodes, std::move(oPressureFromDrives), fOut));
          break;
     default:
          ASSERT(0);
     }

     out << '\n';

     return pElem;
}

template PressureLoadElem* ReadPressureLoad<Quadrangle4, Gauss2x2>(DataManager*, MBDynParser&, unsigned int);
template PressureLoadElem* ReadPressureLoad<Quadrangle8, Gauss3x3>(DataManager*, MBDynParser&, unsigned int);
template PressureLoadElem* ReadPressureLoad<Triangle6h, CollocTria6h>(DataManager*, MBDynParser&, unsigned int);
