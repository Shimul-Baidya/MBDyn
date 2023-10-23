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

template <sp_grad::index_type iNumDrives>
class SurfLoadFromDrives {
public:
     SurfLoadFromDrives() {}
     SurfLoadFromDrives(SurfLoadFromDrives&& oLoadTmp)
          :rgDrives(std::move(oLoadTmp.rgDrives)) {
     }

     template <typename T>
     inline void
     GetNodalLoad(sp_grad::SpMatrix<T, 3, iNumDrives>& F,
                  doublereal dCoef,
                  sp_grad::SpFunctionCall func) const {
          using namespace sp_grad;

          for (index_type j = 1; j <= iNumDrives; ++j) {
               const Vec3 Fj = rgDrives[j - 1]->Get();

               for (index_type i = 1; i <= 3; ++i) {
                    F(i, j) = Fj(i);
               }
          }
     }

     void SetDrive(sp_grad::index_type i, std::unique_ptr<TplDriveCaller<Vec3>>&& pDrive) {
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
     }
private:
     std::array<std::unique_ptr<TplDriveCaller<Vec3>>, iNumDrives> rgDrives;
};

template <typename ElementType, typename CollocationType, typename PressureSource>
class SurfaceLoad: public SurfaceLoadElem {
public:
     static constexpr sp_grad::index_type iNumNodes = ElementType::iNumNodes;
     static constexpr sp_grad::index_type iNumEvalPoints = CollocationType::iNumEvalPoints;
     static constexpr sp_grad::index_type iNumDof = iNumNodes * 3;

     SurfaceLoad(unsigned uLabel,
                 const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                 PressureSource&& oPressureTmp,
                 flag fOut);
     virtual ~SurfaceLoad();

     virtual int
     GetNumConnectedNodes() const override;

     virtual void
     GetConnectedNodes(std::vector<const Node*>& connectedNodes) const override;

protected:
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

     std::array<const StructDispNodeAd*, iNumNodes> rgNodes;
     PressureSource oPressure;

     struct CollocData {
          static constexpr sp_grad::index_type iNumNodes = SurfaceLoad::iNumNodes;

          sp_grad::SpColVectorA<doublereal, iNumNodes> HA;
          sp_grad::SpMatrixA<doublereal, 3, 3 * iNumNodes> Hf, dHf_dr, dHf_ds;
     };

     template <typename CollocDataType>
     static inline void InitCollocData(std::array<CollocDataType, iNumEvalPoints>& rgCollocData);
};

template <typename ElementType, typename CollocationType, typename PressureSource>
class PressureLoad: public SurfaceLoad<ElementType, CollocationType, PressureSource> {
     typedef SurfaceLoad<ElementType, CollocationType, PressureSource> BaseType;
public:
     using BaseType::iNumNodes;
     using BaseType::iNumEvalPoints;
     using BaseType::iNumDof;

     PressureLoad(unsigned uLabel,
                  const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                  PressureSource&& oPressureTmp,
                  flag fOut);
     virtual ~PressureLoad();

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

protected:
     template <typename T>
     inline void
     AssPressureLoad(sp_grad::SpColVector<T, iNumNodes * 3>& f,
                     doublereal dCoef,
                     enum sp_grad::SpFunctionCall func);

     std::array<typename BaseType::CollocData, iNumEvalPoints> rgCollocData;
};

enum class SurfaceTractionType {
     RELATIVE,
     ABSOLUTE
};

template <SurfaceTractionType eType>
struct SurfaceTractionHelper;

template <>
struct SurfaceTractionHelper<SurfaceTractionType::RELATIVE> {
     template <typename T, sp_grad::index_type iNumCols, typename ElementType>
     static inline void
     AssSurfaceTraction(ElementType& oElem,
                        sp_grad::SpColVector<T, iNumCols>& f,
                        doublereal dCoef,
                        enum sp_grad::SpFunctionCall func) {
          oElem.AssSurfaceTractionRel(f, dCoef, func);
     }

     template <typename ElementType, size_t uNumCols>
     static inline void InitCollocData(ElementType& oElem, const std::array<Mat3x3, uNumCols>& Rf) {
          oElem.InitCollocDataRel(Rf);
     }
};

template <>
struct SurfaceTractionHelper<SurfaceTractionType::ABSOLUTE> {
     template <typename T, sp_grad::index_type iNumCols, typename ElementType>
     static inline void
     AssSurfaceTraction(ElementType& oElem,
                        sp_grad::SpColVector<T, iNumCols>& f,
                        doublereal dCoef,
                        enum sp_grad::SpFunctionCall func) {
          oElem.AssSurfaceTractionAbs(f, dCoef, func);
     }

     template <typename ElementType, size_t uNumCols>
     static inline void InitCollocData(ElementType& oElem, const std::array<Mat3x3, uNumCols>& Rf) {
          oElem.InitCollocDataAbs(Rf);
     }
};

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
class SurfaceTraction: public SurfaceLoad<ElementType, CollocationType, PressureSource> {
     typedef SurfaceLoad<ElementType, CollocationType, PressureSource> BaseType;
public:
     using BaseType::iNumNodes;
     using BaseType::iNumEvalPoints;
     using BaseType::iNumDof;

     SurfaceTraction(unsigned uLabel,
                     const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                     PressureSource&& oPressureTmp,
                     const std::array<Mat3x3, iNumEvalPoints>& Rf,
                     flag fOut);
     virtual ~SurfaceTraction();

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

     template <typename T>
     inline void
     AssSurfaceTractionRel(sp_grad::SpColVector<T, iNumNodes * 3>& f,
                           doublereal dCoef,
                           enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssSurfaceTractionAbs(sp_grad::SpColVector<T, iNumNodes * 3>& f,
                           doublereal dCoef,
                           enum sp_grad::SpFunctionCall func);

     inline void InitCollocDataRel(const std::array<Mat3x3, iNumEvalPoints>& Rf);
     inline void InitCollocDataAbs(const std::array<Mat3x3, iNumEvalPoints>& Rf);
private:
     struct CollocData: public BaseType::CollocData {
          sp_grad::SpMatrixA<doublereal, 3, 3> Rrel;
          doublereal dA;
     };

     std::array<CollocData, iNumEvalPoints> rgCollocData;
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

SurfaceLoadElem::SurfaceLoadElem(unsigned uLabel,
                                 flag fOut)
     :Elem(uLabel, fOut),
      InitialAssemblyElem(uLabel, fOut),
      Ftot(::Zero3)
{
}

SurfaceLoadElem::~SurfaceLoadElem()
{
}

Elem::Type SurfaceLoadElem::GetElemType() const
{
     return Elem::SURFACE_LOAD;
}

void
SurfaceLoadElem::SetValue(DataManager *pDM,
                          VectorHandler& X, VectorHandler& XP,
                          SimulationEntity::Hints *ph)
{
}

std::ostream& SurfaceLoadElem::Restart(std::ostream& out) const
{
     out << "## pressure load element: Restart not implemented yet\n";

     return out;
}

unsigned int SurfaceLoadElem::iGetInitialNumDof() const
{
     return 0;
}

bool SurfaceLoadElem::bIsDeformable() const
{
     return false;
}

void SurfaceLoadElem::Output(OutputHandler& OH) const
{
     using namespace sp_grad;

     if (bToBeOutput() && OH.UseText(OutputHandler::SURFACE_LOADS)) {
          if (OH.UseText(OutputHandler::SURFACE_LOADS)) {
               std::ostream& of = OH.SurfaceLoads();

               of << std::setw(8) << GetLabel() << ' ' << Ftot << '\n';
          }
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource>
SurfaceLoad<ElementType, CollocationType, PressureSource>::SurfaceLoad(unsigned uLabel,
                                                                       const std::array<const StructDispNodeAd*, iNumNodes>& rgNodesTmp,
                                                                       PressureSource&& oPressureTmp,
                                                                       flag fOut)
     :Elem(uLabel, fOut),
      SurfaceLoadElem(uLabel, fOut),
      rgNodes(rgNodesTmp),
      oPressure(std::move(oPressureTmp))
{
}

template <typename ElementType, typename CollocationType, typename PressureSource>
SurfaceLoad<ElementType, CollocationType, PressureSource>::~SurfaceLoad()
{
}

template <typename ElementType, typename CollocationType, typename PressureSource>
template <typename T>
inline void
SurfaceLoad<ElementType, CollocationType, PressureSource>::AssVector(sp_grad::SpGradientAssVec<T>& WorkVec,
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
int
SurfaceLoad<ElementType, CollocationType, PressureSource>::GetNumConnectedNodes() const
{
     return iNumNodes + oPressure.GetNumConnectedNodes();
}

template <typename ElementType, typename CollocationType, typename PressureSource>
void
SurfaceLoad<ElementType, CollocationType, PressureSource>::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
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
SurfaceLoad<ElementType, CollocationType, PressureSource>::GetNodalPosition(sp_grad::SpColVector<T, 3 * iNumNodes>& x,
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
SurfaceLoad<ElementType, CollocationType, PressureSource>::UpdateTotalForce(const sp_grad::SpColVector<doublereal, iNumDof>& R)
{
     using namespace sp_grad;

     Ftot = ::Zero3;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          for (index_type i = 1; i <= 3; ++i) {
               Ftot(i) += R((j - 1) * 3 + i);
          }
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource>
template <typename CollocDataType>
inline void
SurfaceLoad<ElementType, CollocationType, PressureSource>::InitCollocData(std::array<CollocDataType, iNumEvalPoints>& rgCollocData)
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
PressureLoad<ElementType, CollocationType, PressureSource>::PressureLoad(unsigned uLabel,
                                                                         const std::array<const StructDispNodeAd*, iNumNodes>& rgNodesTmp,
                                                                         PressureSource&& oPressureTmp,
                                                                         flag fOut)
:Elem(uLabel, fOut),
 BaseType(uLabel, rgNodesTmp, std::move(oPressureTmp), fOut)
{
     BaseType::InitCollocData(rgCollocData);
}

template <typename ElementType, typename CollocationType, typename PressureSource>
PressureLoad<ElementType, CollocationType, PressureSource>::~PressureLoad()
{
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

     constexpr index_type iCapacity = iNumDof + PressureSource::GetNumConnectedNodes();

     SpColVector<T, iNumNodes * 3> R(iNumDof, iCapacity);

     AssPressureLoad(R, dCoef, func);

     ASSERT(R.iGetMaxSize() <= iCapacity);

     this->AssVector(WorkVec, R, &StructDispNodeAd::iGetFirstMomentumIndex);
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

     this->GetNodalPosition(x, dCoef, func);
     this->oPressure.GetNodalPressure(p, dCoef, func);

     sp_grad::SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(x);
     oDofMap.GetDofStat(p);
     oDofMap.Reset();
     oDofMap.InsertDof(x);
     oDofMap.InsertDof(p);
     oDofMap.InsertDone();

     for (auto& g: R) {
          oDofMap.InitDofMap(g);
     }

     for (index_type i = 0; i < iNumEvalPoints; ++i) {
          const doublereal alpha = CollocationType::dGetWeight(i);

          p_i = Dot(rgCollocData[i].HA, p);
          n1.MapAssign(rgCollocData[i].dHf_dr * x, oDofMap);
          n2.MapAssign(rgCollocData[i].dHf_ds * x, oDofMap);
          n.MapAssign(Cross(n1, n2), oDofMap);
          F_i.MapAssign(n * (-alpha * p_i), oDofMap);

          for (index_type k = 1; k <= iNumDof; ++k) {
               for (index_type j = 1; j <= 3; ++j) {
                    oDofMap.Add(R(k), rgCollocData[i].Hf(j, k) * F_i(j));
               }
          }
     }
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

     this->AssVector(WorkVec, R, &StructDispNodeAd::iGetFirstPositionIndex);
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

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::SurfaceTraction(unsigned uLabel,
                                                                                      const std::array<const StructDispNodeAd*, iNumNodes>& rgNodesTmp,
                                                                                      PressureSource&& oPressureTmp,
                                                                                      const std::array<Mat3x3, iNumEvalPoints>& Rf,
                                                                                      flag fOut)
:Elem(uLabel, fOut),
 BaseType(uLabel, rgNodesTmp, std::move(oPressureTmp), fOut)
{
     using namespace sp_grad;

     BaseType::InitCollocData(rgCollocData);

     SurfaceTractionHelper<eType>::InitCollocData(*this, Rf);
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
void SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::InitCollocDataRel(const std::array<Mat3x3, iNumEvalPoints>& Rf)
{
     using namespace sp_grad;

     SpColVector<doublereal, 3 * iNumNodes> x(3 * iNumNodes, 0);

     this->GetNodalPosition(x, 1., SpFunctionCall::REGULAR_RES);

     for (index_type i = 0; i < iNumEvalPoints; ++i) {
          SpMatrix<doublereal, 3, 3> Relem(3, 3, 0);
          SpColVector<doublereal, 3> e1 = rgCollocData[i].dHf_dr * x;
          SpColVector<doublereal, 3> e2 = rgCollocData[i].dHf_ds * x;
          SpColVector<doublereal, 3> e3 = Cross(e1, e2);

          e2 = Cross(e3, e1);

          const doublereal Norm_e1 = Norm(e1);
          const doublereal Norm_e2 = Norm(e2);
          const doublereal Norm_e3 = Norm(e3);

          if (Norm_e1 == 0. || Norm_e2 == 0. || Norm_e3 == 0.) {
               silent_cerr("traction(" << this->GetLabel() << "): orientation matrix is singular ("
                           << Norm_e1 << "," << Norm_e2 << ", " << Norm_e3 << ")\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          for (index_type j = 1; j <= 3; ++j) {
               Relem(j, 1) = e1(j) / Norm_e1;
               Relem(j, 2) = e2(j) / Norm_e2;
               Relem(j, 3) = e3(j) / Norm_e3;
          }

          rgCollocData[i].dA = Norm_e3;
          rgCollocData[i].Rrel = Transpose(Relem) * Rf[i];
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
void SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::InitCollocDataAbs(const std::array<Mat3x3, iNumEvalPoints>& Rf)
{
     using namespace sp_grad;

     SpColVector<doublereal, 3 * iNumNodes> x(3 * iNumNodes, 0);

     this->GetNodalPosition(x, 1., SpFunctionCall::REGULAR_RES);

     for (index_type i = 0; i < iNumEvalPoints; ++i) {
          SpColVector<doublereal, 3> e1 = rgCollocData[i].dHf_dr * x;
          SpColVector<doublereal, 3> e2 = rgCollocData[i].dHf_ds * x;
          SpColVector<doublereal, 3> e3 = Cross(e1, e2);

          const doublereal Norm_e3 = Norm(e3);

          if (Norm_e3 == 0.) {
               silent_cerr("traction(" << this->GetLabel() << "): orientation matrix is singular ("
                           << Norm_e3 << ")\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          rgCollocData[i].dA = Norm_e3;
          rgCollocData[i].Rrel = Rf[i];
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::~SurfaceTraction()
{
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
void SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = iNumDof;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
template <typename T>
inline void
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                             doublereal dCoef,
                                                                             const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                             const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                                                             enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     constexpr index_type iCapacity = iNumDof * (eType == SurfaceTractionType::RELATIVE);

     SpColVector<T, iNumNodes * 3> R(iNumDof, iCapacity);

     SurfaceTractionHelper<eType>::AssSurfaceTraction(*this, R, dCoef, func);

     ASSERT(R.iGetMaxSize() <= iCapacity);

     this->AssVector(WorkVec, R, &StructDispNodeAd::iGetFirstMomentumIndex);
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
template <typename T>
inline void
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::AssSurfaceTractionRel(sp_grad::SpColVector<T, iNumNodes * 3>& R,
                                                                                            doublereal dCoef,
                                                                                            enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     T Norm_e1, Norm_e2, dA;
     SpColVector<doublereal, 3> fl_i(3, 0), frel_i(3, 0);
     SpMatrix<T, 3, 3> Relem_i(3, 3, iNumDof);
     SpColVector<T, 3 * iNumNodes> x(3 * iNumNodes, 1);
     SpMatrix<doublereal, 3, iNumNodes> fl_n(3, iNumNodes, 0);
     SpColVector<T, 3> e1(3, iNumDof), e2(3, iNumDof), e3(3, iNumDof), Fg_i(3, iNumDof + iNumNodes);

     this->GetNodalPosition(x, dCoef, func);
     this->oPressure.GetNodalLoad(fl_n, dCoef, func);

     sp_grad::SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(x);
     oDofMap.Reset();
     oDofMap.InsertDof(x);
     oDofMap.InsertDone();

     for (auto& g: R) {
          oDofMap.InitDofMap(g);
     }

     for (index_type i = 0; i < iNumEvalPoints; ++i) {
          const doublereal alpha = CollocationType::dGetWeight(i);

          for (index_type j = 1; j <= 3; ++j) {
               fl_i(j) = Dot(rgCollocData[i].HA, Transpose(fl_n.GetRow(j))) * alpha;
          }

          frel_i = rgCollocData[i].Rrel * fl_i;

          e1.MapAssign(rgCollocData[i].dHf_dr * x, oDofMap);
          e2.MapAssign(rgCollocData[i].dHf_ds * x, oDofMap);

          e3.MapAssign(Cross(e1, e2), oDofMap);
          e2.MapAssign(Cross(e3, e1), oDofMap);

          oDofMap.MapAssign(Norm_e1, Norm(e1));
          oDofMap.MapAssign(Norm_e2, Norm(e2));
          oDofMap.MapAssign(dA, Norm(e3));

          if (Norm_e1 == 0. || Norm_e2 == 0. || dA == 0.) {
               silent_cerr("traction(" << this->GetLabel() << "): orientation matrix is singular\n");
               throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
          }

          for (index_type k = 1; k <= 3; ++k) {
               oDofMap.MapAssign(Relem_i(k, 1), e1(k) * dA / Norm_e1);
               oDofMap.MapAssign(Relem_i(k, 2), e2(k) * dA / Norm_e2);
               oDofMap.MapAssign(Relem_i(k, 3), e3(k));
          }

          Fg_i.MapAssign(Relem_i * frel_i, oDofMap);

          for (index_type k = 1; k <= iNumDof; ++k) {
               for (index_type j = 1; j <= 3; ++j) {
                    oDofMap.Add(R(k), rgCollocData[i].Hf(j, k) * Fg_i(j));
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
template <typename T>
inline void
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::AssSurfaceTractionAbs(sp_grad::SpColVector<T, iNumNodes * 3>& R,
                                                                                            doublereal dCoef,
                                                                                            enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpColVector<doublereal, 3> fl_i(3, 0);
     SpMatrix<doublereal, 3, iNumNodes> fl_n(3, iNumNodes, 0);
     SpColVector<doublereal, 3> Fg_i(3, iNumDof + iNumNodes);

     this->oPressure.GetNodalLoad(fl_n, dCoef, func);

     for (index_type i = 0; i < iNumEvalPoints; ++i) {
          const doublereal alpha = CollocationType::dGetWeight(i);

          for (index_type j = 1; j <= 3; ++j) {
               fl_i(j) = Dot(rgCollocData[i].HA, Transpose(fl_n.GetRow(j))) * alpha;
          }

          Fg_i = (rgCollocData[i].Rrel * fl_i) * rgCollocData[i].dA;

          for (index_type k = 1; k <= iNumDof; ++k) {
               for (index_type j = 1; j <= 3; ++j) {
                    R(k) += rgCollocData[i].Hf(j, k) * Fg_i(j);
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
SubVectorHandler&
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::AssRes(SubVectorHandler& WorkVec,
                                                                             doublereal dCoef,
                                                                             const VectorHandler& XCurr,
                                                                             const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SurfaceTraction::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
VariableSubMatrixHandler&
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::AssJac(VariableSubMatrixHandler& WorkMat,
                                                                             doublereal dCoef,
                                                                             const VectorHandler& XCurr,
                                                                             const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SurfaceTraction::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
void
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::AssJac(VectorHandler& JacY,
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

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
template <typename T>
void
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                    const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                                    enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpColVector<T, iNumNodes * 3> R(iNumDof, (iNumDof + iNumNodes) * iNumEvalPoints);

     SurfaceTractionHelper<eType>::AssSurfaceTraction(*this, R, 1., func);

     this->AssVector(WorkVec, R, &StructDispNodeAd::iGetFirstPositionIndex);
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
VariableSubMatrixHandler&
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                                                                    const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::InitialAssJac(this,
                                                 WorkMat.SetSparseGradient(),
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_JAC);

     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
SubVectorHandler&
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::InitialAssRes(SubVectorHandler& WorkVec,
                                                                                    const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                 WorkVec,
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename PressureSource, SurfaceTractionType eType>
void
SurfaceTraction<ElementType, CollocationType, PressureSource, eType>::InitialWorkSpaceDim(integer* piNumRows,
                                                                                          integer* piNumCols) const
{
     *piNumRows = iNumDof;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType>
SurfaceLoadElem*
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

     const flag fOut = pDM->fReadOutput(HP, Elem::SURFACE_LOAD);

     std::ostream& out = pDM->GetLogFile();

     out << ElementType::ElementName() << ": " << uLabel;

     for (sp_grad::index_type i = 0; i < iNumNodes; ++i) {
          out << ' ' << rgNodes[i]->GetLabel();
     }

     SurfaceLoadElem* pElem = nullptr;

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

     if (HP.IsArg()) {
          silent_cerr("semicolon expected "
                      "at line " << HP.GetLineData() << std::endl);
          throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     return pElem;
}

template SurfaceLoadElem* ReadPressureLoad<Quadrangle4, Gauss2x2>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadPressureLoad<Quadrangle8, Gauss3x3>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadPressureLoad<Quadrangle9, Gauss3x3>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadPressureLoad<Quadrangle8r, Gauss3x3>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadPressureLoad<Triangle6h, CollocTria6h>(DataManager*, MBDynParser&, unsigned int);

template <typename ElementType, typename CollocationType>
SurfaceLoadElem*
ReadTractionLoad(DataManager* const pDM, MBDynParser& HP, const unsigned int uLabel)
{
     DEBUGCOUTFNAME("ReadTractionLoad");

     using namespace sp_grad;

     constexpr index_type iNumNodes = ElementType::iNumNodes;
     constexpr index_type iNumEvalPoints = CollocationType::iNumEvalPoints;

     typedef SurfaceTraction<ElementType, CollocationType, SurfLoadFromDrives<iNumNodes>, SurfaceTractionType::RELATIVE> SurfaceTractionRel;
     typedef SurfaceTraction<ElementType, CollocationType, SurfLoadFromDrives<iNumNodes>, SurfaceTractionType::ABSOLUTE> SurfaceTractionAbs;

     const SurfaceTractionType eTractionType = HP.IsKeyWord("absolute") ? SurfaceTractionType::ABSOLUTE : SurfaceTractionType::RELATIVE;

     std::array<const StructDispNodeAd*, iNumNodes> rgNodes;

     SurfLoadFromDrives<iNumNodes> oSurfLoadFromDrives;

     for (index_type i = 0; i < iNumNodes; ++i) {
          rgNodes[i] = pDM->ReadNode<const StructDispNodeAd, Node::STRUCTURAL>(HP);
     }

     if (HP.IsKeyWord("from" "drives")) {
          std::unique_ptr<TplDriveCaller<Vec3>> pDrive;

          for (index_type i = 0; i < iNumNodes; ++i) {
               pDrive.reset(ReadDC3D(pDM, HP));

               oSurfLoadFromDrives.SetDrive(i, std::move(pDrive));
          }
     } else {
          silent_cerr("keyword \"from drives\" expected at line " << HP.GetLineData() << "\n");
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     ReferenceFrame oGlobalRefFrame;
     std::array<Mat3x3, iNumEvalPoints> Rf;
     const bool bReadOrientation = HP.IsKeyWord("orientation");

     for (Mat3x3& Rf_i: Rf) {
          Rf_i = bReadOrientation ? HP.GetRotAbs(oGlobalRefFrame) : Eye3;
     }

     const flag fOut = pDM->fReadOutput(HP, Elem::SURFACE_LOAD);

     std::ostream& out = pDM->GetLogFile();

     if (HP.IsArg()) {
          silent_cerr("semicolon expected "
                      "at line " << HP.GetLineData() << std::endl);
          throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     out << ElementType::ElementName() << ": " << uLabel;

     for (sp_grad::index_type i = 0; i < iNumNodes; ++i) {
          out << ' ' << rgNodes[i]->GetLabel();
     }

     SurfaceLoadElem* pElem = nullptr;

     oSurfLoadFromDrives.PrintLogFile(out);

     switch (eTractionType) {
     case SurfaceTractionType::RELATIVE:
          SAFENEWWITHCONSTRUCTOR(pElem,
                                 SurfaceTractionRel,
                                 SurfaceTractionRel(uLabel, rgNodes, std::move(oSurfLoadFromDrives), Rf, fOut));
          break;
     case SurfaceTractionType::ABSOLUTE:
          SAFENEWWITHCONSTRUCTOR(pElem,
                                 SurfaceTractionAbs,
                                 SurfaceTractionAbs(uLabel, rgNodes, std::move(oSurfLoadFromDrives), Rf, fOut));
          break;
     }

     out << '\n';

     return pElem;
}

template SurfaceLoadElem* ReadTractionLoad<Quadrangle4, Gauss2x2>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadTractionLoad<Quadrangle8, Gauss3x3>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadTractionLoad<Quadrangle9, Gauss3x3>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadTractionLoad<Quadrangle8r, Gauss3x3>(DataManager*, MBDynParser&, unsigned int);
template SurfaceLoadElem* ReadTractionLoad<Triangle6h, CollocTria6h>(DataManager*, MBDynParser&, unsigned int);
