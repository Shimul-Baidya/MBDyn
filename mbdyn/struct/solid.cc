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
#include <type_traits>

#include <ac/lapack.h>
#include "nodead.h"
#include "strnodead.h"
#include "sp_matvecass.h"
#include "constltp.h"

#include "solid.h"
#include "solidinteg.h"
#include "solidshape.h"
#include "dataman.h"

template <typename CSL>
class SolidConstLaw {
public:
     typedef std::unique_ptr<CSL> ConstLawPtr;
     static constexpr sp_grad::index_type iDim = CSL::iDim;

     void SetConstLaw(ConstLawPtr&& pConstLawTmp) {
          pConstLaw = std::move(pConstLawTmp);
          ASSERT(pConstLaw != nullptr);
     }

     ConstLawType* operator->() {
          ASSERT(pConstLaw != nullptr);

          return pConstLaw.get();
     }

     void AfterConvergence() {
          ASSERT(pConstLaw != nullptr);

          pConstLaw->AfterConvergence(pConstLaw->GetEpsilon(), pConstLaw->GetEpsilonPrime());
     }
protected:
     ConstLawPtr pConstLaw;
};

class SolidConstLawElastic: public SolidConstLaw<ConstitutiveLaw6D> {
public:
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::ELASTIC;

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, sp_grad::SpColVector<T, 6>& sigma, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          using namespace sp_grad;

          const SpColVector<T, 6> eps{G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1)};

          pConstLaw->Update(eps, sigma, oDofMap);
     }
};

class SolidConstLawElasticIncompr: public SolidConstLaw<ConstitutiveLaw7D> {
public:
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::ELASTICINCOMPR;

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, sp_grad::SpColVector<T, 7>& sigma, const T& ptilde, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          using namespace sp_grad;

          const SpColVector<T, 7> eps{G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1), ptilde};

          pConstLaw->Update(eps, sigma, oDofMap);
     }
};

class SolidConstLawViscoelastic: public SolidConstLaw<ConstitutiveLaw6D> {
public:
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::VISCOELASTIC;

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpMatrix<T, 3, 3>& GP, sp_grad::SpColVector<T, 6>& sigma, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          using namespace sp_grad;

          const SpColVector<T, 6> eps{G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1)};
          const SpColVector<T, 6> epsP{GP(1, 1), GP(2, 2), GP(3, 3), 2. * GP(1, 2), 2. * GP(2, 3), 2. * GP(3, 1)};

          pConstLaw->Update(eps, epsP, sigma, oDofMap);
     }
};

template <bool bEnableUPC, sp_grad::index_type iNumNodes>
class IncomprSolidElemStatic;

template <sp_grad::index_type iNumNodesPressureTpl>
class IncomprSolidElemStatic<false, iNumNodesPressureTpl> {
public:
     static constexpr sp_grad::index_type iNumNodesPressure = iNumNodesPressureTpl;

     explicit IncomprSolidElemStatic(const std::array<const ScalarNodeAd*, iNumNodesPressure>&) {
     }

     template <typename T>
     inline void
     GetNodalPressuresUPC(sp_grad::SpColVector<T, iNumNodesPressure>&,
                       doublereal,
                       sp_grad::SpFunctionCall) const = delete;

     sp_grad::index_type iGetEquIndexPressureUPC(sp_grad::index_type iNode) const = delete;
};

template <sp_grad::index_type iNumNodesPressureTpl>
class IncomprSolidElemStatic<true, iNumNodesPressureTpl> {
public:
     static constexpr sp_grad::index_type iNumNodesPressure = iNumNodesPressureTpl;

     static_assert(iNumNodesPressure > 0);

     explicit IncomprSolidElemStatic(const std::array<const ScalarNodeAd*, iNumNodesPressure>& rgNodesTmp)
          :rgNodesPressure(rgNodesTmp) {
     }

     template <typename T>
     inline void
     GetNodalPressuresUPC(sp_grad::SpColVector<T, iNumNodesPressure>& p,
                       doublereal dCoef,
                       sp_grad::SpFunctionCall func) const {
          using namespace sp_grad;

          for (index_type i = 1; i <= iNumNodesPressure; ++i) {
               rgNodesPressure[i - 1]->GetX(p(i), dCoef, func);
          }
     }

     sp_grad::index_type iGetEquIndexPressureUPC(sp_grad::index_type iNodeIndex) const {
          return rgNodesPressure[iNodeIndex - 1]->iGetFirstRowIndex() + 1;
     }
private:
     const std::array<const ScalarNodeAd*, iNumNodesPressure> rgNodesPressure;
};

template <bool bEnableUPC, typename ElementType>
struct CollocDataPressureUPC {
     static inline void
     Init(const sp_grad::SpColVector<doublereal, 3>& r, const sp_grad::SpMatrix<doublereal, 3, ElementType::iNumNodes>& x0, const SolidElem* pElem) {}
};

template <typename ElementType>
struct CollocDataPressureUPC<true, ElementType> {
     typedef typename ElementType::ElemTypePressureUPC ElemTypeUPC;

     inline void
     Init(const sp_grad::SpColVector<doublereal, 3>& r, const sp_grad::SpMatrix<doublereal, 3, ElementType::iNumNodes>& x0, const SolidElem* pElem) {
          ElemTypeUPC::ShapeFunction(r, g);
     }

     sp_grad::SpColVectorA<doublereal, ElementType::iNumNodesPressure> g;
};

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType = StructDispNodeAd>
class SolidElemStatic: public SolidElem, public IncomprSolidElemStatic<SolidCSLType::eConstLawType == ConstLawType::ELASTICINCOMPR, ElementType::iNumNodesPressure> {
public:
     typedef IncomprSolidElemStatic<SolidCSLType::eConstLawType == ConstLawType::ELASTICINCOMPR, ElementType::iNumNodesPressure> IncomprSolidElemType;

     static constexpr ConstLawType::Type eConstLawType = SolidCSLType::eConstLawType;
     static constexpr sp_grad::index_type iNumNodes = ElementType::iNumNodes;
     static constexpr sp_grad::index_type iNumNodesPressure = IncomprSolidElemType::iNumNodesPressure;
     static constexpr sp_grad::index_type iNumNodesExtrap = ElementType::iNumNodesExtrap;
     static constexpr sp_grad::index_type iNumEvalPointsStiffness = CollocationType::iNumEvalPointsStiffness;
     static constexpr sp_grad::index_type iNumEvalPointsMassLumped = CollocationType::iNumEvalPointsMassLumped;
     static constexpr sp_grad::index_type iNumDof = iNumNodes * 3;
     static constexpr sp_grad::index_type iNumDofPressure = iNumNodesPressure;

     static_assert(ElementType::eElemFlags == SolidElemFlags::DISPLACEMENT || ElementType::eElemFlags == SolidElemFlags::DISPLACEMENT_PRESSURE);
     static_assert(iNumNodesPressure <= iNumNodes);
     static_assert(eConstLawType == ConstLawType::ELASTICINCOMPR ? iNumNodesPressure > 0 : iNumNodesPressure == 0);

     SolidElemStatic(unsigned uLabel,
                     const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                     const std::array<const ScalarNodeAd*, iNumNodesPressure>& rgNodesPressure,
                     const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                     std::array<typename SolidCSLType::ConstLawPtr, iNumEvalPointsStiffness>&& rgMaterialData,
                     const RigidBodyKinematics* pRBK,
                     flag fOut);
     virtual ~SolidElemStatic();

     virtual void Output(OutputHandler& OH) const override;

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                   doublereal dCoef,
                   enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResElasticIncompr(sp_grad::SpGradientAssVec<T>& WorkVec,
                          doublereal dCoef,
                          enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                        doublereal dCoef,
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

     virtual void
     AfterConvergence(const VectorHandler& X,
                      const VectorHandler& XP) override;

     virtual int
     GetNumConnectedNodes() const override;

     virtual void
     GetConnectedNodes(std::vector<const Node*>& connectedNodes) const override;

     virtual doublereal dGetM() const override;
     virtual Vec3 GetS_int() const override;
     virtual Mat3x3 GetJ_int() const override;
protected:
     template <typename T>
     inline void
     GetNodalPositions(sp_grad::SpMatrix<T, 3, iNumNodes>& x,
                       doublereal dCoef,
                       sp_grad::SpFunctionCall func) const;
     template <typename T>
     inline void
     GetNodalDeformations(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                          doublereal dCoef,
                          sp_grad::SpFunctionCall func) const;

     template <typename T>
     inline void
     GetNodalVelocities(sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                        doublereal dCoef,
                        sp_grad::SpFunctionCall func) const;

     inline void
     ComputeStressStrain(sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& epsilon,
                         sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& tau) const;

     template <typename T>
     inline void
     AssStiffnessVecElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                            sp_grad::SpColVector<T, iNumDof>& R,
                            doublereal dCoef,
                            sp_grad::SpFunctionCall func,
                            const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     template <typename T>
     inline void
     AssStiffnessVecElasticIncompr(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                   const sp_grad::SpColVector<T, iNumDofPressure>& p,
                                   sp_grad::SpColVector<T, iNumDof>& FU,
                                   sp_grad::SpColVector<T, iNumDofPressure>& FP,
                                   doublereal dCoef,
                                   sp_grad::SpFunctionCall func,
                                   const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);


     template <typename T>
     inline void
     AssStiffnessVecViscoElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                 const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                 sp_grad::SpColVector<T, iNumDof>& R,
                                 doublereal dCoef,
                                 sp_grad::SpFunctionCall func,
                                 const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     template <typename T>
     inline void
     AssGravityLoadVec(sp_grad::SpColVector<T, iNumDof>& R,
                       doublereal dCoef,
                       sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssVector(sp_grad::SpGradientAssVec<T>& WorkVec,
               sp_grad::SpColVector<T, iNumDof>& R,
               integer (StructDispNode::*pfnGetFirstIndex)(void) const);


     template <typename T>
     inline void
     AssVectorIncompr(sp_grad::SpGradientAssVec<T>& WorkVec,
                      sp_grad::SpColVector<T, iNumDofPressure>& FP);

     template <sp_grad::index_type iNumComp>
     inline void
     GaussToNodal(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& taun,
                  const sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, iNumComp>& taug) const;

     template <typename T>
     inline void
     AssInertiaVecRBK(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                      sp_grad::SpColVector<T, iNumDof>& R,
                      const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     struct CollocData: CollocDataPressureUPC<eConstLawType == ConstLawType::ELASTICINCOMPR, ElementType> {
          typedef CollocDataPressureUPC<eConstLawType == ConstLawType::ELASTICINCOMPR, ElementType> CollocDataUPC;

          static constexpr sp_grad::index_type iNumNodes = SolidElemStatic::iNumNodes;
          static constexpr sp_grad::index_type iNumDof = SolidElemStatic::iNumDof;

          inline void
          Init(sp_grad::index_type iColloc,
               const sp_grad::SpMatrix<doublereal, 3, iNumNodes>& x0,
               typename SolidCSLType::ConstLawPtr&& pConstLawTmp,
               const SolidElemStatic* pElem);

          template <typename T>
          inline void
          AddInternalForceVector(const sp_grad::SpMatrix<T, 3, 3>& F,
                                 const sp_grad::SpColVector<T, SolidCSLType::iDim>& sigma,
                                 const doublereal alpha,
                                 sp_grad::SpColVector<T, iNumDof>& R,
                                 const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) const;

          template <typename T>
          inline void
          ComputeStressElastic(const sp_grad::SpMatrix<T, 3, 3>& G,
                               const sp_grad::SpMatrix<T, 3, 3>& F,
                               sp_grad::SpColVector<T, SolidCSLType::iDim>& sigma,
                               const sp_grad::SpGradExpDofMapHelper<T>& oDofMap,
                               const SolidElemStatic* pElem) {
               oConstLaw.Update(G, sigma, oDofMap);
               UpdateStressStrain(G, sigma, F, pElem);
          }

          template <typename T>
          inline void
          ComputeStressElasticIncompr(const sp_grad::SpMatrix<T, 3, 3>& G,
                                      const sp_grad::SpMatrix<T, 3, 3>& F,
                                      sp_grad::SpColVector<T, SolidCSLType::iDim>& sigma,
                                      const T& ptilde,
                                      const sp_grad::SpGradExpDofMapHelper<T>& oDofMap,
                                      const SolidElemStatic* pElem) {
               oConstLaw.Update(G, sigma, ptilde, oDofMap);
               UpdateStressStrain(G, sigma, F, pElem);
          }

          template <typename T>
          inline void
          ComputeStressViscoElastic(const sp_grad::SpMatrix<T, 3, 3>& G,
                                    const sp_grad::SpMatrix<T, 3, 3>& GP,
                                    const sp_grad::SpMatrix<T, 3, 3>& F,
                                    sp_grad::SpColVector<T, SolidCSLType::iDim>& sigma,
                                    const sp_grad::SpGradExpDofMapHelper<T>& oDofMap,
                                    const SolidElemStatic* pElem) {
               oConstLaw.Update(G, GP, sigma, oDofMap);
               UpdateStressStrain(G, sigma, F, pElem);
          }

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& G_tmp,
                             const sp_grad::SpColVector<doublereal, SolidCSLType::iDim>& sigma_tmp,
                             const sp_grad::SpMatrix<doublereal, 3, 3>& F_tmp,
                             const SolidElemStatic* pElem);

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& G,
                             const sp_grad::SpColVector<sp_grad::SpGradient, SolidCSLType::iDim>& sigma,
                             const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& F,
                             const SolidElemStatic* pElem) {
          }

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& G,
                             const sp_grad::SpColVector<sp_grad::GpGradProd, SolidCSLType::iDim>& sigma,
                             const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& F,
                             const SolidElemStatic* pElem) {
          }

          void AfterConvergence() {
               oConstLaw.AfterConvergence();
          }

          SolidCSLType oConstLaw;
          sp_grad::SpColVectorA<doublereal, iNumNodes> h;
          sp_grad::SpMatrixA<doublereal, iNumNodes, 3> h0d;
          doublereal detJ;
          sp_grad::SpMatrixA<doublereal, 3, 3> G;
          sp_grad::SpMatrixA<doublereal, 3, 3> F;
          sp_grad::SpColVectorA<doublereal, 6> sigma;
     };

     sp_grad::SpMatrixA<doublereal, 3, iNumNodes> x0;
     const std::array<const StructNodeType*, iNumNodes> rgNodes;
     const sp_grad::SpColVectorA<doublereal, iNumNodes> rhon;
     std::array<CollocData, iNumEvalPointsStiffness> rgCollocData;
     const RigidBodyKinematics* const pRBK;
};

enum class MassMatrixType {
     CONSISTENT,
     LUMPED
};

template <MassMatrixType eMassMatrix>
struct MassMatrixHelper;

template <>
struct MassMatrixHelper<MassMatrixType::CONSISTENT> {
     template <sp_grad::index_type iNumDof, sp_grad::index_type iNumNodes>
     struct MassMatrix {
          typedef sp_grad::SpMatrixA<doublereal, iNumNodes, iNumNodes> Type;
     };

     template <typename ElemType>
     static void AssMassMatrix(ElemType& oElem, sp_grad::SpMatrix<doublereal, ElemType::iNumNodes, ElemType::iNumNodes>& Mcon) {
          oElem.AssMassMatrixConsistent(Mcon);
     }

     template <typename ElemType>
     static void AddInertia(ElemType& oElem, const sp_grad::SpMatrix<doublereal, ElemType::iNumNodes, ElemType::iNumNodes>& Mcon) {
          DEBUGCERR("warning: AddInertia() not implemented for consistent mass matrix\n");
     }

     template <typename ElemType, typename T>
     static void
     AssInertiaVec(ElemType& oElem,
                   const sp_grad::SpMatrix<doublereal, ElemType::iNumNodes, ElemType::iNumNodes>& Mcon,
                   const sp_grad::SpMatrix<T, 3, ElemType::iNumNodes>& uP,
                   sp_grad::SpColVector<T, ElemType::iNumDof>& R,
                   const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          oElem.AssInertiaVecConsistent(Mcon, uP, R, oDofMap);
     }
};

template <>
struct MassMatrixHelper<MassMatrixType::LUMPED> {
     template <sp_grad::index_type iNumDof, sp_grad::index_type iNumNodes>
     struct MassMatrix {
          typedef sp_grad::SpColVectorA<doublereal, iNumNodes> Type;
     };

     template <typename ElemType>
     static void AssMassMatrix(ElemType& oElem, sp_grad::SpColVector<doublereal, ElemType::iNumNodes>& diagM) {
          oElem.AssMassMatrixLumped(diagM);
     }

     template <typename ElemType>
     static void AddInertia(ElemType& oElem, sp_grad::SpColVector<doublereal, ElemType::iNumNodes>& diagM) {
          oElem.AddInertia(diagM);
     }

     template <typename ElemType, typename T>
     static void
     AssInertiaVec(ElemType& oElem,
                   const sp_grad::SpColVector<doublereal, ElemType::iNumNodes>& diagM,
                   const sp_grad::SpMatrix<T, 3, ElemType::iNumNodes>& uP,
                   sp_grad::SpColVector<T, ElemType::iNumDof>& R,
                   const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          oElem.AssInertiaVecLumped(diagM, uP, R, oDofMap);
     }
};

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
class SolidElemDynamic: public SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructDispNodeAd> {
public:
     typedef SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructDispNodeAd> SolidElemStaticType;
     using SolidElemStaticType::eConstLawType;
     using SolidElemStaticType::iNumNodes;
     using SolidElemStaticType::iNumNodesPressure;
     using SolidElemStaticType::iNumEvalPointsStiffness;
     using SolidElemStaticType::iNumEvalPointsMassLumped;
     using SolidElemStaticType::iNumDof;
     using SolidElemStaticType::iNumDofPressure;

     SolidElemDynamic(unsigned uLabel,
                      const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                      const std::array<const ScalarNodeAd*, iNumNodesPressure>& rgNodesPressure,
                      const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                      std::array<typename SolidCSLType::ConstLawPtr, iNumEvalPointsStiffness>&& rgMaterialData,
                      const RigidBodyKinematics* pRBK,
                      flag fOut);
     virtual ~SolidElemDynamic();

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     template <typename T>
     inline void
     AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                   doublereal dCoef,
                   enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResElasticIncompr(sp_grad::SpGradientAssVec<T>& WorkVec,
                          doublereal dCoef,
                          enum sp_grad::SpFunctionCall func);

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

     virtual Vec3 GetB_int() const override;
     virtual Vec3 GetG_int() const override;

     template <typename T>
     inline void
     AssInertiaVecConsistent(const sp_grad::SpMatrix<doublereal, iNumNodes, iNumNodes>& Mcon,
                             const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                             sp_grad::SpColVector<T, iNumDof>& R,
                             const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     template <typename T>
     inline void
     AssInertiaVecLumped(const sp_grad::SpColVector<doublereal, iNumNodes>& Mlumped,
                         const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                         sp_grad::SpColVector<T, iNumDof>& R,
                         const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     inline void
     AssMassMatrixConsistent(sp_grad::SpMatrix<doublereal, iNumNodes, iNumNodes>& Mcon);

     inline void
     AssMassMatrixLumped(sp_grad::SpColVector<doublereal, iNumNodes>& Mlumped);

     inline void
     AddInertia(const sp_grad::SpColVector<doublereal, iNumNodes>& diagM);

     virtual void
     SetValue(DataManager *pDM,
              VectorHandler& X, VectorHandler& XP,
              SimulationEntity::Hints *ph) override;

     virtual doublereal dGetE() const override;

private:
     typename MassMatrixHelper<eMassMatrix>::template MassMatrix<iNumDof, iNumNodes>::Type M;
};

SolidElem::SolidElem(unsigned uLabel,
                     flag fOut)
     :Elem{uLabel, fOut},
      ElemGravityOwner{uLabel, fOut},
      InitialAssemblyElem(uLabel, fOut)
{
}

SolidElem::~SolidElem()
{
}

void
SolidElem::SetValue(DataManager *pDM,
                    VectorHandler& X, VectorHandler& XP,
                    SimulationEntity::Hints *ph)
{
}

std::ostream&
SolidElem::Restart(std::ostream& out) const
{
     out << "## solid element: Restart not implemented yet\n";

     return out;
}

Elem::Type SolidElem::GetElemType() const
{
     return SOLID;
}

bool SolidElem::bIsDeformable() const
{
     return true;
}

unsigned int SolidElem::iGetInitialNumDof() const
{
     return 0u;
}

unsigned int SolidElem::iGetNumPrivData() const
{
     return 1;
}

unsigned int SolidElem::iGetPrivDataIdx(const char *s) const
{
     if (strcmp(s, "E") == 0) {
          return 1u;
     }

     throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal SolidElem::dGetPrivData(unsigned int i) const
{
     switch (i) {
     case 1u:
          return dGetE();

     default:
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }
}

doublereal SolidElem::dGetE() const
{
     return 0.; // Used for all static elements
}

template <ConstLawType::Type SolidCSLType>
struct SolidElemCSLHelper;

template <>
struct SolidElemCSLHelper<ConstLawType::ELASTIC> {
     template <typename SolidElementType, typename T>
     static inline void AssRes(SolidElementType* pElem,
                               sp_grad::SpGradientAssVec<T>& WorkVec,
                               doublereal dCoef,
                               enum sp_grad::SpFunctionCall func) {
          pElem->AssResElastic(WorkVec, dCoef, func);
     }
};

template <>
struct SolidElemCSLHelper<ConstLawType::ELASTICINCOMPR> {
     template <typename SolidElementType, typename T>
     static inline void AssRes(SolidElementType* pElem,
                               sp_grad::SpGradientAssVec<T>& WorkVec,
                               doublereal dCoef,
                               enum sp_grad::SpFunctionCall func) {
          pElem->AssResElasticIncompr(WorkVec, dCoef, func);
     }
};

template <>
struct SolidElemCSLHelper<ConstLawType::VISCOELASTIC> {
     template <typename SolidElementType, typename T>
     static inline void AssRes(SolidElementType* pElem,
                               sp_grad::SpGradientAssVec<T>& WorkVec,
                               doublereal dCoef,
                               enum sp_grad::SpFunctionCall func) {
          pElem->AssResViscoElastic(WorkVec, dCoef, func);
     }
};

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumNodes;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumEvalPointsStiffness;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumEvalPointsMassLumped;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumDof;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::SolidElemStatic(unsigned uLabel,
                                                                                             const std::array<const StructDispNodeAd*, iNumNodes>& rgNodesDisp,
                                                                                             const std::array<const ScalarNodeAd*, iNumNodesPressure>& rgNodesPressure,
                                                                                             const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                                                                                             std::array<typename SolidCSLType::ConstLawPtr, iNumEvalPointsStiffness>&& rgMaterialData,
                                                                                             const RigidBodyKinematics* pRBK,
                                                                                             flag fOut)
: Elem{uLabel, fOut},
  SolidElem{uLabel, fOut},
  IncomprSolidElemType{rgNodesPressure},
  rgNodes{rgNodesDisp},
  rhon{rhon},
  pRBK{pRBK}
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumNodes; ++i) {
          //rgNodes[i - 1] = rgNodesDisp[i - 1];

          const Vec3& x0i = rgNodes[i - 1]->GetXCurr();

          for (index_type j = 1; j <= 3; ++j) {
               x0(j, i) = x0i(j);
          }
     }

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          rgCollocData[iColloc].Init(iColloc, x0, std::move(rgMaterialData[iColloc]), this);
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::~SolidElemStatic()
{
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::Output(OutputHandler& OH) const
{
     using namespace sp_grad;

     if (bToBeOutput() && OH.UseText(OutputHandler::SOLIDS)) {
          sp_grad::SpMatrixA<doublereal, iNumEvalPointsStiffness, 6> epsilone;
          sp_grad::SpMatrixA<doublereal, iNumEvalPointsStiffness, 6> taue;
          sp_grad::SpMatrixA<doublereal, iNumNodes, 6> epsilonn;
          sp_grad::SpMatrixA<doublereal, iNumNodes, 6> taun;

          ComputeStressStrain(epsilone, taue);
          GaussToNodal(epsilonn, epsilone);
          GaussToNodal(taun, taue);

          if (OH.UseText(OutputHandler::SOLIDS)) {
               std::ostream& of = OH.Solids();

               of << std::setw(8) << GetLabel();

               for (index_type j = 1; j <= iNumNodes; ++j) {
                    for (index_type i = 1; i <= 6; ++i) {
                         of << ' ' << epsilonn(j, i);
                    }

                    for (index_type i = 1; i <= 6; ++i) {
                         of << ' ' << taun(j, i);
                    }
               }

               of << '\n';
          }
     }
}



template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = iNumDof + iNumDofPressure;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                    doublereal dCoef,
                                                                                    const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                                    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                                                                    enum sp_grad::SpFunctionCall func)
{
     return SolidElemCSLHelper<eConstLawType>::AssRes(this, WorkVec, dCoef, func);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                           const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                                           enum sp_grad::SpFunctionCall func)
{
     if constexpr(eConstLawType == ConstLawType::ELASTICINCOMPR) {
          silent_cerr("solid(" << GetLabel() << "): "
                      "initial assembly with incompressible constitutive laws is not yet implemented!\n");
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS); // because of hydraulic nodes
     } else {
          return SolidElemCSLHelper<eConstLawType>::AssRes(this, WorkVec, 1., func);
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                           doublereal dCoef,
                                                                                           enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1);

     GetNodalDeformations(u, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDone();

     SpColVector<T, iNumDof> R(iNumDof, oDofMap);

     AssStiffnessVecElastic(u, R, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (pRBK) {
          AssInertiaVecRBK(u, R, oDofMap);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (pGravity) {
          AssGravityLoadVec(R, dCoef, func);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssResElasticIncompr(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                                  doublereal dCoef,
                                                                                                  enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1);
     SpColVector<T, iNumNodesPressure> p(iNumNodesPressure, 1);

     GetNodalDeformations(u, dCoef, func);
     this->GetNodalPressuresUPC(p, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(p);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(p);
     oDofMap.InsertDone();

     SpColVector<T, iNumDof> FU(iNumDof, oDofMap);
     SpColVector<T, iNumDofPressure> FP(iNumDofPressure, oDofMap);

     AssStiffnessVecElasticIncompr(u, p, FU, FP, dCoef, func, oDofMap);

     ASSERT(FU.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (pRBK) {
          AssInertiaVecRBK(u, FU, oDofMap);
     }

     ASSERT(FU.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (pGravity) {
          AssGravityLoadVec(FU, dCoef, func);
     }

     ASSERT(FU.iGetMaxSize() == oDofMap.iGetLocalSize());

     AssVector(WorkVec, FU, &StructDispNode::iGetFirstMomentumIndex);

     AssVectorIncompr(WorkVec, FP);

     ASSERT(FP.iGetMaxSize() == oDofMap.iGetLocalSize());
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                                doublereal dCoef,
                                                                                                enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1), uP(3, iNumNodes, 1);

     GetNodalDeformations(u, dCoef, func);
     GetNodalVelocities(uP, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(uP);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(uP);
     oDofMap.InsertDone();

     SpColVector<T, iNumDof> R(iNumDof, oDofMap);

     AssStiffnessVecViscoElastic(u, uP, R, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (pRBK) {
          AssInertiaVecRBK(u, R, oDofMap);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (pGravity) {
          AssGravityLoadVec(R, dCoef, func);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);
}

template <typename T>
inline void GreenLagrangeStrain(const sp_grad::SpMatrix<T, 3, 3>& F, sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     // G:0.5 * (transpose(F).F - ident(3));
     // exploit symmetry of the strain tensor

     oDofMap.MapAssign(G(1,1), 5.0E-1*(F(3,1)*F(3,1)+F(2,1)*F(2,1)+F(1,1)*F(1,1)-1));
     oDofMap.MapAssign(G(1,2), 5.0E-1*(F(3,1)*F(3,2)+F(2,1)*F(2,2)+F(1,1)*F(1,2)));
     oDofMap.MapAssign(G(1,3), 5.0E-1*(F(3,1)*F(3,3)+F(2,1)*F(2,3)+F(1,1)*F(1,3)));
     G(2,1) = G(1, 2);
     oDofMap.MapAssign(G(2,2), 5.0E-1*(F(3,2)*F(3,2)+F(2,2)*F(2,2)+F(1,2)*F(1,2)-1));
     oDofMap.MapAssign(G(2,3), 5.0E-1*(F(3,2)*F(3,3)+F(2,2)*F(2,3)+F(1,2)*F(1,3)));
     G(3,1) = G(1, 3);
     G(3,2) = G(2, 3);
     oDofMap.MapAssign(G(3,3), 5.0E-1*(F(3,3)*F(3,3)+F(2,3)*F(2,3)+F(1,3)*F(1,3)-1));
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssStiffnessVecElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                                    sp_grad::SpColVector<T, iNumDof>& R,
                                                                                                    const doublereal dCoef,
                                                                                                    const sp_grad::SpFunctionCall func,
                                                                                                    const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     SpMatrix<T, 3, 3> F(3, 3, oDofMap), G(3, 3, oDofMap);
     SpColVector<T, 6> sigma(6, oDofMap);

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          const doublereal alpha = CollocationType::dGetWeightStiffness(iColloc);
          const auto& h0d = rgCollocData[iColloc].h0d;

          F.MapAssign(u * h0d, oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               F(i, i) += 1.;
          }

          ASSERT(F.iGetMaxSize() <= iNumDof);

          GreenLagrangeStrain(F, G, oDofMap);

          ASSERT(G.iGetMaxSize() <= iNumDof);

          rgCollocData[iColloc].ComputeStressElastic(G, F, sigma, oDofMap, this);

          rgCollocData[iColloc].AddInternalForceVector(F, sigma, alpha, R, oDofMap);

          ASSERT(R.iGetMaxSize() <= iNumDof);
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssStiffnessVecElasticIncompr(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                                           const sp_grad::SpColVector<T, iNumDofPressure>& p,
                                                                                                           sp_grad::SpColVector<T, iNumDof>& FU,
                                                                                                           sp_grad::SpColVector<T, iNumDofPressure>& FP,
                                                                                                           doublereal dCoef,
                                                                                                           sp_grad::SpFunctionCall func,
                                                                                                           const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     SpMatrix<T, 3, 3> F(3, 3, oDofMap), G(3, 3, oDofMap);
     SpColVector<T, SolidCSLType::iDim> sigma(SolidCSLType::iDim, oDofMap);

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          const doublereal alpha = CollocationType::dGetWeightStiffness(iColloc);
          const auto& h0d = rgCollocData[iColloc].h0d;
          const doublereal detJg = rgCollocData[iColloc].detJ;

          F.MapAssign(u * h0d, oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               F(i, i) += 1.;
          }

          ASSERT(F.iGetMaxSize() <= iNumDof + iNumDofPressure);

          GreenLagrangeStrain(F, G, oDofMap);

          ASSERT(G.iGetMaxSize() <= iNumDof + iNumDofPressure);

          const SpColVector<doublereal, iNumDofPressure>& g = rgCollocData[iColloc].g;

          const T ptilde = Dot(g, p, oDofMap);

          rgCollocData[iColloc].ComputeStressElasticIncompr(G, F, sigma, ptilde, oDofMap, this);

          rgCollocData[iColloc].AddInternalForceVector(F, sigma, alpha, FU, oDofMap);

          ASSERT(FU.iGetMaxSize() <= iNumDof + iNumDofPressure);

          for (index_type i = 1; i <= iNumDofPressure; ++i) {
               oDofMap.Add(FP(i), -alpha * detJg * g(i) * sigma(SolidCSLType::iDim));
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssStiffnessVecViscoElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                                         const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                                                                                         sp_grad::SpColVector<T, iNumDof>& R,
                                                                                                         const doublereal dCoef,
                                                                                                         const sp_grad::SpFunctionCall func,
                                                                                                         const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     SpColVector<T, 6> sigma(6, oDofMap);
     SpMatrix<T, 3, 3> F(3, 3, oDofMap), FP(3, 3, oDofMap), C(3, 3, oDofMap), invC(3, 3, oDofMap), G(3, 3, oDofMap), FP_Tr_F(3, 3, oDofMap);
     SpMatrix<T, 3, 3> GP_detF(3, 3, oDofMap), GP_scaled(3, 3, oDofMap);

     T detC;

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          const auto& h0d = rgCollocData[iColloc].h0d;

          F.MapAssign(u * h0d, oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               F(i, i) += 1.;
          }

          ASSERT(F.iGetMaxSize() == oDofMap.iGetLocalSize());

          FP.MapAssign(uP * h0d, oDofMap);

          ASSERT(FP.iGetMaxSize() == oDofMap.iGetLocalSize());

          C.MapAssign(Transpose(F) * F, oDofMap);

          ASSERT(C.iGetMaxSize() == oDofMap.iGetLocalSize());

          InvSymm(C, invC, detC, oDofMap);

          ASSERT(invC.iGetMaxSize() == oDofMap.iGetLocalSize());

          G.MapAssign(0.5 * (C - Eye3), oDofMap);

          ASSERT(G.iGetMaxSize() == oDofMap.iGetLocalSize());

          FP_Tr_F.MapAssign(Transpose(FP) * F, oDofMap);

          ASSERT(FP_Tr_F.iGetMaxSize() == oDofMap.iGetLocalSize());

          GP_detF.MapAssign((0.5 * (FP_Tr_F + Transpose(FP_Tr_F)) * Det(F)), oDofMap);

          ASSERT(GP_detF.iGetMaxSize() == oDofMap.iGetLocalSize());

          GP_scaled.MapAssign(invC * GP_detF * invC, oDofMap); // Lars Kuebler 2005, equation 2.92, page 38

          ASSERT(GP_scaled.iGetMaxSize() == oDofMap.iGetLocalSize());

          rgCollocData[iColloc].ComputeStressViscoElastic(G, GP_scaled, F, sigma, oDofMap, this);

          ASSERT(sigma.iGetMaxSize() == oDofMap.iGetLocalSize());

          const doublereal alpha = CollocationType::dGetWeightStiffness(iColloc);

          rgCollocData[iColloc].AddInternalForceVector(F, sigma, alpha, R, oDofMap);

          ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssGravityLoadVec(sp_grad::SpColVector<T, iNumDof>& R,
                                                                                               doublereal dCoef,
                                                                                               sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 X, fg;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);

          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);

          X = Zero3;

          for (index_type i = 1; i <= iNumNodes; ++i) {
               const Vec3& Xi = rgNodes[i - 1]->GetXCurr();

               for (index_type j = 1; j <= 3; ++j) {
                    X(j) += Xi(j) * h(i);
               }
          }

          bGetGravity(X, fg); // FIXME: no contribution to the Jacobian but it would be required only for the CentralGravity

          const doublereal rho = Dot(h, rhon);

          fg *= rho * alpha * Det(J);

          for (index_type i = 1; i <= iNumNodes; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    R((i - 1) * 3 + j) += h(i) * fg(j);
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssVector(sp_grad::SpGradientAssVec<T>& WorkVec,
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
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssVectorIncompr(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                              sp_grad::SpColVector<T, iNumDofPressure>& FP)
{
     using namespace sp_grad;

     static_assert(iNumNodesPressure == iNumDofPressure);

     for (index_type i = 1; i <= iNumNodesPressure; ++i) {
          WorkVec.AddItem(this->iGetEquIndexPressureUPC(i), FP(i));
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SubVectorHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssRes(SubVectorHandler& WorkVec,
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

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
VariableSubMatrixHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssJac(VariableSubMatrixHandler& WorkMat,
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

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssJac(VectorHandler& JacY,
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

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
VariableSubMatrixHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                                                                           const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::InitialAssJac(this,
                                                 WorkMat.SetSparseGradient(),
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_JAC);

     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SubVectorHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialAssRes(SubVectorHandler& WorkVec,
                                                                                           const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                 WorkVec,
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialWorkSpaceDim(integer* piNumRows,
                                                                                                 integer* piNumCols) const
{
     *piNumRows = iNumDof + iNumDofPressure;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AfterConvergence(const VectorHandler& X,
                                                                                              const VectorHandler& XP)
{
     for (auto& oColloc: rgCollocData) {
          oColloc.AfterConvergence();
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
int
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetNumConnectedNodes() const
{
     return iNumNodes;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
     using namespace sp_grad;

     connectedNodes.resize(iNumNodes);

     for (index_type i = 0; i < iNumNodes; ++i) {
          connectedNodes[i] = rgNodes[i];
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
doublereal
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::dGetM() const
{
     using namespace sp_grad;

     doublereal dm = 0.;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, rhon); // interpolate from nodes to collocation points
          dm += rho * detJ * alpha;
     }

     return dm;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
Vec3
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetS_int() const
{
     using namespace sp_grad;

     Vec3 dS(::Zero3);

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     SpMatrixA<doublereal, 3, iNumNodes> x;

     GetNodalPositions(x, 1., SpFunctionCall::REGULAR_RES);

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);
          const SpColVector<doublereal, 3> X = x * h;
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, rhon); // interpolate from nodes to collocation points

          for (index_type j = 1; j <= 3; ++j) {
               dS(j) += rho * detJ * alpha * X(j);
          }
     }

     return dS;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
Mat3x3
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetJ_int() const
{
     using namespace sp_grad;

     Mat3x3 dJ(::Zero3x3);

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     SpMatrixA<doublereal, 3, iNumNodes> x;

     GetNodalPositions(x, 1., SpFunctionCall::REGULAR_RES);

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);
          const SpColVector<doublereal, 3> X = x * h;
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, rhon); // interpolate from nodes to collocation points

          dJ -= Mat3x3(MatCrossCross, Vec3(X.begin()), Vec3(X.begin()) * (rho * detJ * alpha));
     }

     return dJ;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetNodalPositions(sp_grad::SpMatrix<T, 3, iNumNodes>& x,
                                                                                               const doublereal dCoef,
                                                                                               const sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     SpColVectorA<T, 3, 1> Xj;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          rgNodes[j - 1]->GetXCurr(Xj, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               x(i, j) = std::move(Xj(i));
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetNodalDeformations(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                                  const doublereal dCoef,
                                                                                                  const sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     GetNodalPositions(u, dCoef, func);

     for (index_type j = 1; j <= iNumNodes; ++j) {
          for (index_type i = 1; i <= 3; ++i) {
               u(i, j) -= x0(i, j);
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetNodalVelocities(sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                                                                                const doublereal dCoef,
                                                                                                const sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     SpColVectorA<T, 3, 1> XPj;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          rgNodes[j - 1]->GetVCurr(XPj, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               uP(i, j) = std::move(XPj(i));
          }
     }
}

namespace {
     template <int N>
     struct log2int {
          static constexpr int value = (N / 2 > 0) ? (log2int<N / 2>::value + 1) : 0;
     };

     template <>
     struct log2int<1> {
          static constexpr int value = 0;
     };

     template <>
     struct log2int<0> {
          static constexpr int value = std::numeric_limits<int>::min();
     };

     static_assert(log2int<0>::value < 0, "logarithm unit test failed");
     static_assert(log2int<1>::value == 0, "logarithm unit test failed");
     static_assert(log2int<2>::value == 1, "logarithm unit test failed");
     static_assert(log2int<4>::value == 2, "logarithm unit test failed");
     static_assert(log2int<8>::value == 3, "logarithm unit test failed");
     static_assert(log2int<16>::value == 4, "logarithm unit test failed");
     static_assert(log2int<1024>::value == 10, "logarithm unit test failed");
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <sp_grad::index_type iNumComp>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GaussToNodal(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& taun,
                                                                                          const sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, iNumComp>& taug) const
{
#ifdef HAVE_DGELSD
     using namespace sp_grad;

     static_assert(iNumNodesExtrap <= iNumNodes, "invalid number of nodes");
     static_assert(iNumNodesExtrap <= iNumEvalPointsStiffness, "invalid number of nodes");

     SpMatrixA<doublereal, iNumEvalPointsStiffness, iNumNodesExtrap> H;
     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodesExtrap> h;

     for (index_type i = 1; i <= iNumEvalPointsStiffness; ++i) {
          CollocationType::GetPositionStiffness(i - 1, r);
          ElementType::ShapeFunctionExtrap(r, h);

          for (index_type j = 1; j <= iNumNodesExtrap; ++j) {
               H(i, j) = h(j); // FIXME: select only a subset of available nodes
          }
     }

     constexpr integer M = iNumEvalPointsStiffness;
     constexpr integer N = iNumNodesExtrap;
     constexpr integer MINMN = M < N ? M : N;
     constexpr integer NRHS = iNumComp;
     constexpr integer LDB = M;
     constexpr integer SMLSIZ = 25; // FIXME: this value should be determined by configure
     constexpr integer NLVLTMP = log2int<static_cast<integer>(MINMN / (SMLSIZ + 1.0) + 0.5)>::value + 1;
     constexpr integer NLVL = NLVLTMP > 0 ? NLVLTMP : 0;
     constexpr integer LIWORK = 3 * MINMN * NLVL + 11 * MINMN;
     constexpr integer LWORK = M >= N
          ? 12 * N + 2 * N * SMLSIZ + 8 * N * NLVL + N * NRHS + (SMLSIZ + 1) * (SMLSIZ + 1)
          : 12 * M + 2 * M * SMLSIZ + 8 * M * NLVL + M * NRHS + (SMLSIZ + 1) * (SMLSIZ + 1);
     doublereal RCOND = -1., WORK[LWORK], S[MINMN];
     integer RANK, IWORK[LIWORK], INFO;

     SpMatrixA<doublereal, LDB, NRHS> B;

     for (index_type j = 1; j <= NRHS; ++j) {
          for (index_type i = 1; i <= iNumEvalPointsStiffness; ++i) {
               B(i, j) = taug(i, j);
          }
     }

     __FC_DECL__(dgelsd)(&M, &N, &NRHS, H.begin(), &M, B.begin(), &LDB, S, &RCOND, &RANK, WORK, &LWORK, IWORK, &INFO);

     if (INFO != 0) {
          silent_cerr("solid(" << GetLabel() << "): dgelsd failed with status " << INFO << "\n");
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     ElementType::GaussToNodalInterp(taun, B);
#else
     silent_cerr("Output of solid element data is not available because LAPACK's dgelsd function was not found.\n"
                 "It may be disabled by means of a \"default output:\" statement within the control data section.\n");
     throw ErrNotAvailableYet(MBDYN_EXCEPT_ARGS);
#endif
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssInertiaVecRBK(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                              sp_grad::SpColVector<T, iNumDof>& R,
                                                                                              const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(pRBK != nullptr);

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     SpColVector<T, 3> frbk(3, oDofMap.iGetLocalSize());
     SpColVector<T, 3> s(3, iNumNodes);

     const Mat3x3 WxWx_WPx = Mat3x3(MatCrossCross, pRBK->GetW(), pRBK->GetW()) + Mat3x3(MatCross, pRBK->GetWP());
     const Vec3 XPP = pRBK->GetXPP();

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);

          const doublereal rho = Dot(h, rhon);
          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const doublereal dm = rho * alpha * Det(J);

          for (index_type i = 1; i <= 3; ++i) {
               SpGradientTraits<T>::ResizeReset(s(i), 0., iNumNodes);
          }

          for (index_type j = 1; j <= iNumNodes; ++j) {
               const doublereal hjdm = h(j) * dm;

               for (index_type i = 1; i <= 3; ++i) {
                    s(i) += (u(i, j) + x0(i, j)) * hjdm;
               }
          }

          for (index_type i = 1; i <= 3; ++i) {
               SpGradientTraits<T>::ResizeReset(frbk(i), XPP(i) * dm, oDofMap.iGetLocalSize());

               oDofMap.InitDofMap(frbk(i));

               for (index_type k = 1; k <= 3; ++k) {
                    oDofMap.Add(frbk(i), WxWx_WPx(i, k) * s(k));
               }
          }

          ASSERT(frbk.iGetMaxSize() <= iNumDof);

          for (index_type i = 1; i <= iNumNodes; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    oDofMap.Sub(R((i - 1) * 3 + j), h(i) * frbk(j));
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::ComputeStressStrain(sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& epsilon,
                                                                                                 sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& tau) const
{
     using namespace sp_grad;

     static constexpr index_type idxr[6] = {1, 2, 3, 1, 2, 3};
     static constexpr index_type idxc[6] = {1, 2, 3, 2, 3, 1};

     SpMatrixA<doublereal, 3, 3> S;

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
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
               epsilon(iColloc + 1, i) = sqrt(1. + 2. * G(i, i)) - 1.;
          }

          for (index_type i = 4; i <= 6; ++i) {
               epsilon(iColloc + 1, i) = 2. * G(idxr[i - 1], idxc[i - 1]) / ((1. + epsilon(iColloc + 1, idxr[i - 1])) * (1. + epsilon(iColloc + 1, idxc[i - 1])));
          }

          for (index_type i = 1; i <= 6; ++i) {
               tau(iColloc + 1, i) = taui(idxr[i - 1], idxc[i - 1]);
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::Init(const sp_grad::index_type iColloc,
                                                                                              const sp_grad::SpMatrix<doublereal, 3, iNumNodes>& x0,
                                                                                              typename SolidCSLType::ConstLawPtr&& pConstLawTmp,
                                                                                              const SolidElemStatic* const pElem)
{
     using namespace sp_grad;

     oConstLaw.SetConstLaw(std::move(pConstLawTmp));

     SpColVectorA<doublereal, 3> r;
     SpMatrixA<doublereal, ElementType::iNumNodes, 3> hd;

     CollocationType::GetPositionStiffness(iColloc, r);
     ElementType::ShapeFunction(r, h);
     ElementType::ShapeFunctionDeriv(r, hd);
     CollocDataUPC::Init(r, x0, pElem);

     const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);
     SpMatrixA<doublereal, 3, 3> invJ;

     Inv(J, invJ, detJ);

     if (detJ <= 0.) {
          silent_cerr("solid(" << pElem->GetLabel() << "): Jacobian is singular: det(J)=" << detJ << "\n");
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     h0d = hd * Transpose(invJ);

     F = ::Eye3; // Avoid singular matrix if "abort after: input;" is used
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::UpdateStressStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& G_tmp,
                                                                                                            const sp_grad::SpColVector<doublereal, SolidCSLType::iDim>& sigma_tmp,
                                                                                                            const sp_grad::SpMatrix<doublereal, 3, 3>& F_tmp,
                                                                                                            const SolidElemStatic* const pElem)
{
     using namespace sp_grad;

     G = G_tmp;
     sigma = SubColVector<1,1,6>(sigma_tmp);
     F = F_tmp;

     const doublereal detF = Det(F);

     if (detF <= 0.) {
          silent_cerr("solid(" << pElem->GetLabel() << "): element became excessively distorted: Det(F) = " << detF << "\n");
          throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::AddInternalForceVector(const sp_grad::SpMatrix<T, 3, 3>& F,
                                                                                                                const sp_grad::SpColVector<T, SolidCSLType::iDim>& sigma,
                                                                                                                const doublereal alpha,
                                                                                                                sp_grad::SpColVector<T, iNumDof>& R,
                                                                                                                const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) const
{
     using namespace sp_grad;

     const doublereal c1 = 0.5 * alpha * detJ;

     static constexpr index_type idxS[3][3] = {
          {1, 4, 6},
          {4, 2, 5},
          {6, 5, 3}
     };

     T a4, a5;

     for (index_type i = 1; i <= 3; ++i) {
          for (index_type j = i; j <= 3; ++j) { // exploit symmetry of stress and strain tensor
               const T& Sij = sigma(idxS[i - 1][j - 1]);
               const doublereal c2 = (i == j) ? c1 : 2. * c1; // because we are exploiting symmetry

               for (index_type k = 1; k <= 3; ++k) {
                    oDofMap.MapAssign(a4, c2 * F(k, j) * Sij);
                    oDofMap.MapAssign(a5, c2 * F(k, i) * Sij);

                    for (index_type l = 1; l <= iNumNodes; ++l) {
                         oDofMap.Sub(R((l - 1) * 3 + k), a4 * h0d(l, i) + a5 * h0d(l, j));
                    }
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::SolidElemDynamic(unsigned uLabel,
                                                                                            const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                                                                                            const std::array<const ScalarNodeAd*, iNumNodesPressure>& rgNodesPressure,
                                                                                            const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                                                                                            std::array<typename SolidCSLType::ConstLawPtr, iNumEvalPointsStiffness>&& rgMaterialData,
                                                                                            const RigidBodyKinematics* const pRBK,
                                                                                            flag fOut)
:Elem{uLabel, fOut},
 SolidElemStaticType{uLabel, rgNodes, rgNodesPressure, rhon, std::move(rgMaterialData), pRBK, fOut}
{
     MassMatrixHelper<eMassMatrix>::AssMassMatrix(*this, M);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::~SolidElemDynamic()
{
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
void SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssMassMatrixConsistent(sp_grad::SpMatrix<doublereal, iNumNodes, iNumNodes>& Mcon)
{
     using namespace sp_grad;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << this->GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, this->rhon); // interpolate from nodes to collocation points
          const doublereal dm = rho * detJ * alpha;

          for (index_type k = 1; k <= iNumNodes; ++k) {
               for (index_type j = 1; j <= iNumNodes; ++j) {
                    Mcon(j, k) += dm * h(j) * h(k);
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssMassMatrixLumped(sp_grad::SpColVector<doublereal, iNumNodes>& Mlumped)
{
     using namespace sp_grad;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     doublereal mtot = 0.;

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);

          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);
          const doublereal detJ = Det(J);
          const doublereal rho = Dot(h, this->rhon); // interpolate from nodes to collocation points
          const doublereal dm = rho * detJ * alpha;

          if (detJ <= 0.) {
               silent_cerr("solid(" << this->GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          mtot += dm;

          for (index_type k = 1; k <= iNumNodes; ++k) {
               Mlumped(k) += dm * std::pow(h(k), 2);
          }
     }

     doublereal mdiag = 0.;

     for (index_type k = 1; k <= iNumNodes; ++k) {
          mdiag += Mlumped(k);
     }

     for (index_type k = 1; k <= iNumNodes; ++k) {
          Mlumped(k) *= mtot / mdiag;
     }

     DEBUGCERR("Mlumped=" << Mlumped << "\n");
     DEBUGCERR("mtot=" << mtot << "\n");
     DEBUGCERR("mdiag=" << mdiag << "\n");
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AddInertia(const sp_grad::SpColVector<doublereal, iNumNodes>& diagM)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumNodes; ++i) {
          ASSERT(dynamic_cast<const DynamicStructDispNodeAd*>(this->rgNodes[i - 1]));

          static_cast<const DynamicStructDispNodeAd*>(this->rgNodes[i - 1])->AddInertia(diagM(i));
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::SetValue(DataManager *pDM,
                                                                                    VectorHandler& X, VectorHandler& XP,
                                                                                    SimulationEntity::Hints *ph)
{
     using namespace sp_grad;

     SpMatrix<doublereal, 3, iNumNodes> uP(3, iNumNodes, 1);

     this->GetNodalVelocities(uP, 1., SpFunctionCall::REGULAR_RES);

     SpGradExpDofMapHelper<doublereal> oDofMap;
     SpColVectorA<doublereal, iNumDof> beta;

     MassMatrixHelper<eMassMatrix>::AssInertiaVec(*this, M, uP, beta, oDofMap);

     for (index_type i = 1; i <= iNumNodes; ++i) {
          ASSERT(dynamic_cast<const DynamicStructDispNodeAd*>(this->rgNodes[i - 1]));

          const index_type iFirstIndex = this->rgNodes[i - 1]->iGetFirstMomentumIndex();

          for (index_type j = 1; j <= 3; ++j) {
               X.DecCoef(iFirstIndex + j, beta((i - 1) * 3 + j));
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
void SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 2 * iNumDof + iNumDofPressure;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
template <typename T>
inline void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                         doublereal dCoef,
                                                                                         enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1), uP(3, iNumNodes, 1);

     this->GetNodalDeformations(u, dCoef, func);
     this->GetNodalVelocities(uP, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(uP);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(uP);
     oDofMap.InsertDone();

     SpColVector<T, iNumDof> R(iNumDof, oDofMap);

     this->AssStiffnessVecElastic(u, R, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (this->pGravity) {
          this->AssGravityLoadVec(R, dCoef, func);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (this->pRBK) {
          this->AssInertiaVecRBK(u, R, oDofMap);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);

     MassMatrixHelper<eMassMatrix>::AssInertiaVec(*this, M, uP, R, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstPositionIndex);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
template <typename T>
inline void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssResElasticIncompr(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                                doublereal dCoef,
                                                                                                enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1), uP(3, iNumNodes, 1);
     SpColVector<T, iNumNodesPressure> p(iNumNodesPressure, 1);

     this->GetNodalDeformations(u, dCoef, func);
     this->GetNodalVelocities(uP, dCoef, func);
     this->GetNodalPressuresUPC(p, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(uP);
     oDofMap.GetDofStat(p);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(uP);
     oDofMap.InsertDof(p);
     oDofMap.InsertDone();

     SpColVector<T, iNumDof> R(iNumDof, oDofMap);
     SpColVector<T, iNumDofPressure> FP(iNumDofPressure, oDofMap);

     this->AssStiffnessVecElasticIncompr(u, p, R, FP, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (this->pGravity) {
          this->AssGravityLoadVec(R, dCoef, func);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (this->pRBK) {
          this->AssInertiaVecRBK(u, R, oDofMap);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);

     this->AssVectorIncompr(WorkVec, FP);

     ASSERT(FP.iGetMaxSize() == oDofMap.iGetLocalSize());

     MassMatrixHelper<eMassMatrix>::AssInertiaVec(*this, M, uP, R, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstPositionIndex);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
template <typename T>
inline void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                              doublereal dCoef,
                                                                                              enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1), uP(3, iNumNodes, 1);

     this->GetNodalDeformations(u, dCoef, func);
     this->GetNodalVelocities(uP, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(uP);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(uP);
     oDofMap.InsertDone();

     SpColVector<T, iNumDof> R(iNumDof, oDofMap);

     this->AssStiffnessVecViscoElastic(u, uP, R, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (this->pGravity) {
          this->AssGravityLoadVec(R, dCoef, func);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     if (this->pRBK) {
          this->AssInertiaVecRBK(u, R, oDofMap);
     }

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     MassMatrixHelper<eMassMatrix>::AssInertiaVec(*this, M, uP, R, oDofMap);

     ASSERT(R.iGetMaxSize() == oDofMap.iGetLocalSize());

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstPositionIndex);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
template <typename T>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssInertiaVecConsistent(const sp_grad::SpMatrix<doublereal, iNumNodes, iNumNodes>& Mcon,
                                                                                                   const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                                                                                   sp_grad::SpColVector<T, iNumDof>& R,
                                                                                                   const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumDof; ++i) {
          SpGradientTraits<T>::ResizeReset(R(i), 0., oDofMap.iGetLocalSize());

          oDofMap.InitDofMap(R(i));
     }

     for (index_type l = 1; l <= iNumNodes; ++l) {
          for (index_type k = 1; k <= iNumNodes; ++k) {
               for (index_type i = 1; i <= 3; ++i) {
                    oDofMap.Sub(R((k - 1) * 3 + i), Mcon(k, l) * uP(i, l));
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
template <typename T>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssInertiaVecLumped(const sp_grad::SpColVector<doublereal, iNumNodes>& Mlumped,
                                                                                               const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                                                                               sp_grad::SpColVector<T, iNumDof>& R,
                                                                                               const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumNodes; ++i) {
          for (index_type j = 1; j <= 3; ++j) {
               oDofMap.MapAssign(R((i - 1) * 3 + j), -Mlumped(i) * uP(j, i));
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
template <typename T>
inline void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                  doublereal dCoef,
                                                                                  const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                                  const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                                                                  enum sp_grad::SpFunctionCall func)
{
     return SolidElemCSLHelper<eConstLawType>::AssRes(this, WorkVec, dCoef, func);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
SubVectorHandler&
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssRes(SubVectorHandler& WorkVec,
                                                                                  doublereal dCoef,
                                                                                  const VectorHandler& XCurr,
                                                                                  const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemDynamic::AssRes");

     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     MassMatrixHelper<eMassMatrix>::AddInertia(*this, M);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
VariableSubMatrixHandler&
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssJac(VariableSubMatrixHandler& WorkMat,
                                                                                  doublereal dCoef,
                                                                                  const VectorHandler& XCurr,
                                                                                  const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemDynamic::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::AssJac(VectorHandler& JacY,
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

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
Vec3
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::GetB_int() const
{
     using namespace sp_grad;

     Vec3 dBeta(::Zero3);

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     SpMatrixA<doublereal, 3, iNumNodes> v;

     this->GetNodalVelocities(v, 1., SpFunctionCall::REGULAR_RES);

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);
          const SpColVector<doublereal, 3> V = v * h;
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << this->GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, this->rhon); // interpolate from nodes to collocation points

          for (index_type j = 1; j <= 3; ++j) {
               dBeta(j) += rho * detJ * alpha * V(j);
          }
     }

     return dBeta;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
doublereal
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::dGetE() const
{
     using namespace sp_grad;

     doublereal dE = 0.;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     SpMatrixA<doublereal, 3, iNumNodes> v;

     this->GetNodalVelocities(v, 1., SpFunctionCall::REGULAR_RES);

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);
          const SpColVector<doublereal, 3> V = v * h;
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << this->GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, this->rhon); // interpolate from nodes to collocation points
          const doublereal dm = rho * detJ * alpha;

          dE += 0.5 * dm * Dot(V, V);
     }

     return dE;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, MassMatrixType eMassMatrix>
Vec3
SolidElemDynamic<ElementType, CollocationType, SolidCSLType, eMassMatrix>::GetG_int() const
{
     using namespace sp_grad;

     Vec3 dGamma(::Zero3);

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     SpMatrixA<doublereal, 3, iNumNodes> x, v;

     this->GetNodalPositions(x, 1., SpFunctionCall::REGULAR_RES);
     this->GetNodalVelocities(v, 1., SpFunctionCall::REGULAR_RES);

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);
          const SpColVector<doublereal, 3> X = x * h;
          const SpColVector<doublereal, 3> V = v * h;
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << this->GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, this->rhon); // interpolate from nodes to collocation points

          const SpColVector<doublereal, 3> XtildeV = Cross(X, V);

          for (index_type j = 1; j <= 3; ++j) {
               dGamma(j) += rho * detJ * alpha * XtildeV(j);
          }
     }

     return dGamma;
}

template <typename ElementType, typename CollocationType>
SolidElem*
ReadSolid(DataManager* const pDM, MBDynParser& HP, const unsigned int uLabel)
{
     DEBUGCOUTFNAME("ReadSolid");

     using namespace sp_grad;

     typedef SolidElemStatic<ElementType, CollocationType, SolidConstLawElastic> SolidStaticElemTypeElastic;
     typedef SolidElemStatic<ElementType, CollocationType, SolidConstLawElasticIncompr> SolidStaticElemTypeElasticIncompr;
     typedef SolidElemStatic<ElementType, CollocationType, SolidConstLawViscoelastic> SolidStaticElemTypeViscoElastic;

     typedef SolidElemDynamic<ElementType, CollocationType, SolidConstLawElastic, MassMatrixType::CONSISTENT>      SolidDynamicElemTypeElasticConMass;
     typedef SolidElemDynamic<ElementType, CollocationType, SolidConstLawViscoelastic, MassMatrixType::CONSISTENT> SolidDynamicElemTypeViscoElasticConMass;
     typedef SolidElemDynamic<ElementType, CollocationType, SolidConstLawElasticIncompr, MassMatrixType::CONSISTENT> SolidDynamicElemTypeElasticIncomprConMass;

     typedef SolidElemDynamic<ElementType, CollocationType, SolidConstLawElastic, MassMatrixType::LUMPED>      SolidDynamicElemTypeElasticLumpedMass;
     typedef SolidElemDynamic<ElementType, CollocationType, SolidConstLawViscoelastic, MassMatrixType::LUMPED> SolidDynamicElemTypeViscoElasticLumpedMass;
     typedef SolidElemDynamic<ElementType, CollocationType, SolidConstLawElasticIncompr, MassMatrixType::LUMPED> SolidDynamicElemTypeElasticIncomprLumpedMass;

     constexpr SolidElemFlags eElemFlags = ElementType::eElemFlags;
     constexpr index_type iNumNodes = ElementType::iNumNodes;

     constexpr index_type iNumNodesPressure = ElementType::iNumNodesPressure;
     constexpr index_type iNumEvalPointsStiffness = CollocationType::iNumEvalPointsStiffness;

     std::array<const StructDispNodeAd*, iNumNodes> rgNodes;
     std::array<const ScalarNodeAd*, iNumNodesPressure> rgNodesPressure;

     SpColVectorA<doublereal, iNumNodes> rhon;
     const bool bStaticModel = HP.IsKeyWord("static") ? true : pDM->bIsStaticModel();
     const MassMatrixType eMassMatrixType = HP.IsKeyWord("lumped" "mass") ? MassMatrixType::LUMPED : MassMatrixType::CONSISTENT;

     for (index_type i = 0; i < iNumNodes; ++i) {
          rgNodes[i] = bStaticModel
               ? pDM->ReadNode<const StructDispNodeAd, Node::STRUCTURAL>(HP)
               : pDM->ReadNode<const DynamicStructDispNodeAd, Node::STRUCTURAL>(HP);
     }

     if constexpr(eElemFlags == SolidElemFlags::DISPLACEMENT_PRESSURE) {
          static_assert(iNumNodesPressure > 0);

          for (index_type i = 0; i < iNumNodesPressure; ++i) {
               rgNodesPressure[i] = pDM->ReadNode<const ScalarNodeAd, Node::HYDRAULIC>(HP);
          }
     } else {
          static_assert(iNumNodesPressure == 0);
     }

     for (index_type i = 1; i <= iNumNodes; ++i) {
          rhon(i) = HP.GetReal();

          if (rhon(i) <= 0.) {
               silent_cerr("solid(" << uLabel << ") density must be greater than zero at line "
                           << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }
     }

     constexpr bool bUseDisplacementPressure = eElemFlags == SolidElemFlags::DISPLACEMENT_PRESSURE;

     typedef typename std::conditional<bUseDisplacementPressure, ConstitutiveLaw7D, ConstitutiveLaw6D>::type ConstLawTplType;

     std::array<std::unique_ptr<ConstLawTplType>, iNumEvalPointsStiffness> rgMaterialData;

     ConstLawType::Type eConstLawType = ConstLawType::UNKNOWN;

     for (index_type i = 0; i < iNumEvalPointsStiffness; ++i) {
          if (i > 0 && HP.IsKeyWord("same")) {
               rgMaterialData[i].reset(rgMaterialData[i - 1]->pCopy());
               continue;
          }

          ConstLawType::Type CLType_I = ConstLawType::UNKNOWN;

          if constexpr(eElemFlags == SolidElemFlags::DISPLACEMENT_PRESSURE) {
               rgMaterialData[i].reset(HP.GetConstLaw7D(CLType_I));
          } else {
               rgMaterialData[i].reset(HP.GetConstLaw6D(CLType_I));
          }

          if (rgMaterialData[i]->iGetNumDof() != 0) {
               silent_cerr("line " << HP.GetLineData()
                           << ": solid(" << uLabel << ") "
                           "does not support dynamic constitutive laws yet"
                           << std::endl);
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          if (eConstLawType == ConstLawType::UNKNOWN) {
               switch (CLType_I) {
               case ConstLawType::ELASTIC:
               case ConstLawType::VISCOELASTIC:
                    if constexpr(eElemFlags == SolidElemFlags::DISPLACEMENT_PRESSURE) {
                         silent_cerr("solid(" << uLabel << ") does not support non incompressible constitutive law types "
                                     "at line " << HP.GetLineData() << "\n");
                         throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
                    }
                    break;
               case ConstLawType::ELASTICINCOMPR:
                    if constexpr(eElemFlags != SolidElemFlags::DISPLACEMENT_PRESSURE) {
                         silent_cerr("solid(" << uLabel << ") does not support incompressible constitutive law types "
                                     "at line " << HP.GetLineData() << "\n");
                         throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
                    }
                    break;
               default:
                    silent_cerr("solid(" << uLabel << ") does not support this constitutive law type "
                                "at line " << HP.GetLineData() << "\n");
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }
               eConstLawType = CLType_I;
          } else if (eConstLawType != CLType_I) {
               silent_cerr("solid(" << uLabel << ") all constitutive laws must have the same type "
                           "at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }
     }

     const flag fOut = pDM->fReadOutput(HP, Elem::SOLID);

     std::ostream& out = pDM->GetLogFile();

     out << ElementType::ElementName() << ": " << uLabel;

     for (sp_grad::index_type i = 0; i < iNumNodes; ++i) {
          out << ' ' << rgNodes[i]->GetLabel();
     }

     out << '\n';

     SolidElem* pEl = nullptr;
     const RigidBodyKinematics* const pRBK = pDM->pGetRBK();

     if (bStaticModel) {
          switch (eConstLawType) {
          case ConstLawType::ELASTIC:
               if constexpr(eElemFlags != SolidElemFlags::DISPLACEMENT_PRESSURE) {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidStaticElemTypeElastic,
                                           SolidStaticElemTypeElastic(uLabel,
                                                                      rgNodes,
                                                                      rgNodesPressure,
                                                                      rhon,
                                                                      std::move(rgMaterialData),
                                                                      pRBK,
                                                                      fOut));
               } else {
                    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
               }
               break;
          case ConstLawType::VISCOELASTIC:
               if constexpr(eElemFlags != SolidElemFlags::DISPLACEMENT_PRESSURE) {
                         SAFENEWWITHCONSTRUCTOR(pEl,
                                                SolidStaticElemTypeViscoElastic,
                                                SolidStaticElemTypeViscoElastic(uLabel,
                                                                                rgNodes,
                                                                                rgNodesPressure,
                                                                                rhon,
                                                                                std::move(rgMaterialData),
                                                                                pRBK,
                                                                                fOut));
               } else {
                    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
               }
               break;
          case ConstLawType::ELASTICINCOMPR:
               if constexpr(eElemFlags == SolidElemFlags::DISPLACEMENT_PRESSURE) {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidStaticElemTypeElasticIncompr,
                                           SolidStaticElemTypeElasticIncompr(uLabel,
                                                                             rgNodes,
                                                                             rgNodesPressure,
                                                                             rhon,
                                                                             std::move(rgMaterialData),
                                                                             pRBK,
                                                                             fOut));
               } else {
                    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
               }
               break;
          default:
               ASSERT(0);
               throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
          }
     } else {
          switch (eConstLawType) {
          case ConstLawType::ELASTIC:
               if constexpr(eElemFlags != SolidElemFlags::DISPLACEMENT_PRESSURE) {
                    switch (eMassMatrixType) {
                    case MassMatrixType::CONSISTENT:
                         SAFENEWWITHCONSTRUCTOR(pEl,
                                                SolidDynamicElemTypeElasticConMass,
                                                SolidDynamicElemTypeElasticConMass(uLabel,
                                                                                   rgNodes,
                                                                                   rgNodesPressure,
                                                                                   rhon,
                                                                                   std::move(rgMaterialData),
                                                                                   pRBK,
                                                                                   fOut));
                    break;
                    case MassMatrixType::LUMPED:
                         SAFENEWWITHCONSTRUCTOR(pEl,
                                                SolidDynamicElemTypeElasticLumpedMass,
                                                SolidDynamicElemTypeElasticLumpedMass(uLabel,
                                                                                      rgNodes,
                                                                                      rgNodesPressure,
                                                                                      rhon,
                                                                                      std::move(rgMaterialData),
                                                                                      pRBK,
                                                                                      fOut));
                         break;
                    default:
                         ASSERT(0);
                    }
               } else {
                    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
               }
               break;
          case ConstLawType::VISCOELASTIC:
               if constexpr(eElemFlags != SolidElemFlags::DISPLACEMENT_PRESSURE) {
                    switch (eMassMatrixType) {
                    case MassMatrixType::CONSISTENT:
                         SAFENEWWITHCONSTRUCTOR(pEl,
                                                SolidDynamicElemTypeViscoElasticConMass,
                                                SolidDynamicElemTypeViscoElasticConMass(uLabel,
                                                                                        rgNodes,
                                                                                        rgNodesPressure,
                                                                                        rhon,
                                                                                        std::move(rgMaterialData),
                                                                                        pRBK,
                                                                                        fOut));
                         break;
                    case MassMatrixType::LUMPED:
                         SAFENEWWITHCONSTRUCTOR(pEl,
                                                SolidDynamicElemTypeViscoElasticLumpedMass,
                                                SolidDynamicElemTypeViscoElasticLumpedMass(uLabel,
                                                                                           rgNodes,
                                                                                           rgNodesPressure,
                                                                                           rhon,
                                                                                           std::move(rgMaterialData),
                                                                                           pRBK,
                                                                                           fOut));
                         break;
                    default:
                         ASSERT(0);
                    }
               } else {
                    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
               }
               break;
          case ConstLawType::ELASTICINCOMPR:
               if constexpr(eElemFlags == SolidElemFlags::DISPLACEMENT_PRESSURE) {
                    switch (eMassMatrixType) {
                    case MassMatrixType::CONSISTENT:
                         SAFENEWWITHCONSTRUCTOR(pEl,
                                                SolidDynamicElemTypeElasticIncomprConMass,
                                                SolidDynamicElemTypeElasticIncomprConMass(uLabel,
                                                                                          rgNodes,
                                                                                          rgNodesPressure,
                                                                                          rhon,
                                                                                          std::move(rgMaterialData),
                                                                                          pRBK,
                                                                                          fOut));
                         break;
                    case MassMatrixType::LUMPED:
                         SAFENEWWITHCONSTRUCTOR(pEl,
                                                SolidDynamicElemTypeElasticIncomprLumpedMass,
                                                SolidDynamicElemTypeElasticIncomprLumpedMass(uLabel,
                                                                                             rgNodes,
                                                                                             rgNodesPressure,
                                                                                             rhon,
                                                                                             std::move(rgMaterialData),
                                                                                             pRBK,
                                                                                             fOut));
                         break;
                    default:
                         ASSERT(0);
                    }
               } else {
                    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
               }
               break;
          default:
               ASSERT(0);
               throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
          }
     }

     if (HP.IsArg()) {
          silent_cerr("semicolon expected "
                      "at line " << HP.GetLineData() << std::endl);
          throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     return pEl;
}

// Displacement based elements
template SolidElem* ReadSolid<Hexahedron8, Gauss2x2x2>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Hexahedron20, Gauss3x3x3>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Hexahedron27, Gauss3x3x3>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Hexahedron20r, GaussH20r>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Pentahedron15, CollocPenta15>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Tetrahedron10h, CollocTet10h>(DataManager*, MBDynParser&, unsigned int);

// Displacement/pressure formulation
template SolidElem* ReadSolid<Hexahedron8upc, Gauss2x2x2>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Hexahedron20upc, Gauss3x3x3>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Hexahedron20upcr, GaussH20r>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Pentahedron15upc, CollocPenta15>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Tetrahedron10upc, CollocTet10h>(DataManager*, MBDynParser&, unsigned int);
