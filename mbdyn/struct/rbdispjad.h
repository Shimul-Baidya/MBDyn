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

#ifndef RIGIDBODYDISPJOINTAD_H
#define RIGIDBODYDISPJOINTAD_H

#include "joint.h"
#include <sp_gradient.h>
#include <sp_matrix_base.h>
#include <sp_matvecass.h>
#include "strnodead.h"

class RigidBodyDispJointAd: public Joint {
public:
     struct SlaveNodeData {
          SlaveNodeData(const StructDispNodeAd* pNode,
                        const Vec3& offset,
                        doublereal weight)
               :pNode{pNode},
                offset{offset},
                weight{weight} {
                }
          const StructDispNodeAd* pNode;
          Vec3 offset;
          doublereal weight;
     };

     RigidBodyDispJointAd(unsigned int uL,
                          const DofOwner* pD,
                          const StructNodeAd* pNodeMaster,
                          std::vector<SlaveNodeData>&& rgNodesSlave,
                          flag fOut);
     virtual ~RigidBodyDispJointAd();

     virtual void Output(OutputHandler& OH) const override;
     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
     virtual unsigned int iGetNumDof(void) const override;
     virtual DofOrder::Order GetDofType(unsigned int i) const override;
     virtual DofOrder::Order GetEqType(unsigned int i) const override;
     virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const override;
     virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const override;
     VariableSubMatrixHandler&
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
     SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;
     unsigned int iGetNumPrivData(void) const override;
     virtual unsigned int iGetPrivDataIdx(const char *s) const override;
     virtual doublereal dGetPrivData(unsigned int i) const override;
     int GetNumConnectedNodes(void) const override;
     void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override;
     void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                   SimulationEntity::Hints *ph) override;
     std::ostream& Restart(std::ostream& out) const override;
     virtual unsigned int iGetInitialNumDof(void) const override;
     virtual void
     InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
     VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
                   const VectorHandler& XCurr) override;
     SubVectorHandler&
     InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr) override;

     virtual const OutputHandler::Dimensions
     GetEquationDimension(integer index) const override;

     virtual Type GetJointType() const override;

     using Elem::AssRes;

     template <typename T>
     void AssRes(sp_grad::SpGradientAssVec<T>& WorkMat,
                 doublereal dCoef,
                 const sp_grad::SpGradientVectorHandler<T>& XCurr,
                 const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                 sp_grad::SpFunctionCall func);

     template <typename T>
     void InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkMat,
                        const sp_grad::SpGradientVectorHandler<T>& XCurr,
                        sp_grad::SpFunctionCall func);
private:
     void SaveReactionForce(const sp_grad::SpColVector<doublereal, 3>& Fm,
                            const sp_grad::SpColVector<doublereal, 3>& Mm);
     void SaveReactionForce(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&,
                            const sp_grad::SpColVector<sp_grad::SpGradient, 3>&) {}
     void SaveReactionForce(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&,
                            const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&) {}
     const StructNodeAd* const pNodeMaster;
     const std::vector<SlaveNodeData> rgNodesSlave;
     Vec3 FmTmp;
     Vec3 MmTmp;
};

#endif
