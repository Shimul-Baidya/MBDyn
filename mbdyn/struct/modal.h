/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
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

/* Elemento modale */

/* 
 * Copyright 1999-2023 Felice Felippone <ffelipp@tin.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

/* 
 * Copyright 1999-2023 Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 * Modified by Pierangelo Masarati
 */

#ifndef MODAL_H
#define MODAL_H

#include <array>
#include <fstream>
#include <joint.h>

#if 0
#define MODAL_USE_INV9
#endif

/* Modal - begin */

/* 
 * ATTENZIONE! 
 * per ora e' derivato da Joint; 
 * puo' darsi che venga creata una classe apposta
 */

class Modal : virtual public Elem, public Joint {
public:
	struct StrNodeData;
protected:
	const ModalNode* const pModalNode;
	const unsigned iRigidOffset;		/* 0 iff pModalNode == 0; else 12 */

	/* configuration of reference point;
	 * from ModalNode iff pModalNode == 0 */
	mutable Vec3	x;
	mutable Mat3x3	R;
	mutable Mat3x3	RT;

	const unsigned int NModes;
	const unsigned int NStrNodes;

	const unsigned int NFEMNodes; // number of FEM nodes, common
	const std::vector<std::string> IdFEMNodes; // ID of FEM nodes, common
	const Mat3xN oXYZFEMNodes; // local position of FEM nodes, common
	const doublereal dMass; // mass, common
	const Vec3 Inv2; // undeformed static moment, common
	const Mat3x3 Inv7; // undeformed inertia moment, common

	const std::vector<unsigned int> uModeNumber;
	const MatNxN oModalMass;
	const MatNxN oModalStiff;
	const MatNxN oModalDamp;
	const Mat3xN oPHIt;
	const Mat3xN oPHIr;
   
	const Mat3xN oModeShapest;
	const Mat3xN oModeShapesr;

	Mat3xN oCurrXYZ;
	Mat3xN oCurrXYZVel;

	const Mat3xN oInv3;
	const Mat3xN oInv4;
	const Mat3xN oInv5;
	const Mat3xN oInv8;
	const Mat3xN oInv9;

	const Mat3xN oInv10;
	const Mat3xN oInv11;

	Vec3   Inv3jaj;
	Vec3   Inv3jaPj;
	Mat3x3 Inv8jaj;

	Mat3x3 Inv8jaPj;
	Mat3xN Inv5jaj;
	Mat3xN Inv5jaPj;
     
	Mat3x3 Inv9jkajak;
	Mat3x3 Inv9jkajaPk;
     
	VecN a, a0;
	VecN aPrime, aPrime0;
	VecN b;
	VecN bPrime;

public:
        template <unsigned N>
        class StressStiffIndex {
        public:
             static_assert(N > 0);

             StressStiffIndex()
                  :uSize(0u) {
#ifdef DEBUG
                  std::fill(std::begin(rgIndexMat), std::end(rgIndexMat), SS_UNUSED);
                  std::fill(std::begin(rgIndexVec), std::end(rgIndexVec), SS_UNUSED);
#endif
             }

             void Insert(unsigned uIndexMat, unsigned uIndexVec) {
                  ASSERT(uSize >= 0);
                  ASSERT(uSize < N);

                  ASSERT(rgIndexMat[uSize] == SS_UNUSED);
                  ASSERT(rgIndexVec[uSize] == SS_UNUSED);

                  rgIndexMat[uSize] = uIndexMat;
                  rgIndexVec[uSize] = uIndexVec;

                  ++uSize;

                  ASSERT(uSize <= N);
             }

             unsigned uGetSize() const {
                  ASSERT(uSize >= 0);
                  ASSERT(uSize <= N);

                  return uSize;
             }

             unsigned uGetIndexMat(unsigned iIndex) const {
                  ASSERT(iIndex >= 0);
                  ASSERT(iIndex < uSize);
                  ASSERT(uSize <= N);
                  ASSERT(rgIndexMat[iIndex] != SS_UNUSED);
                  return rgIndexMat[iIndex];
             }

             unsigned uGetIndexVec(unsigned iIndex) const {
                  ASSERT(iIndex >= 0);
                  ASSERT(iIndex < uSize);
                  ASSERT(uSize <= N);
                  ASSERT(rgIndexVec[iIndex] != SS_UNUSED);

                  return rgIndexVec[iIndex];
             }

        private:
             enum: unsigned { SS_UNUSED = ~0u };
             unsigned uSize;
             std::array<unsigned, N> rgIndexMat;
             std::array<unsigned, N> rgIndexVec;
        };
     
	struct StrNodeData {
                StrNodeData();
		// constant, defined once for all at input
		const StructNode *pNode;
                const class StructNodeAd *pNodeAd;
		std::string FEMNode;
		Vec3 OffsetFEM;
		Vec3 OffsetMB;
		Mat3x3 RotMB;

		// variable, constructed during analysis
		Vec3 d1tot;
		Mat3x3 R1tot;
		Mat3x3 R2;

		// variable, constructed during analysis
		Vec3 F;
		Vec3 M;
                StressStiffIndex<3> oStressStiffIndexF;
                StressStiffIndex<3> oStressStiffIndexM;
                bool bOut;
        };

protected:
	std::vector<StrNodeData> SND;

	/* from gravity.h */
	/* momento statico */
	Vec3 GetS_int(void) const;

	/* momento d'inerzia */
	Mat3x3 GetJ_int(void) const;
 
	Vec3 GetB_int(void) const;
	Vec3 GetG_int(void) const;

public:
	/* Costruttore non banale */
	Modal(unsigned int uL,
			const ModalNode* pModalNodeTmp, 
			const Vec3& x0,
			const Mat3x3& R0,
			const DofOwner* pDO,
			unsigned int N,
			unsigned int NS,
			unsigned int NFN,
			doublereal dMass,
			const Vec3& STmp,
			const Mat3x3& JTmp,
			std::vector<unsigned int>&& uModeNumber,
			MatNxN&& oGenMass,
			MatNxN&& oGenStiff,
			MatNxN&& oGenDamp,
			std::vector<std::string>&& IdFEMNodes,
			Mat3xN&& oN,
			std::vector<Modal::StrNodeData>&& snd,
			Mat3xN&& oPHIt,
			Mat3xN&& oPHIr,
			Mat3xN&& oModeShapest,
			Mat3xN&& oModeShapesr,
			Mat3xN&& oInv3,
			Mat3xN&& oInv4,
			Mat3xN&& oInv5,
			Mat3xN&& oInv8,
			Mat3xN&& oInv9,
			Mat3xN&& oInv10,
			Mat3xN&& oInv11,
			VecN&& a,
			VecN&& aP,
			flag fOut);

	/* Distruttore */
	~Modal(void);
   
	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const;

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual unsigned int iGetNumDof(void) const;
	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType(unsigned int i) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	void Output(OutputHandler& OH) const;

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);   
	/* Contributo al residuo durante l'assemblaggio iniziale */   
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
 
	/* Setta il valore iniziale delle proprie variabili */
	void SetInitialValue(VectorHandler& /* X */ );

	void SetValue(DataManager *pDM,
			VectorHandler& /* X */ , VectorHandler& /* XP */ ,
			SimulationEntity::Hints *ph = 0);

#if 0
	/* Aggiorna dati durante l'iterazione fittizia iniziale */
	virtual void DerivativesUpdate(const VectorHandler& X,
		const VectorHandler& XP);
#endif

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;   

	/* Funzioni che restituiscono dati che possono servire ad
	 * altri elementi (ad es. agli elementi aerodinamici modali)
	 */

	const Mat3xN& pGetPHIt(void) const {
		return oModeShapest;
	};

	const Mat3xN& pGetPHIr(void) const {
		return oModeShapesr;
	};

	// NOTE: not 'const' because modify internal storage
	const Mat3xN& GetCurrFEMNodesPosition(void);
	const Mat3xN& GetCurrFEMNodesVelocity(void);

	integer uGetNModes(void) const {
		return NModes;
	};

	const std::vector<unsigned int>& GetModeList(void) const {
		return uModeNumber;
	};

	const VecN& GetA(void) const {
		return a;
	};

	const VecN& GetAP(void) const {
		return aPrime;
	};

	const VecN& GetB(void) const {
		return b;
	};

	const VecN& GetBP(void) const {
		return bPrime;
	};

	integer uGetNFEMNodes(void) {
		return NFEMNodes;
	};

	integer iGetModalIndex(void) const {
		return iGetFirstIndex();
	};

	const ModalNode* pGetModalNode(void) const {
		return pModalNode;
	};

	/* from gravity.h */
	/* massa totale */
	doublereal dGetM(void) const;

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(NStrNodes + (pModalNode ? 1 : 0));
		for (unsigned int j = 0; j < NStrNodes; j++) {
			connectedNodes[j] = SND[j].pNode;
		}
		if (pModalNode) {
			connectedNodes[NStrNodes] = pModalNode;
		}
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* Modal - end */

class DataManager;
class MBDynParser;

extern Joint *
ReadModal(DataManager* pDM, MBDynParser& HP, const DofOwner* pD0,
		unsigned int uLabel);

#endif /* MODAL_H */

