/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

/* Cerniera deformabile */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "vehj2.h"

/* DeformableDispJoint - begin */

/* Costruttore non banale */
DeformableDispJoint::DeformableDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw3DOwner(pCL),
pNode1(pN1), pNode2(pN2),
tilde_f1(tilde_f1), tilde_f2(tilde_f2),
tilde_R1h(tilde_R1h), tilde_R2h(tilde_R2h),
tilde_d(0.), tilde_dPrime(0.),
bFirstRes(true), F(0.)
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

/* Distruttore */
DeformableDispJoint::~DeformableDispJoint(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
DeformableDispJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", deformable displacement joint, "
		<< pNode1->GetLabel() << ", reference, node, ",
	tilde_f1.Write(out, ", ") << ", hinge, reference, node, 1, ",
	(tilde_R1h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (tilde_R1h.GetVec(2)).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", reference, node, ",
	tilde_f2.Write(out, ", ") << ", hinge, reference, node, 1, ",
	(tilde_R2h.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (tilde_R2h.GetVec(2)).Write(out, ", ") << ", ";

	return pGetConstLaw()->Restart(out) << ';' << std::endl;
}

void
DeformableDispJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Vec3 v(GetF());
		Joint::Output(OH.Joints(), "DeformableDispJoint", GetLabel(),
	    			v, Zero3, pNode1->GetRCurr()*v, Zero3) << std::endl;
	}
}

void
DeformableDispJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					Mat3x3 R1t(pNode1->GetRCurr().Transpose());
					Vec3 f2(pNode2->GetRCurr()*tilde_f2);
  	 
					tilde_f1 = R1t*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());
	
				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					Mat3x3 R2t(pNode2->GetRCurr().Transpose());
					Vec3 f1(pNode1->GetRCurr()*tilde_f1);
  	 
					tilde_f2 = R2t*(pNode1->GetXCurr() + f1 - pNode2->GetXCurr());
	
				} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
					tilde_R1h = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*tilde_R2h;
	
				} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
					tilde_R2h = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*tilde_R1h;
	
				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}

				continue;
			}

			ConstitutiveLaw3DOwner::SetValue(pDM, X, XP, ph);
		}
	}
}

Hint *
DeformableDispJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /* } */ , sizeof("offset{" /* } */ ) - 1) == 0)
	{
		s += sizeof("offset{" /* } */ ) - 1;

		if (strcmp(&s[1], /* { */ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '2':
			return new Joint::OffsetHint<2>;
		}

	} else if (strncasecmp(s, "hinge{" /* } */, sizeof("hinge{" /* } */) - 1) == 0) {
		s += sizeof("hinge{" /* } */) - 1;

		if (strcmp(&s[1], /* { */ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}
	}

	return ConstitutiveLaw3DOwner::ParseHint(pDM, s);
}

unsigned int
DeformableDispJoint::iGetNumPrivData(void) const
{
	return 9 + ConstitutiveLaw3DOwner::iGetNumPrivData();
}

unsigned int
DeformableDispJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned idx = 0;

	switch (s[0]) {
	case 'd':
		break;

	case 'v':
		idx += 3;
		break;

	case 'F':
		idx += 6;
		break;

	default:
	{
		size_t l = sizeof("constitutiveLaw.") - 1;
		if (strncmp(s, "constitutiveLaw.", l) == 0) {
			return 9 + ConstitutiveLaw3DOwner::iGetPrivDataIdx(s + l);
		}
		return 0;
	}
	}

	switch (s[1]) {
	case 'x':
		idx += 1;
		break;
	case 'y':
		idx += 2;
		break;
	case 'z':
		idx += 3;
		break;
	default:
		return 0;
	}

	if (s[2] != '\0') {
		return 0;
	}

	return idx;
}

doublereal
DeformableDispJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());

	switch (i) {
	case 1:
	case 2:
	case 3:
	{
		/* FIXME: allows simplifications by using only the column
		 * and the components of tilde_R1h that is actually required */
		Vec3 f1(pNode1->GetRCurr()*tilde_f1);
		Vec3 f2(pNode2->GetRCurr()*tilde_f2);
		Mat3x3 R1hT((pNode1->GetRCurr()*tilde_R1h).Transpose());
		Vec3 tilde_d(R1hT*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1));
		return tilde_d(i);
	}

	case 4:
	case 5:
	case 6:
	{
		Vec3 f2(pNode2->GetRCurr()*tilde_f2);
		Mat3x3 R1hT(pNode1->GetRCurr().Transpose());
		Vec3 tilde_dPrime(R1hT*(pNode2->GetVCurr() - pNode1->GetVCurr()
					+ (pNode2->GetXCurr() - pNode1->GetXCurr()).Cross(pNode1->GetWCurr())
					- f2.Cross(pNode2->GetWCurr() - pNode1->GetWCurr())));

		return tilde_dPrime(i - 3);
	}

	case 7:
	case 8:
	case 9:
		return GetF()(i - 6);

	default:
		return ConstitutiveLaw3DOwner::dGetPrivData(i - 9);
	}
}

/* DeformableDispJoint - end */


/* ElasticDispJoint - begin */

ElasticDispJoint::ElasticDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableDispJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	FDE = R1h*GetFDE()*R1h.Transpose();
}

ElasticDispJoint::~ElasticDispJoint(void)
{
	NO_OP;
}

void
ElasticDispJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ElasticDispJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMat(WM, dCoef);

	return WorkMat;
}

/* assemblaggio jacobiano */
void
ElasticDispJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	WorkMatB.SetNullMatrix();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMat(WMA, 1.);
}

void
ElasticDispJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Vec3 f1(pNode1->GetRRef()*tilde_f1);
	Vec3 f2(pNode2->GetRRef()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());

	/* F/d */
	Mat3x3 DTmp(FDE*dCoef);

	/* Force equations */

	/* delta x1 */
	WM.Add(1, 1, DTmp);
	WM.Sub(6 + 1, 1, DTmp);

	/* delta x2 */
	WM.Sub(1, 6 + 1, DTmp);
	WM.Add(6 + 1, 6 + 1, DTmp);

	Vec3 FTmp(F*dCoef);

	/* delta g1 */
	Mat3x3 MTmp(DTmp*Mat3x3(d1) - Mat3x3(FTmp));
	WM.Sub(1, 4, MTmp);
	WM.Add(6 + 1, 4, MTmp);

	/* delta g2 */
	MTmp = DTmp*Mat3x3(f2);
	WM.Add(1, 6 + 4, MTmp);
	WM.Sub(6 + 1, 6 + 4, MTmp);

	/* Moment equation on node 1 */
	MTmp = Mat3x3(f1)*DTmp;

	/* delta x1 */
	WM.Add(4, 1, MTmp);

	/* delta x2 */
	WM.Sub(4, 6 + 1, MTmp);

	/* delta g1 */
	WM.Sub(4, 4, MTmp*Mat3x3(d1) - Mat3x3(f1.Cross(FTmp)));

	/* delta g2 */
	WM.Add(4, 6 + 4, MTmp*Mat3x3(f2));

	/* Moment equation on node 2 */
	MTmp = Mat3x3(f2)*DTmp;

	/* delta x1 */
	WM.Sub(6 + 4, 1, MTmp);

	/* delta x2 */
	WM.Add(6 + 4, 6 + 1, MTmp);

	/* delta g1 */
	WM.Add(6 + 4, 4, MTmp*Mat3x3(d1) - Mat3x3(f2, FTmp));

	/* delta g2 */
	WM.Sub(6 + 4, 6 + 4, MTmp*Mat3x3(f2) - Mat3x3(FTmp, f2));
}

/* assemblaggio residuo */
SubVectorHandler&
ElasticDispJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

void
ElasticDispJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		tilde_d = R1h.Transpose()*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1);

		ConstitutiveLaw3DOwner::Update(tilde_d);
	}

	F = R1h*GetF();

	WorkVec.Add(1, F);
	WorkVec.Add(4, f1.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, f2.Cross(F));
}

void
ElasticDispJoint::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);

	tilde_d = R1hT*(pNode2->GetXCurr() + f2 - pNode1->GetXCurr() - f1);

	ConstitutiveLaw3DOwner::Update(tilde_d);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDE = R1h*GetFDE()*R1hT;

	bFirstRes = true;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ElasticDispJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMat(WM, 1.);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ElasticDispJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* ElasticDispJoint - end */


/* ViscousDispJoint - begin */

ViscousDispJoint::ViscousDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableDispJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut)
{
#if 0
	/* Temporary */
	silent_cerr("DeformableHingeJoint(" << GetLabel() << "): "
			"this element is not implemented yet" << std::endl);
	throw ErrNotImplementedYet();
#endif

	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	FDEPrime = R1h*GetFDEPrime()*R1h.Transpose();
}

ViscousDispJoint::~ViscousDispJoint(void)
{
	NO_OP;
}

void
ViscousDispJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(Zero3, tilde_dPrime);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ViscousDispJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, dCoef);

	return WorkMat;
}

/* assemblaggio jacobiano */
void
ViscousDispJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);
	WMB.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);

		WMB.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMB.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMB.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMB.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WMA, WMB, 1.);
}

void
ViscousDispJoint::AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef)
{
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);

	/* F/dot{d} */
	WMB.Add(1, 1, FDEPrime);
	WMB.Sub(6 + 1, 1, FDEPrime);
	WMB.Sub(1, 6 + 1, FDEPrime);
	WMB.Add(6 + 1, 6 + 1, FDEPrime);

	Mat3x3 MTmp(Mat3x3(f1)*FDEPrime);

	WMB.Sub(4, 1, MTmp);
	WMB.Add(4, 6 + 1, MTmp);

	MTmp = Mat3x3(f2)*FDEPrime;

	WMB.Add(6 + 4, 1, MTmp);
	WMB.Sub(6 + 4, 6 + 1, MTmp);

	/* F/dot{d} * [ ( w1 * dCoef ) x ] */
	Mat3x3 CTmp = FDEPrime*Mat3x3(pNode1->GetWRef()*dCoef);

	WMA.Sub(1, 1, CTmp);
	WMA.Add(6 + 1, 1, CTmp);
	WMA.Add(1, 6 + 1, CTmp);
	WMA.Sub(6 + 1, 6 + 1, CTmp);

	MTmp = Mat3x3(f1)*CTmp;

	WMA.Add(4, 1, MTmp);
	WMA.Sub(4, 6 + 1, MTmp);

	MTmp = Mat3x3(f2)*CTmp;

	WMA.Sub(6 + 4, 1, MTmp);
	WMA.Add(6 + 4, 6 + 1, MTmp);

	/* F/dot{d} * [ f2 x ] */
	MTmp = FDEPrime*Mat3x3(f2);

	WMB.Add(1, 6 + 4, MTmp);
	WMB.Sub(6 + 1, 6 + 4, MTmp);

	WMB.Add(4, 6 + 4, Mat3x3(f1)*MTmp);
	WMB.Sub(6 + 4, 6 + 4, Mat3x3(f2)*MTmp);

	/* F/dot{d} * [ (x2 + f2 - x) x ] */
	Vec3 d1(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());
	MTmp = FDEPrime*Mat3x3(d1);

	WMB.Sub(1, 4, MTmp);
	WMB.Add(6 + 1, 4, MTmp);

	WMB.Sub(4, 4, Mat3x3(f1)*MTmp);
	WMB.Add(6 + 4, 4, Mat3x3(f2)*MTmp);

	/* F/dot{d} * ( [ ( f2 x w2 ) x ] + [ w1 x ] [ f2 x ] ) * dCoef */
	MTmp = FDEPrime*(Mat3x3(f2.Cross(pNode2->GetWCurr()*dCoef))
			+ Mat3x3(pNode1->GetWCurr(), f2*dCoef));
	WMA.Sub(1, 6 + 4, MTmp);
	WMA.Add(6 + 1, 6 + 4, MTmp);

	WMA.Sub(4, 6 + 4, Mat3x3(f1)*MTmp);
	WMA.Add(6 + 4, 6 + 4, Mat3x3(f2)*MTmp);

	/* F/dot{d} * ( [ d1Prime x ] - [ w1 x ] [ d1 x ] ) * dCoef */
	Vec3 d1Prime(pNode2->GetVCurr() - f2.Cross(pNode2->GetWCurr()) - pNode1->GetVCurr());
	MTmp = FDEPrime*(Mat3x3(d1Prime) - Mat3x3(pNode1->GetWCurr(), d1));
	WMB.Sub(1, 4, MTmp);
	WMB.Add(6 + 1, 4, MTmp);

	WMB.Sub(4, 6 + 4, Mat3x3(f1)*MTmp);
	WMB.Add(6 + 4, 6 + 4, Mat3x3(f2)*MTmp);

	/* - [ F x ] * dCoef */
	MTmp = Mat3x3(F*dCoef);
	WMA.Add(1, 4, MTmp);
	WMA.Sub(6 + 1, 4, MTmp);

	WMA.Add(4, 6 + 4, Mat3x3(f1)*MTmp);
	WMA.Sub(6 + 4, 6 + 4, Mat3x3(f2)*MTmp);

	/* [ F x ] [ fi x ] * dCoef */
	WMA.Sub(4, 4, Mat3x3(F*dCoef, f1));
	WMA.Add(6 + 4, 6 + 4, Mat3x3(F*dCoef, f2));
}

/* assemblaggio residuo */
SubVectorHandler&
ViscousDispJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

void
ViscousDispJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		tilde_dPrime = R1h.Transpose()*(pNode2->GetVCurr() - pNode1->GetVCurr()
				-f2.Cross(pNode2->GetWCurr())
				- (pNode2->GetXCurr() + f2 - pNode1->GetXCurr()).Cross(pNode1->GetWCurr()));

		IncrementalUpdate(Zero3, tilde_dPrime);
	}

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Add(4, f1.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, f2.Cross(F));
}

void
ViscousDispJoint::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);

	tilde_dPrime = R1h.Transpose()*(pNode2->GetVCurr() - pNode1->GetVCurr()
			-f2.Cross(pNode2->GetWCurr())
			- (pNode2->GetXCurr() + f2 - pNode1->GetXCurr()).Cross(pNode1->GetWCurr()));

	ConstitutiveLaw3DOwner::IncrementalUpdate(Zero3, tilde_dPrime);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDEPrime = R1h*GetFDEPrime()*R1hT;

	bFirstRes = true;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscousDispJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6+iCnt, iNode1FirstVelIndex + iCnt);

		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Vec3 Omega1(pNode1->GetWRef());
	Vec3 Omega2(pNode2->GetWRef());

	Vec3 F(R1h*GetF());
	Mat3x3 FDEPrime(R1h*GetFDEPrime()*R1h.Transpose());
	Mat3x3 Tmp(Mat3x3(F) - FDEPrime*Mat3x3(Omega2 - Omega1));

	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(1, 4, FDEPrime);
	WM.Add(4, 6 + 4, FDEPrime);

	WM.Sub(1, 6 + 4, FDEPrime);
	WM.Sub(4, 4, FDEPrime);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ViscousDispJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());

	Vec3 g1(pNode1->GetgCurr());
	Vec3 g2(pNode2->GetgCurr());

	/* Aggiornamento: tilde_d += R1h^T(G(g2)*g2-G(g1)*g1) */
	tilde_d = R1h.Transpose()*(g2*(4./(4. + g2.Dot()))
			- g1*(4./(4. + g1.Dot())));
	tilde_dPrime = R1h.Transpose()*(Omega2 - Omega1);
	IncrementalUpdate(tilde_d, tilde_dPrime);

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Sub(4, F);

	return WorkVec;
}

/* ViscousDispJoint - end */


/* ViscoElasticDispJoint - begin */

ViscoElasticDispJoint::ViscoElasticDispJoint(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw3D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& tilde_f1,
		const Vec3& tilde_f2,
		const Mat3x3& tilde_R1h,
		const Mat3x3& tilde_R2h,
		flag fOut)
: Elem(uL, fOut),
DeformableDispJoint(uL, pDO, pCL, pN1, pN2, tilde_f1, tilde_f2, tilde_R1h, tilde_R2h, fOut)
{
#if 0
	/* Temporary */
	silent_cerr("DeformableHingeJoint(" << GetLabel() << "): "
			"this element is not implemented yet" << std::endl);
	throw ErrNotImplementedYet();
#endif

	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	FDE = R1h*GetFDE()*R1hT;
	FDEPrime = R1h*GetFDEPrime()*R1hT;
}

ViscoElasticDispJoint::~ViscoElasticDispJoint(void)
{
	NO_OP;
}

void
ViscoElasticDispJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw3DOwner::AfterConvergence(tilde_d, tilde_dPrime);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ViscoElasticDispJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WM, WM, dCoef);

	return WorkMat;
}

/* assemblaggio jacobiano */
void
ViscoElasticDispJoint::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);
	WMB.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMA.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMA.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);

		WMB.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WMB.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WMB.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WMB.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssMats(WMA, WMB, 1.);
}

/* assemblaggio jacobiano */
void
ViscoElasticDispJoint::AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef)
{

	Vec3 f1(pNode1->GetRRef()*tilde_f1);
	Vec3 f2(pNode2->GetRRef()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());
	Vec3 d1Prime(pNode2->GetVCurr() - f2.Cross(pNode2->GetWCurr())
			- pNode1->GetVCurr());
 
	/* F/d */
	Mat3x3 DTmp(FDE*dCoef);

	/* Force equations */

	/* delta x1 */
	WMA.Add(1, 1, DTmp);
	WMA.Sub(6 + 1, 1, DTmp);

	/* delta x2 */
	WMA.Sub(1, 6 + 1, DTmp);
	WMA.Add(6 + 1, 6 + 1, DTmp);

	/* delta g1 */
	Mat3x3 MTmp(DTmp*Mat3x3(d1));
	WMA.Sub(1, 4, MTmp);
	WMA.Add(6 + 1, 4, MTmp);

	/* delta g2 */
	MTmp = DTmp*Mat3x3(f2);
	WMA.Add(1, 6 + 4, MTmp);
	WMA.Sub(6 + 1, 6 + 4, MTmp);

	/* Moment equation on node 1 */
	MTmp = Mat3x3(f1)*DTmp;

	/* delta x1 */
	WMA.Add(4, 1, MTmp);

	/* delta x2 */
	WMA.Sub(4, 6 + 1, MTmp);

	/* delta g1 */
	WMA.Sub(4, 4, MTmp*Mat3x3(d1));

	/* delta g2 */
	WMA.Add(4, 6 + 4, MTmp*Mat3x3(f2));

	/* Moment equation on node 2 */
	MTmp = Mat3x3(f2)*DTmp;

	/* delta x1 */
	WMA.Sub(6 + 4, 1, MTmp);

	/* delta x2 */
	WMA.Add(6 + 4, 6 + 1, MTmp);

	/* delta g1 */
	WMA.Add(6 + 4, 4, MTmp*Mat3x3(d1));

	/* delta g2 */
	WMA.Sub(6 + 4, 6 + 4, MTmp*Mat3x3(f2));

	/* F/dot{d} */
	WMB.Add(1, 1, FDEPrime);
	WMB.Sub(6 + 1, 1, FDEPrime);
	WMB.Sub(1, 6 + 1, FDEPrime);
	WMB.Add(6 + 1, 6 + 1, FDEPrime);

	MTmp = Mat3x3(f1)*FDEPrime;

	WMB.Add(4, 1, MTmp);
	WMB.Sub(4, 6 + 1, MTmp);

	MTmp = Mat3x3(f2)*FDEPrime;

	WMB.Sub(6 + 4, 1, MTmp);
	WMB.Add(6 + 4, 6 + 1, MTmp);

	/* F/dot{d} * [ ( w1 * dCoef ) x ] */
	Mat3x3 CTmp = FDEPrime*Mat3x3(pNode1->GetWRef()*dCoef);

	WMA.Sub(1, 1, CTmp);
	WMA.Add(6 + 1, 1, CTmp);
	WMA.Add(1, 6 + 1, CTmp);
	WMA.Sub(6 + 1, 6 + 1, CTmp);

	MTmp = Mat3x3(f1)*CTmp;

	WMA.Sub(4, 1, MTmp);
	WMA.Add(4, 6 + 1, MTmp);

	MTmp = Mat3x3(f2)*CTmp;

	WMA.Add(6 + 4, 1, MTmp);
	WMA.Sub(6 + 4, 6 + 1, MTmp);

	/* F/dot{d} * [ f2 x ] */
	MTmp = FDEPrime*Mat3x3(f2);

	WMB.Add(1, 6 + 4, MTmp);
	WMB.Sub(6 + 1, 6 + 4, MTmp);

	WMB.Add(4, 6 + 4, Mat3x3(f1)*MTmp);
	WMB.Sub(6 + 4, 6 + 4, Mat3x3(f2)*MTmp);

	/* F/dot{d} * [ (x2 + f2 - x) x ] */
	MTmp = FDEPrime*Mat3x3(d1);

	WMB.Sub(1, 4, MTmp);
	WMB.Add(6 + 1, 4, MTmp);

	WMB.Sub(4, 4, Mat3x3(f1)*MTmp);
	WMB.Add(6 + 4, 4, Mat3x3(f2)*MTmp);

	/* F/dot{d} * ( [ ( f2 x w2 ) x ] + [ w1 x ] [ f2 x ] ) * dCoef */
	MTmp = FDEPrime*(Mat3x3(f2.Cross(pNode2->GetWCurr()*dCoef))
			+ Mat3x3(pNode1->GetWCurr(), f2*dCoef));
	WMA.Sub(1, 6 + 4, MTmp);
	WMA.Add(6 + 1, 6 + 4, MTmp);

	WMA.Sub(4, 6 + 4, Mat3x3(f1)*MTmp);
	WMA.Add(6 + 4, 6 + 4, Mat3x3(f2)*MTmp);

	/* F/dot{d} * ( [ d1Prime x ] - [ w1 x ] [ d1 x ] ) * dCoef */
	MTmp = FDEPrime*(Mat3x3(d1Prime*dCoef) - Mat3x3(pNode1->GetWCurr(), d1*dCoef));
	WMA.Sub(1, 4, MTmp);
	WMA.Add(6 + 1, 4, MTmp);

	WMA.Sub(4, 6 + 4, Mat3x3(f1)*MTmp);
	WMA.Add(6 + 4, 6 + 4, Mat3x3(f2)*MTmp);

	Vec3 FTmp(F*dCoef);

	/* - [ F x ] * dCoef */
	MTmp = Mat3x3(FTmp);
	WMA.Add(1, 4, MTmp);
	WMA.Sub(6 + 1, 4, MTmp);

	WMA.Add(4, 4, Mat3x3(f1.Cross(FTmp)));
	WMA.Sub(6 + 4, 4, Mat3x3(f2)*MTmp);

	/* [ F x ] [ fi x ] * dCoef */
	WMA.Add(6 + 4, 6 + 4, Mat3x3(FTmp, f2));
}

/* assemblaggio residuo */
SubVectorHandler&
ViscoElasticDispJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);
	
	return WorkVec;
}

void
ViscoElasticDispJoint::AssVec(SubVectorHandler& WorkVec)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R1hT(R1h.Transpose());
		Vec3 d1(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());
		Vec3 d1Prime(pNode2->GetVCurr() - f2.Cross(pNode2->GetWCurr())
				- pNode1->GetVCurr());

		tilde_d = R1hT*(d1 - f1);
		tilde_dPrime = R1hT*(d1Prime - d1.Cross(pNode1->GetWCurr()));

		ConstitutiveLaw3DOwner::Update(tilde_d, tilde_dPrime);
	}

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Add(4, f1.Cross(F));
	WorkVec.Sub(6 + 1, F);
	WorkVec.Sub(6 + 4, f2.Cross(F));
}

void
ViscoElasticDispJoint::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 f1(pNode1->GetRCurr()*tilde_f1);
	Vec3 f2(pNode2->GetRCurr()*tilde_f2);
	Vec3 d1(pNode2->GetXCurr() + f2 - pNode1->GetXCurr());
	Vec3 d1Prime(pNode2->GetVCurr() - f2.Cross(pNode2->GetWCurr())
			- pNode1->GetVCurr());
	
	tilde_d = R1hT*(d1 - f1);
	tilde_dPrime = R1hT*(d1Prime - d1.Cross(pNode1->GetWCurr()));

	ConstitutiveLaw3DOwner::Update(tilde_d, tilde_dPrime);

	/* FIXME: we need to be able to regenerate FDE and FDEPrime
	 * if the constitutive law throws ChangedEquationStructure */
	FDE = R1h*GetFDE()*R1hT;
	FDEPrime = R1h*GetFDEPrime()*R1hT;

	bFirstRes = true;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscoElasticDispJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);

		WM.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRRef()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 Omega1(pNode1->GetWRef());
	Vec3 Omega2(pNode2->GetWRef());

	Vec3 F(R1h*GetF());
	Mat3x3 FDE(R1h*GetFDE()*R1hT);
	Mat3x3 FDEPrime(R1h*GetFDEPrime()*R1hT);

	Mat3x3 Tmp(Mat3x3(F) - FDEPrime*Mat3x3(Omega2-Omega1)+FDE);

	WM.Add(1, 1, Tmp);
	WM.Sub(4, 1, Tmp);

	WM.Add(4, 6 + 1, FDE);
	WM.Sub(1, 6 + 1, FDE);

	WM.Add(1, 4, FDEPrime);
	WM.Add(4, 6 + 4, FDEPrime);

	WM.Sub(1, 6 + 4, FDEPrime);
	WM.Sub(4, 4, FDEPrime);

	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ViscoElasticDispJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(3 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	Mat3x3 R1h(pNode1->GetRCurr()*tilde_R1h);
	Mat3x3 R1hT(R1h.Transpose());
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());

	Vec3 g1(pNode1->GetgCurr());
	Vec3 g2(pNode2->GetgCurr());

	/* Aggiornamento: tilde_d += R1h^T(G(g2)*g2-G(g1)*g1) */
	tilde_d = R1hT*(g2*(4./(4.+g2.Dot())) - g1*(4./(4.+g1.Dot())));
	tilde_dPrime = R1hT*(Omega2-Omega1);
	IncrementalUpdate(tilde_d, tilde_dPrime);

	Vec3 F(R1h*GetF());

	WorkVec.Add(1, F);
	WorkVec.Sub(4, F);

	return WorkVec;
}

/* ViscoElasticDispJoint - end */
