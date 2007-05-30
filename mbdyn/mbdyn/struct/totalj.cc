/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

/* TotalJoint
 * Authors: Alessandro Fumagalli, Pierangelo Masarati
 *
 * */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/iostream>
#include <fstream>

#include "totalj.h"
#include "Rot.hh"
#include "hint_impl.h"

static const char idx2xyz[] = { 'x', 'y', 'z' };

/* TotalJoint - begin */

TotalJoint::TotalJoint(unsigned int uL, const DofOwner *pDO,
	bool bPos[3],
	TplDriveCaller<Vec3> *const pDCPos[3],
	bool bRot[3],
	TplDriveCaller<Vec3> *const pDCRot[3],
	const StructNode *pN1,
	const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp,
	const StructNode *pN2,
	const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2),
f1(f1Tmp), R1h(R1hTmp), R1hr(R1hrTmp),
f2(f2Tmp), R2h(R2hTmp), R2hr(R2hrTmp),
XDrv(pDCPos[0]), XPDrv(pDCPos[1]), XPPDrv(pDCPos[2]),
ThetaDrv(pDCRot[0]), OmegaDrv(pDCRot[1]), OmegaPDrv(pDCRot[2]),
nConstraints(0), nPosConstraints(0), nRotConstraints(0),
tilde_f1(R1h.Transpose()*f1),
M(0.), F(0.), ThetaDelta(0.), ThetaDeltaPrev(0.)
{
	/* Equations 1->3: Positions
	 * Equations 4->6: Rotations */

	for (unsigned int i = 0; i < 3; i++) {
		bPosActive[i] = bPos[i];
		bRotActive[i] = bRot[i];
		if (bPosActive[i]) {
			iPosIncid[nPosConstraints] = i + 1;
			nPosConstraints++;
		}
		if (bRotActive[i]) {
			iRotIncid[nRotConstraints] = i + 1;
			nRotConstraints++;
		}
	}
	nConstraints = nPosConstraints + nRotConstraints;
}

TotalJoint::~TotalJoint(void)
{
	NO_OP;
};

std::ostream&
TotalJoint::DescribeDof(std::ostream& out,
	char *prefix, bool bInitial, int i) const
{
	integer iIndex = iGetFirstIndex();

	if (i >= 0) {
		silent_cerr("TotalJoint(" << GetLabel() << "): "
			"DescribeDof(" << i << ") "
			"not implemented yet" << std::endl);
		throw ErrGeneric();
	}

	if (nPosConstraints > 1) {
		out << prefix << iIndex + 1;
		out << "->" << iIndex + nPosConstraints;
		out << ": ";
	}
	out << "reaction force(s) [";

	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bPosActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "F" << idx2xyz[i];
		}
	}
	out << "]" << std::endl;
	

	if (nRotConstraints > 1) {
		out << prefix << iIndex + nPosConstraints + 1;
		out << "->" << iIndex + nConstraints ;
		out << ": ";
	}
	out << "reaction couple(s) [";

	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bRotActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "m" << idx2xyz[i];
		}
	}
	out << "]" << std::endl;

	if (bInitial) {
		iIndex += nConstraints;

		if (nPosConstraints > 1) {
			out << prefix << iIndex + 1;
			out << "->" << iIndex + nPosConstraints;
			out << ": ";
		}
		out << "reaction force(s) derivative(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bPosActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "FP" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
		

		if (nRotConstraints > 1) {
			out << prefix << iIndex + nPosConstraints + 1;
			out << "->" << iIndex + nConstraints;
			out << ": ";
		}
		out << "reaction couple(s) derivative(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bRotActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "mP" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
	}
	return out;
}

std::ostream&
TotalJoint::DescribeEq(std::ostream& out,
	char *prefix, bool bInitial, int i) const
{
	integer iIndex = iGetFirstIndex();

	if (i >= 0) {
		silent_cerr("TotalJoint(" << GetLabel() << "): "
			"DescribeEq(" << i << ") "
			"not implemented yet" << std::endl);
		throw ErrGeneric();
	}

	if (nPosConstraints > 1) {
		out << prefix << iIndex + 1;
		out << "->" << iIndex + nPosConstraints;
		out << ": ";
	}
	
	if (nPosConstraints > 0) {
		out << "position constraint(s) [";
	}
	
	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bPosActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "P" << idx2xyz[i] << "1=P" << idx2xyz[i] << "2";
		}
	}
	
	if (nPosConstraints > 0) {
		out << "]" << std::endl;
	}

	if (nRotConstraints > 1) {
		out << prefix << iIndex + nPosConstraints + 1;
		out << "->" << iIndex + nConstraints ;
		out << ": ";
	}
	if (nRotConstraints > 0) {
		out << "orientation constraint(s) [";
	}
	
	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bRotActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "g" << idx2xyz[i] << "1=g" << idx2xyz[i] << "2";
		}
	}
	if (nRotConstraints > 0) {
		out << "]" << std::endl;
	}

	if (bInitial) {
		iIndex += nConstraints;

		if (nPosConstraints > 1) {
			out << prefix << iIndex + 1;
			out << "->" << iIndex + nPosConstraints;
			out << ": ";
		}
		
		if (nPosConstraints > 0) {
			out << "velocity constraint(s) [";
		}
		
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bPosActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "v" << idx2xyz[i] << "1=v" << idx2xyz[i] << "2";
			}
		}
		if (nPosConstraints > 0) {
			out << "]" << std::endl;
		}
		
		if (nRotConstraints > 1) {
			out << prefix << iIndex + nPosConstraints + 1;
			out << "->" << iIndex + nConstraints ;
			out << ": ";
		}
		if (nRotConstraints > 0) {
			out << "angular velocity constraint(s) [";
		}
		
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bRotActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "w" << idx2xyz[i] << "1=w" << idx2xyz[i] << "2";
			}
		}
		if (nRotConstraints > 0) {
			out << "]" << std::endl;
		}	
	}

	return out;
}

void
TotalJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned int i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {

				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					Mat3x3 R1T(pNode1->GetRCurr().Transpose());
					Vec3 fTmp2(pNode2->GetRCurr()*f2);

					f1 = R1T*(pNode2->GetXCurr() + fTmp2 - pNode1->GetXCurr());
					tilde_f1 = R1h.Transpose()*f1;
					

				} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
					Mat3x3 R2T(pNode2->GetRCurr().Transpose());
					Vec3 fTmp1(pNode1->GetRCurr()*f1);

					f2 = R2T*(pNode1->GetXCurr() + fTmp1 - pNode2->GetXCurr());

				} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<1> *>(pjh)) {
						R1h = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*R2h;

					} else if (dynamic_cast<Joint::OrientationHingeHint<1> *>(pjh)) {
						R1hr = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*R2hr;
					}

				} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<2> *>(pjh)) {
						R2h = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*R1h;

					} else if (dynamic_cast<Joint::OrientationHingeHint<2> *>(pjh)) {
						R2hr = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*R1hr;
					}

				} else if (dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh)) {
					Joint::JointDriveHint<Vec3> *pjdh
						= dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh);
					pedantic_cout("TotalJoint(" << uLabel << "): "
						"creating drive from hint[" << i << "]..." << std::endl);

					TplDriveCaller<Vec3> *pDC = pjdh->pTDH->pCreateDrive(pDM);
					if (pDC == 0) {
						silent_cerr("TotalJoint(" << uLabel << "): "
							"unable to create drive "
							"after hint #" << i << std::endl);
						throw ErrGeneric();
					}

					if (dynamic_cast<Joint::PositionDriveHint<Vec3> *>(pjdh)) {
						XDrv.Set(pDC);

					} else if (dynamic_cast<Joint::VelocityDriveHint<Vec3> *>(pjdh)) {
						XPDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AccelerationDriveHint<Vec3> *>(pjdh)) {
						XPPDrv.Set(pDC);

					} else if (dynamic_cast<Joint::OrientationDriveHint<Vec3> *>(pjdh)) {
						ThetaDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AngularVelocityDriveHint<Vec3> *>(pjdh)) {
						OmegaDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AngularAccelerationDriveHint<Vec3> *>(pjdh)) {
						OmegaPDrv.Set(pDC);

					} else {
						delete pDC;
					}

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}
		}
	}
}

Hint *
TotalJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/ , STRLENOF("offset{" /*}*/ )) == 0)
	{
		s += STRLENOF("offset{" /*}*/ );

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '2':
			return new Joint::OffsetHint<2>;
		}

	} else if (strncasecmp(s, "position-hinge{" /*}*/, STRLENOF("position-hinge{" /*}*/)) == 0) {
		s += STRLENOF("position-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}

	} else if (strncasecmp(s, "position-drive3{" /*}*/, STRLENOF("position-drive3{" /*}*/)) == 0) {
		s += STRLENOF("position-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new PositionDriveHint<Vec3>(pTDH);
			}
		}
		return 0;

	} else if (strncasecmp(s, "orientation-hinge{" /*}*/, STRLENOF("orientation-hinge{" /*}*/)) == 0) {
		s += STRLENOF("orientation-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}

	} else if (strncasecmp(s, "orientation-drive3{" /*}*/, STRLENOF("orientation-drive3{" /*}*/)) == 0) {
		s += STRLENOF("orientation-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new OrientationDriveHint<Vec3>(pTDH);
			}
		}
		return 0;
	}

	return 0;
}

void
TotalJoint::AfterConvergence(const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */ )
{
	ThetaDeltaPrev = Unwrap(ThetaDeltaPrev, ThetaDelta);
}

void
TotalJoint::AfterConvergence(const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */, const VectorHandler& /* XP */)
{
	ThetaDeltaPrev = Unwrap(ThetaDeltaPrev, ThetaDelta);
}

/* Contributo al file di restart */
std::ostream&
TotalJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", total joint, "
		<< pNode1->GetLabel() << ", "
			<< "position, " << f1.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", R1h.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R1h.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", R1hr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R1hr.GetVec(2).Write(out, ", ") << ", "
		<< pNode2->GetLabel() << ", "
			<< "position, " << f2.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", R2h.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R2h.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", R2hr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", R2hr.GetVec(2).Write(out, ", ");

	if (bPosActive[0] || bPosActive[1] || bPosActive[2]) {
		out << ", position constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bPosActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}
	}

	if (bRotActive[0] || bRotActive[1] || bRotActive[2]) {
		out << ", orientation constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bRotActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}
	}

	return out << std::endl;
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
TotalJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	/*
	 * Constraint Equations:
	 * Position: 	R1^T(x2 + R2*f2 -x1 - R1*f1) - d = x^delta
	 * 		==> d(Vec3) = imposed displacement in node 1 local R.F.
	 * 		==> x^delta is used to activate/deactivate the constraint
	 * 			equation along the corresponding direction. If each
	 * 			component is set to 0, all relative displacement
	 * 			are forbidden (or imposed by the drive). If a component
	 * 			of x^delta is left free, the corrsponding equation is
	 * 			dropped
	 *
	 * Orientation:	Theta - Theta0 = ax(exp^-1(R1^T * R2 * R0^T)) = Theta^delta
	 * 		==> Theta = ax(exp^-1(R1^T * R2)) = Relative orientation in node1 R.F.
	 * 		==> Theta0 = Imposed relative orientation = ax(exp^-1(R0))
	 * 		==> Theta^delta is used to activate/deactivate the constraint
	 * 			equation along the corresponding direction. If each
	 * 			component is set to 0, all relative rotation
	 * 			are forbidden (or imposed by the drive). If a component
	 * 			of Theta^delta is left free, the corrsponding equation is
	 * 			dropped
	 *Jacobian Matrix:
	 *       x1  	     g1       	   x2    	g2       	 F	      M
	 * Q1 |  0   	     F1X           0            0              -R1            0	 | | x1 |
	 * G1 |-(F1)X  (b1)X(F1)X+(M1)X  (F1)X     -(F1)X(b2)X       (b1)X(R1)      -R1r | | g1 |
	 * Q2 |  0          -F1X           0    	0   	         R1	      0  | | x2 |
	 * G2 |  0    -(b2)X(F1)X-(M1)X    0        (F1)X(b2)X       (b2)X(R1)       R1r | | g2 |
	 * F  |-c*R1^T  c*R1^T*[(b1)X]   c*R1^T   -c*R1^T*[(b2)X]        0            0	 | | F  |if(bPos)
	 * M  |  0        -c*R1r^T         0         c* R1r^T            0   	      0	 | | M  |if(bRot)
	 *           	                                               if(bPos)    if(bRot)
	 *with: _ b1 = (x2 + R2*f2 - x1)
	 *      _ b2 = (R2*f2)
	 *      _ R1 = R1*R1h
	 *      _ R2 = R2*R2h
	 *      _ R1r = R1*R1hr
	 *      _ F1 = R1*F
	 *      _ M1 = R1*M
	 *      _ X = "Cross" operator
	 *
	 *     */

	DEBUGCOUT("Entering TotalJoint::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	/* Moltiplica il momento e la forza per il coefficiente del metodo */
	Vec3 FTmp(R1*(F*dCoef));
	Vec3 MTmp(R1r*(M*dCoef));

	/* Equilibrium: ((Phi/q)^T*Lambda)/q */

	Mat3x3 Tmp;

	/* [ F x ] */
	Tmp = Mat3x3(FTmp);

	/* Lines 1->3: */
	WM.Add(1, 3 + 1, Tmp);

	/* Lines 4->6: */
	WM.Sub(3 + 1, 1, Tmp);

	WM.Add(3 + 1, 6 + 1, Tmp);

	/* Lines 7->9: */
	WM.Sub(6 + 1, 3 + 1, Tmp);

	/* [ F x ] [ b2 x ] */
	Tmp = Mat3x3(FTmp, b2);

	/* Lines 4->6: */
	WM.Sub(3 + 1, 9 + 1, Tmp);

	/* Lines 10->12: */
	WM.Add(9 + 1, 9 + 1, Tmp);

	/* [ b1 x ] [ F x ] + [ M x ] */

	/* Lines 4->6: */
	WM.Add(3 + 1, 3 + 1, Mat3x3(b1, FTmp) + Mat3x3(MTmp));

	/* [ b2 x ] [ F x ] + [ M x ] */

	/* Lines 10->12: */
	WM.Sub(9 + 1, 3 + 1, Mat3x3(b2, FTmp) + Mat3x3(MTmp));

/* Phi/q and (Phi/q)^T */

	Mat3x3 b1Cross_R1(Mat3x3(b1)*R1); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(Mat3x3(b2)*R1); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(1, 12 + 1 + iCnt, vR1);
      		WM.Sub(3 + 1, 12 + 1 + iCnt, vb1Cross_R1);

		/* Constraint, node 1 */
      		WM.SubT(12 + 1 + iCnt, 1, vR1);
      		WM.SubT(12 + 1 + iCnt, 3 + 1, vb1Cross_R1);

		/* Equilibrium, node 2 */
      		WM.Add(6 + 1, 12 + 1 + iCnt, vR1);
      		WM.Add(9 + 1, 12 + 1 + iCnt, vb2Cross_R1);

		/* Constraint, node 2 */
      		WM.AddT(12 + 1 + iCnt, 6 + 1, vR1);
      		WM.AddT(12 + 1 + iCnt, 9 + 1, vb2Cross_R1);
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(3 + 1, 12 + 1 + nPosConstraints +  iCnt, vR1);

		/* Constraint, node 1 */
      		WM.SubT(12 + 1 + nPosConstraints + iCnt, 3 + 1, vR1);

		/* Equilibrium, node 2 */
      		WM.Add(9 + 1, 12 + 1 + nPosConstraints + iCnt, vR1);

		/* Constraint, node 2 */
      		WM.AddT(12 + 1 + nPosConstraints +  iCnt, 9 + 1, vR1);
	}

	return WorkMat;
}

/* Assemblaggio residuo */
SubVectorHandler&
TotalJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering TotalJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Get constraint reactions */

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
	}


	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	Mat3x3 R1 = pNode1->GetRCurr()*R1h;
	Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
	Mat3x3 R2r = pNode2->GetRCurr()*R2hr;

	Vec3 XDelta = R1.Transpose()*b1 - tilde_f1 - XDrv.Get();

	Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
	Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
	ThetaDelta = RotManip::VecRot(RDelta);

	Vec3 FTmp(R1*F);
	Vec3 MTmp(R1r*M);

	/* Equilibrium, node 1 */
	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, MTmp + b1.Cross(FTmp));

	/* Equilibrium, node 2 */
	WorkVec.Sub(6 + 1, FTmp);
	WorkVec.Sub(9 + 1, MTmp + b2.Cross(FTmp));

	/* Constraint equations are divided by dCoef */
	if (dCoef != 0.) {

		/* Position constraint:  */
		for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
			WorkVec.PutCoef(12 + 1 + iCnt,
				-(XDelta(iPosIncid[iCnt])/dCoef));
		}

		/* Rotation constraints: */
		for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
			WorkVec.PutCoef(12 + 1 + nPosConstraints + iCnt,
				-(ThetaDelta(iRotIncid[iCnt])/dCoef));
		}
	}

	return WorkVec;
}

/* Inverse Dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
TotalJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	/*
	 * identical to regular AssJac's lower-left block
	 */
	DEBUGCOUT("Entering TotalJoint::AssJac()" << std::endl);

	if (iGetNumDof() == 12) {
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* original - nodes, nodes */
	WM.ResizeReset(iNumRows - 12, 12);

	/* Recupera gli indici delle varie incognite */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	Mat3x3 b1Cross_R1(Mat3x3(b1)*R1); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(Mat3x3(b2)*R1); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(1 + iCnt, 1, vR1);
      		WM.SubT(1 + iCnt, 3 + 1, vb1Cross_R1);

		/* Constraint, node 2 */
      		WM.AddT(1 + iCnt, 6 + 1, vR1);
      		WM.AddT(1 + iCnt, 9 + 1, vb2Cross_R1);
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(1 + nPosConstraints + iCnt, 3 + 1, vR1);

		/* Constraint, node 2 */
      		WM.AddT(1 + nPosConstraints +  iCnt, 9 + 1, vR1);
	}

	return WorkMat;
}

/* Inverse Dynamics residual assembly */
SubVectorHandler&
TotalJoint::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr,
	const VectorHandler& /* XPrimePrimeCurr */,
	int iOrder)
{
	DEBUGCOUT("Entering TotalJoint::AssRes(" << iOrder<< ")" << std::endl);

	if (iGetNumDof() == 0) {
		WorkVec.ResizeReset(0);
		return WorkVec;
	}

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	/* original - node equations (6 * 2) */
	WorkVec.ResizeReset(iNumRows - 12);

	/* Indici */
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstReactionIndex + iCnt);
	}

	switch (iOrder) {
	case 0:	// Position - Orientation
		/*
		 * identical to regular AssRes' lower block
		 */
		{ /* need brackets to create a "block" */	
			Vec3 b2(pNode2->GetRCurr()*f2);
			Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
		
			Mat3x3 R1 = pNode1->GetRCurr()*R1h;
			Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
			Mat3x3 R2r = pNode2->GetRCurr()*R2hr;
		
			Vec3 XDelta = R1.Transpose()*b1 - tilde_f1 - XDrv.Get();

			Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
			Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
			Vec3 ThetaDelta = RotManip::VecRot(RDelta);
		
			/* Position constraint  */
			for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
				WorkVec.DecCoef(1 + iCnt, XDelta(iPosIncid[iCnt]));
			}
	
			/* Rotation constraint  */
			for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
				WorkVec.DecCoef(1 + nPosConstraints + iCnt, ThetaDelta(iRotIncid[iCnt]));
			}
		}// end case 0:
		break;

	case 1:	// Velocity
		/*
		 * first derivative of regular AssRes' lower block
		 */
		{ /* need brackets to create a "block" */
			Vec3 Tmp = XPDrv.Get(); 	
			
			/* Position constraint derivative  */
			for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
				WorkVec.PutCoef(1 + iCnt, Tmp(iPosIncid[iCnt]));
			}

			Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
			Mat3x3 R2r = pNode2->GetRCurr()*R2hr;
		
			Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
			Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
	
			/*This name is only for clarity...*/
			Vec3 WDelta = RDelta * OmegaDrv.Get();
	
			/* Rotation constraint derivative */
			for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
				WorkVec.PutCoef(1 + nPosConstraints + iCnt, WDelta(iRotIncid[iCnt]));
			}
		} // end case 1:	
		break;

	case 2:	// Acceleration
		/*
		 * second derivative of regular AssRes' lower block
		 */
		{ /* need brackets to create a "block" */
			Vec3 b2(pNode2->GetRCurr()*f2);
			Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
		
			Mat3x3 R1 = pNode1->GetRCurr()*R1h;
			
			Vec3 b1Prime(pNode2->GetVCurr() + pNode2->GetWCurr().Cross(b2) - pNode1->GetVCurr());  
			
			Vec3 Tmp2 = (- pNode1->GetWCurr().Cross(b1) + b1Prime).Cross(pNode1->GetWCurr());
			Tmp2 +=  (
				  pNode1->GetWCurr().Cross(- pNode2->GetVCurr() 
							   + b2.Cross(pNode2->GetWCurr()) 
							   + pNode1->GetVCurr())
				  -pNode2->GetWCurr().Cross(b2.Cross(pNode2->GetWCurr()))
				);
			Vec3 Tmp = R1.Transpose()*Tmp2 - XPPDrv.Get();
			
		/* Position constraint second derivative  */
			for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
				WorkVec.DecCoef(1 + iCnt, Tmp(iPosIncid[iCnt]));
			}
	
			Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
			Mat3x3 R2r = pNode2->GetRCurr()*R2hr;
			Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
			Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
		
			Tmp = R0T * OmegaDrv.Get();
			Tmp2 = R2r * Tmp;
			Tmp = Mat3x3(pNode2->GetWCurr() - pNode1->GetWCurr()) * Tmp;
			Tmp += Mat3x3(pNode1->GetWCurr())*pNode2->GetWCurr();
			Tmp2 = R1r.Transpose()*Tmp;
			Tmp2 += RDelta * OmegaPDrv.Get();
			
			/* Rotation constraint second derivative */
			for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
				WorkVec.PutCoef(1 + nPosConstraints + iCnt, Tmp2(iRotIncid[iCnt]));
			}
		} // end case 2:
		break;
	
	default:
		/*
		 * higher-order derivatives make no sense right now
		 */

		ASSERT(0);
		throw ErrGeneric();
	}

	return WorkVec;
}

/* Inverse Dynamics update */
void 
TotalJoint::Update(const VectorHandler& XCurr, int /*iOrder*/)
{
	integer iFirstReactionIndex = iGetFirstIndex();
	
	/* Get constraint reactions */

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
			F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
			M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
	}
};

DofOrder::Order
TotalJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(),
		"INDEX ERROR in TotalJoint::GetEqType");

	return DofOrder::ALGEBRAIC;
}

/* Output (da mettere a punto), per ora solo reazioni */
void
TotalJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 R1Tmp(pNode1->GetRCurr()*R1h);
		Mat3x3 R1rTmp(pNode1->GetRCurr()*R1hr);
		Mat3x3 R1rTmpT(R1rTmp.Transpose());
		Mat3x3 R2rTmp(pNode2->GetRCurr()*R2hr);
		Mat3x3 RTmp(R1rTmpT*R2rTmp);

		Joint::Output(OH.Joints(), "TotalJoint", GetLabel(),
			F, M, R1Tmp*F, R1rTmp*M)
			<< " " << MatR2EulerAngles(RTmp)*dRaDegr
			<< " " << R1rTmpT*(pNode2->GetWCurr() - pNode1->GetWCurr())
			<< " " << ThetaDeltaPrev << std::endl;
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
TotalJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{

	/* Per ora usa la matrice piena; eventualmente si puo'
	 * passare a quella sparsa quando si ottimizza */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Equazioni: vedi joints.dvi */

	/*	 equazioni ed incognite
	 * F1					  Delta_x1	1	   
	 * M1					  Delta_g1	3 + 1 
	 * FP1  				  Delta_xP1	6 + 1
	 * MP1  				  Delta_w1	9 + 1 
	 * F2					  Delta_x2	12 + 1
	 * M2					  Delta_g2	15 + 1
	 * FP2  				  Delta_xP2	18 + 1  
	 * MP2  				  Delta_w2	21 + 1 
	 * vincolo spostamento  		  Delta_F	24 + 1 
	 * vincolo rotazione			  Delta_M	24 + nPosConstraints  
	 * derivata vincolo spostamento 	  Delta_FP	24 + nConstraints  
	 * derivata vincolo rotazione		  Delta_MP	24 + nConstraints + nPosConstraints  
	 */
	
	
	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;


	/* Setta gli indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	/* Setta gli indici delle reazioni */
	for (unsigned int iCnt = 1; iCnt <=  nConstraints; iCnt++) {
		WM.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutRowIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
		WM.PutColIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}
	
	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	
	/* F ed M sono gia' state aggiornate da InitialAssRes;
	 * Recupero FPrime e MPrime*/
	Vec3 MPrime(0.);
	Vec3 FPrime(0.);

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		FPrime(iPosIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + iCnt);
	}
	
	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		MPrime(iRotIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + nPosConstraints + iCnt);
	}

	Vec3 FTmp(R1 * F);
	Vec3 MTmp(R1r * M);
	Vec3 FPrimeTmp(R1 * FPrime);
	Vec3 MPrimeTmp(R1r * MPrime);

	Mat3x3 Tmp;

	/* [ F x ] */
	Tmp = Mat3x3(FTmp);

	/* Force, Node 1, Lines 1->3: */
	WM.Add(1, 3 + 1, Tmp);	// * Delta_g1

	/* Moment, Node 1, Lines 4->6: */
	WM.Sub(3 + 1, 1, Tmp);	// * Delta_x1

	WM.Add(3 + 1, 12 + 1, Tmp);	// * Delta_x2

	/* Force, Node 2, Lines 13->15: */
	WM.Sub(12 + 1, 3 + 1, Tmp);	// * Delta_g1
	
	/* [ FP x ] */
	Tmp = Mat3x3(FPrimeTmp);

	/* d/dt(Force), Node 1, Lines 7->9: */
	WM.Add(6 + 1, 9 + 1, Tmp);	// * Delta_W1

	/* d/dt(Moment), Node 1, Lines 10->12: */
	WM.Sub(9 + 1, 6 + 1, Tmp);	// * Delta_x1P

	WM.Add(9 + 1, 18 + 1 + 6, Tmp);	// * Delta_x2P

	/* d/dt(Force), Node 2, Lines 19->21: */
	WM.Sub(18 + 1 , 9 + 1, Tmp);	// * Delta_W1
	
	/* [ F x ] [ b2 x ] */
	Tmp = Mat3x3(FTmp, b2);

	/* Moment, Node1, Lines 4->6: */
	WM.Sub(3 + 1, 15 + 1, Tmp);	// * Delta_g2

	/* Moment, Node2, Lines 16->18: */
	WM.Add(15 + 1, 15 + 1, Tmp);	// * Delta_g2
	
	/* [ FP x ] [ b2 x ] */
	Tmp = Mat3x3(FPrimeTmp, b2);

	/* d/dt(Moment), Node1, Lines 10->12: */
	WM.Sub(9 + 1, 21 + 1, Tmp);	// * Delta_W2

	/* d/dt(Moment), Node2, Lines 22->24: */
	WM.Add(21 + 1, 21 + 1, Tmp);	// * Delta_W2

	/* [ b1 x ] [ F x ] + [ M x ] */

	/* Moment, Node1, Lines 4->6: */
	WM.Add(3 + 1, 3 + 1, Mat3x3(b1, FTmp) + Mat3x3(MTmp));	// *Delta_g1

	/* d/dt(Moment), Node1, Lines 10->12: */
	WM.Add(9 + 1, 9 + 1, Mat3x3(b1, FPrimeTmp) + Mat3x3(MPrimeTmp));// *Delta_W1
	
	/* [ b2 x ] [ F x ] + [ M x ] */

	/* Moment, Node2, Lines 16->18: */
	WM.Sub(15 + 1, 3 + 1, Mat3x3(b2, FTmp) + Mat3x3(MTmp));	// * Delta_g1
	
	/* d/dt(Moment), Node2, Lines 22->24: */
	WM.Sub(21 + 1, 9 + 1, Mat3x3(b2, FTmp) + Mat3x3(MTmp));	// * Delta_W1
	
	/* Constraints: Add only active rows/columns*/	
	
	/* Positions contribution:*/

	Mat3x3 b1Cross_R1(Mat3x3(b1)*R1); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(Mat3x3(b2)*R1); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(1, 24 + 1 + iCnt, vR1);	// * Delta_F
      		WM.Sub(3 + 1, 24 + 1 + iCnt, vb1Cross_R1);	// * Delta_F

		/* Constraint, node 1 */
      		WM.SubT(24 + 1 + iCnt, 1, vR1);	// * Delta_x1
      		WM.SubT(24 + 1 + iCnt, 3 + 1, vb1Cross_R1);	// * Delta_g1

		/* d/dt(Equilibrium), node 1 */
      		WM.Sub(6 + 1, 24 + 1 + nConstraints + iCnt, vR1);	// * Delta_FP
      		WM.Sub(9 + 1, 24 + 1 + nConstraints + iCnt, vb1Cross_R1);	// * Delta_FP

		/* d/dt(Constraint), node 1 */
      		WM.SubT(24 + 1 + nConstraints + iCnt, 6 + 1, vR1);	// * Delta_xP1
      		WM.SubT(24 + 1 + nConstraints + iCnt, 9 + 1, vb1Cross_R1);	// * Delta_W1
		
		/* Equilibrium, node 2 */
      		WM.Add(12 + 1, 24 + 1 + iCnt, vR1);	// * Delta_F
      		WM.Add(15 + 1, 24 + 1 + iCnt, vb2Cross_R1);	// * Delta_F

		/* Constraint, node 2 */
      		WM.AddT(24 + 1 + iCnt, 12 + 1, vR1);	// * Delta_x2
      		WM.AddT(24 + 1 + iCnt, 15 + 1, vb2Cross_R1);	// * Delta_g2
	
		/* d/dt(Equilibrium), node 2 */
      		WM.Add(18 + 1, 24 + 1 + nConstraints + iCnt, vR1);	// * Delta_FP
      		WM.Add(21 + 1, 24 + 1 + nConstraints + iCnt, vb2Cross_R1);	// * Delta_FP

		/* d/dt(Constraint), node 2 */
      		WM.AddT(24 + 1 + nConstraints +  iCnt, 18 + 1, vR1);	// * Delta_xP2
      		WM.AddT(24 + 1 + nConstraints +  iCnt, 21 + 1, vb2Cross_R1);	// * Delta_W2
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Equilibrium, node 1 */
      		WM.Sub(3 + 1, 24 + 1 + nPosConstraints +  iCnt, vR1);	// * Delta_M

		/* Constraint, node 1 */
      		WM.SubT(24 + 1 + nPosConstraints + iCnt, 3 + 1, vR1);	// * Delta_g1

		/* d/dt(Equilibrium), node 1 */
      		WM.Sub(9 + 1, 24 + 1 + nConstraints + nPosConstraints +  iCnt, vR1);	// * Delta_MP

		/* d/dt(Constraint), node 1 */
      		WM.SubT(24 + 1 + nConstraints + nPosConstraints + iCnt, 9 + 1, vR1);	// * Delta_W1

		/* Equilibrium, node 2 */
      		WM.Add(15 + 1, 24 + 1 + nPosConstraints + iCnt, vR1);	// * Delta_M

		/* Constraint, node 2 */
      		WM.AddT(24 + 1 + nPosConstraints +  iCnt, 15 + 1, vR1);	// * Delta_g2

		/* d/dt(Equilibrium), node 2 */
      		WM.Add(21 + 1, 24 + 1 + nConstraints + nPosConstraints + iCnt, vR1);	// * Delta_MP

		/* d/dt(Constraint), node 2 */
      		WM.AddT(24 + 1 + nConstraints + nPosConstraints + iCnt, 21 + 1, vR1);	// * Delta_W2
	}

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
TotalJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{

	DEBUGCOUT("Entering TotalJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++)	{
		WorkVec.PutRowIndex(24 + iCnt, iFirstReactionIndex + iCnt);
		WorkVec.PutRowIndex(24 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}

	/* Recupera i dati */
	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Mat3x3 R2r(pNode2->GetRCurr()*R2hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Mat3x3 Omega2Cross(Omega2);
	Vec3 Omega2Crossb2(Omega2Cross*b2);
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	
	Vec3 FPrime(0.);
	Vec3 MPrime(0.);

	/* Aggiorna F ed M, che restano anche per InitialAssJac */
	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
		FPrime(iPosIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
		MPrime(iRotIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + nPosConstraints + iCnt);
	}

	Vec3 FTmp(R1 * F);
	Vec3 MTmp(R1r * M);
	Vec3 FPrimeTmp(R1 * FPrime);
	Vec3 MPrimeTmp(R1r * MPrime);

	/* Equilibrium, node 1 */
	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, b1.Cross(FTmp) + MTmp);

	/* d/dt(Equilibrium), node 1 */
	WorkVec.Add(6 + 1, FPrimeTmp);
	WorkVec.Add(9 + 1, b1.Cross(FPrimeTmp) + MPrimeTmp);

	/* Equilibrium, node 2 */
	WorkVec.Sub(12 + 1, FTmp);
	WorkVec.Sub(15 + 1, b2.Cross(FTmp) + MTmp);

	/* d/dt( Equilibrium ) , node 2 */
	WorkVec.Sub(18 + 1, FPrimeTmp);
	WorkVec.Sub(21 + 1, b2.Cross(FPrimeTmp) + MPrimeTmp);

	/* Constraint Equations */
	
	Vec3 XDelta = R1.Transpose()*b1 - tilde_f1 - XDrv.Get();
	Vec3 XDeltaPrime = R1.Transpose()*(b1Prime + b1.Cross(Omega1));
	
	if(XDrv.bIsDifferentiable())	{
		XDeltaPrime -= XDrv.GetP();
	}
	
	Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
	Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
	ThetaDelta = RotManip::VecRot(RDelta);
	Vec3 ThetaDeltaPrime = R1r.Transpose()*(Omega2 - Omega1);

	if(ThetaDrv.bIsDifferentiable())	{
		ThetaDeltaPrime -= RDelta*ThetaDrv.GetP();
	}


	/* Position constraint:  */
	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		WorkVec.PutCoef(24 + 1 + iCnt, -XDelta(iPosIncid[iCnt]));
		WorkVec.PutCoef(24 + 1 + nConstraints + iCnt, -XDeltaPrime(iPosIncid[iCnt]));
	}

	/* Rotation constraints: */
	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		WorkVec.PutCoef(24 + 1 + nPosConstraints + iCnt, -ThetaDelta(iRotIncid[iCnt]));
		WorkVec.PutCoef(24 + 1 + nPosConstraints + nConstraints +  iCnt, -ThetaDeltaPrime(iRotIncid[iCnt]));
	}

return WorkVec;
}


unsigned int
TotalJoint::iGetNumPrivData(void) const
{
	return 18;
}

unsigned int
TotalJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned int off = 0;

	switch (s[0]) {
	case 'p':
		/* relative position */
		break;

	case 'r':
		/* relative orientation */
		off += 3;
		break;

	case 'F':
		/* force */
		off += 6;
		break;

	case 'M':
		/* moment */
		off += 9;
		break;

	case 'd':
		/* imposed relative position */
		off += 12;
		break;

	case 't':
		/* imposed relative orientation */
		off += 15;
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		return off + 1;

	case 'y':
		return off + 2;

	case 'z':
		return off + 3;
	}

	return 0;
}

doublereal
TotalJoint::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
	case 2:
	case 3: {
		Vec3 x(pNode1->GetRCurr().Transpose()*(
			pNode2->GetXCurr() + pNode2->GetRCurr()*f2
				- pNode1->GetXCurr()) - f1);
			return R1h.GetVec(i)*x;
		}

	case 4:
	case 5:
	case 6: {
		Vec3 Theta(Unwrap(ThetaDeltaPrev, ThetaDelta) + ThetaDrv.Get());
		return Theta(i - 3);
		}

	case 7:
	case 8:
	case 9:
		return F(i - 6);

	case 10:
	case 11:
	case 12:
		return M(i - 9);

	case 13:
	case 14:
	case 15:
		if (!bPosActive[i - 13]) {
			return 0.;
		}
		return XDrv.Get()(i - 12);

	case 16:
	case 17:
	case 18:
		if (!bRotActive[i - 16]) {
			return 0.;
		}
		return ThetaDrv.Get()(i - 15);

	default:
		ASSERT(0);
	}

	return 0.;
}

/* TotalJoint - end */

/* TotalPinJoint - begin */

TotalPinJoint::TotalPinJoint(unsigned int uL, const DofOwner *pDO,
	bool bPos[3],
	TplDriveCaller<Vec3> *const pDCPos[3],
	bool bRot[3],
	TplDriveCaller<Vec3> *const pDCRot[3],
	const Vec3& XcTmp, const Mat3x3& RchTmp, const Mat3x3& RchrTmp,
	const StructNode *pN,
	const Vec3& fnTmp, const Mat3x3& RnhTmp, const Mat3x3& RnhrTmp,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
pNode(pN),
Xc(XcTmp), Rch(RchTmp), Rchr(RchrTmp),
tilde_fn(fnTmp), tilde_Rnh(RnhTmp), tilde_Rnhr(RnhrTmp),
RchT(Rch.Transpose()), tilde_Xc(RchT*Xc), RchrT(Rchr.Transpose()),
XDrv(pDCPos[0]), XPDrv(pDCPos[1]), XPPDrv(pDCPos[2]),
ThetaDrv(pDCRot[0]), OmegaDrv(pDCRot[1]), OmegaPDrv(pDCRot[2]),
nConstraints(0), nPosConstraints(0), nRotConstraints(0),
M(0.), F(0.), ThetaDelta(0.), ThetaDeltaPrev(0.)
{
	/* Equations 1->3: Positions
	 * Equations 4->6: Rotations */

	for (unsigned int i = 0; i < 3; i++) {
		bPosActive[i] = bPos[i];
		bRotActive[i] = bRot[i];
		if (bPosActive[i]) {
			iPosIncid[nPosConstraints] = i + 1;
			nPosConstraints++;
		}
		if (bRotActive[i]) {
			iRotIncid[nRotConstraints] = i + 1;
			nRotConstraints++;
		}
	}
	nConstraints = nPosConstraints + nRotConstraints;
}

TotalPinJoint::~TotalPinJoint(void)
{
	NO_OP;
};

std::ostream&
TotalPinJoint::DescribeDof(std::ostream& out,
	char *prefix, bool bInitial, int i) const
{
	integer iIndex = iGetFirstIndex();

	if (i >= 0) {
		silent_cerr("TotalPinJoint(" << GetLabel() << "): "
			"DescribeDof(" << i << ") "
			"not implemented yet" << std::endl);
		throw ErrGeneric();
	}

	if (nPosConstraints > 1) {
		out << prefix << iIndex + 1;
		out << "->" << iIndex + nPosConstraints;
		out << ": ";
	}
	out << "reaction force(s) [";

	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bPosActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "F" << idx2xyz[i];
		}
	}
	out << "]" << std::endl;
	

	if (nRotConstraints > 1) {
		out << prefix << iIndex + nPosConstraints + 1;
		out << "->" << iIndex + nConstraints ;
		out << ": ";
	}
	out << "reaction couple(s) [";

	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bRotActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "m" << idx2xyz[i];
		}
	}
	out << "]" << std::endl;

	if (bInitial) {
		iIndex += nConstraints;

		if (nPosConstraints > 1) {
			out << prefix << iIndex + 1;
			out << "->" << iIndex + nPosConstraints;
			out << ": ";
		}
		out << "reaction force(s) derivative(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bPosActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "FP" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
		

		if (nRotConstraints > 1) {
			out << prefix << iIndex + nPosConstraints + 1;
			out << "->" << iIndex + nConstraints;
			out << ": ";
		}
		out << "reaction couple(s) derivative(s) [";

		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bRotActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "mP" << idx2xyz[i];
			}
		}
		out << "]" << std::endl;
	}
	return out;
}

std::ostream&
TotalPinJoint::DescribeEq(std::ostream& out,
	char *prefix, bool bInitial, int i) const
{
	integer iIndex = iGetFirstIndex();

	if (i >= 0) {
		silent_cerr("TotalPinJoint(" << GetLabel() << "): "
			"DescribeEq(" << i << ") "
			"not implemented yet" << std::endl);
		throw ErrGeneric();
	}

	if (nPosConstraints > 1) {
		out << prefix << iIndex + 1;
		out << "->" << iIndex + nPosConstraints;
		out << ": ";
	}
	
	if (nPosConstraints > 0) {
		out << "position constraint(s) [";
	}
	
	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bPosActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "P" << idx2xyz[i] << "1=P" << idx2xyz[i] << "2";
		}
	}
	
	if (nPosConstraints > 0) {
		out << "]" << std::endl;
	}

	if (nRotConstraints > 1) {
		out << prefix << iIndex + nPosConstraints + 1;
		out << "->" << iIndex + nConstraints ;
		out << ": ";
	}
	if (nRotConstraints > 0) {
		out << "orientation constraint(s) [";
	}
	
	for (unsigned int i = 0, cnt = 0; i < 3; i++) {
		if (bRotActive[i]) {
			cnt++;
			if (cnt > 1) {
				out << ",";
			}
			out << "g" << idx2xyz[i] << "1=g" << idx2xyz[i] << "2";
		}
	}
	if (nRotConstraints > 0) {
		out << "]" << std::endl;
	}

	if (bInitial) {
		iIndex += nConstraints;

		if (nPosConstraints > 1) {
			out << prefix << iIndex + 1;
			out << "->" << iIndex + nPosConstraints;
			out << ": ";
		}
		
		if (nPosConstraints > 0) {
			out << "velocity constraint(s) [";
		}
		
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bPosActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "v" << idx2xyz[i] << "1=v" << idx2xyz[i] << "2";
			}
		}
		if (nPosConstraints > 0) {
			out << "]" << std::endl;
		}
		
		if (nRotConstraints > 1) {
			out << prefix << iIndex + nPosConstraints + 1;
			out << "->" << iIndex + nConstraints ;
			out << ": ";
		}
		if (nRotConstraints > 0) {
			out << "angular velocity constraint(s) [";
		}
		
		for (unsigned int i = 0, cnt = 0; i < 3; i++) {
			if (bRotActive[i]) {
				cnt++;
				if (cnt > 1) {
					out << ",";
				}
				out << "w" << idx2xyz[i] << "1=w" << idx2xyz[i] << "2";
			}
		}
		if (nRotConstraints > 0) {
			out << "]" << std::endl;
		}	
	}

	return out;
}

void
TotalPinJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned int i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh) {
				if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
					Mat3x3 RnT(pNode->GetRCurr().Transpose());

					tilde_fn = RnT*(Xc + Rch*XDrv.Get() - pNode->GetXCurr());

				} else if (dynamic_cast<Joint::OffsetHint<0> *>(pjh)) {
					Xc = pNode->GetXCurr() + pNode->GetRCurr()*tilde_fn
						- Rch*XDrv.Get();
					tilde_Xc = RchT*Xc;

				} else if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<1> *>(pjh)) {
						tilde_Rnh = pNode->GetRCurr().Transpose()*Rch;
						/* NOTE: pointless */

					} else if (dynamic_cast<Joint::OrientationHingeHint<1> *>(pjh)) {
						tilde_Rnhr = pNode->GetRCurr().Transpose()*Rchr*RotManip::Rot(ThetaDrv.Get());
					}

				} else if (dynamic_cast<Joint::HingeHint<0> *>(pjh)) {
					if (dynamic_cast<Joint::PositionHingeHint<0> *>(pjh)) {
						Rch = pNode->GetRCurr()*tilde_Rnh;
						RchT = Rch.Transpose();
						tilde_Xc = RchT*Xc;
						/* NOTE: results in constraint violation if XDrv is not 0 */

					} else if (dynamic_cast<Joint::OrientationHingeHint<0> *>(pjh)) {
						Rchr = pNode->GetRCurr()*tilde_Rnhr*RotManip::Rot(-ThetaDrv.Get());
						RchrT = Rchr.Transpose();
					}

				} else if (dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh)) {
					Joint::JointDriveHint<Vec3> *pjdh
						= dynamic_cast<Joint::JointDriveHint<Vec3> *>(pjh);
					pedantic_cout("TotalPinJoint(" << uLabel << "): "
						"creating drive from hint[" << i << "]..." << std::endl);

					TplDriveCaller<Vec3> *pDC = pjdh->pTDH->pCreateDrive(pDM);
					if (pDC == 0) {
						silent_cerr("TotalPinJoint(" << uLabel << "): "
							"unable to create drive "
							"after hint #" << i << std::endl);
						throw ErrGeneric();
					}

					if (dynamic_cast<Joint::PositionDriveHint<Vec3> *>(pjdh)) {
						XDrv.Set(pDC);

					} else if (dynamic_cast<Joint::VelocityDriveHint<Vec3> *>(pjdh)) {
						XPDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AccelerationDriveHint<Vec3> *>(pjdh)) {
						XPPDrv.Set(pDC);

					} else if (dynamic_cast<Joint::OrientationDriveHint<Vec3> *>(pjdh)) {
						ThetaDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AngularVelocityDriveHint<Vec3> *>(pjdh)) {
						OmegaDrv.Set(pDC);

					} else if (dynamic_cast<Joint::AngularAccelerationDriveHint<Vec3> *>(pjdh)) {
						OmegaPDrv.Set(pDC);

					} else {
						delete pDC;
					}

				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
					/* TODO */
				}
				continue;
			}
		}
	}
}

Hint *
TotalPinJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/ , STRLENOF("offset{" /*}*/ )) == 0)
	{
		s += STRLENOF("offset{" /*}*/ );

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '0':
			return new Joint::OffsetHint<0>;
		}

	} else if (strncasecmp(s, "position-hinge{" /*}*/, STRLENOF("position-hinge{" /*}*/)) == 0) {
		s += STRLENOF("position-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '0':
			return new Joint::HingeHint<0>;
		}

	} else if (strncasecmp(s, "position-drive3{" /*}*/, STRLENOF("position-drive3{" /*}*/)) == 0) {
		s += STRLENOF("position-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new PositionDriveHint<Vec3>(pTDH);
			}
		}
		return 0;

	} else if (strncasecmp(s, "orientation-hinge{" /*}*/, STRLENOF("orientation-hinge{" /*}*/)) == 0) {
		s += STRLENOF("orientation-hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '0':
			return new Joint::HingeHint<0>;
		}

	} else if (strncasecmp(s, "orientation-drive3{" /*}*/, STRLENOF("orientation-drive3{" /*}*/)) == 0) {
		s += STRLENOF("orientation-");

		Hint *pH = ::ParseHint(pDM, s);
		if (pH) {
			TplDriveHint<Vec3> *pTDH = dynamic_cast<TplDriveHint<Vec3> *>(pH);
			if (pTDH) {
				return new OrientationDriveHint<Vec3>(pTDH);
			}
		}
		return 0;
	}

	return 0;
}

void
TotalPinJoint::AfterConvergence(const VectorHandler& /* X */ ,
	const VectorHandler& /* XP */ )
{
	ThetaDeltaPrev = Unwrap(ThetaDeltaPrev, ThetaDelta);
}

/* Contributo al file di restart */
std::ostream&
TotalPinJoint::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", total joint, "
		<< pNode->GetLabel() << ", "
			<< "position, " << tilde_fn.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", tilde_Rnh.GetVec(1).Write(out, ", ") << ", "
				"2 , ", tilde_Rnh.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", tilde_Rnhr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", tilde_Rnhr.GetVec(2).Write(out, ", ") << ", "
			<< "position, " << Xc.Write(out, ", ") << ", "
			<< "position orientation, "
				"1 , ", Rch.GetVec(1).Write(out, ", ") << ", "
				"2 , ", Rch.GetVec(2).Write(out, ", ") << ", "
			<< "rotation orientation, "
				"1 , ", Rchr.GetVec(1).Write(out, ", ") << ", "
				"2 , ", Rchr.GetVec(2).Write(out, ", ");

	if (bPosActive[0] || bPosActive[1] || bPosActive[2]) {
		out << ", position constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bPosActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}

		/* TODO: output drive(s) */
		out << ", zero";
	}

	if (bRotActive[0] || bRotActive[1] || bRotActive[2]) {
		out << ", orientation constraint";
		for (unsigned i = 0; i < 3; i++) {
			if (bRotActive[i]) {
				out << ", active";
			} else {
				out << ", inactive";
			}
		}

		/* TODO: output drive(s) */
		out << ", zero";
	}

	return out << ";" << std::endl;
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
TotalPinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	/*
	 * See tecman.pdf
	 */

	DEBUGCOUT("Entering TotalPinJoint::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WM.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Vec3 fn(pNode->GetRCurr()*tilde_fn);

	/* Moltiplica la forza per il coefficiente del metodo */
	Vec3 FTmp(Rch*(F*dCoef));

	/* Equilibrium: ((Phi/q)^T*Lambda)/q */

	/* Lines 4->6: */
	/* [ F x ] [ b2 x ] */
	WM.Add(3 + 1, 3 + 1, Mat3x3(FTmp, fn));

/* Phi/q and (Phi/q)^T */

	Mat3x3 fnCross_Rch(Mat3x3(fn)*Rch); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vRch(Rch.GetVec(iPosIncid[iCnt]));
		Vec3 vfnCross_Rch(fnCross_Rch.GetVec(iPosIncid[iCnt]));

		/* Equilibrium, node */
      		WM.Add(1, 6 + 1 + iCnt, vRch);
      		WM.Add(3 + 1, 6 + 1 + iCnt, vfnCross_Rch);

		/* Constraint, node */
      		WM.AddT(6 + 1 + iCnt, 1, vRch);
      		WM.AddT(6 + 1 + iCnt, 3 + 1, vfnCross_Rch);
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vRchr(Rchr.GetVec(iRotIncid[iCnt]));

		/* Equilibrium, node */
      		WM.Add(3 + 1, 6 + 1 + nPosConstraints + iCnt, vRchr);

		/* Constraint, node */
      		WM.AddT(6 + 1 + nPosConstraints + iCnt, 3 + 1, vRchr);
	}

	return WorkMat;
}

/* Assemblaggio residuo */
SubVectorHandler&
TotalPinJoint::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering TotalPinJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
	}

	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WorkVec.PutRowIndex(6 + iCnt, iFirstReactionIndex + iCnt);
	}

	/* Get constraint reactions */

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
	}

	Vec3 fn(pNode->GetRCurr()*tilde_fn);

	Mat3x3 Rnhr = pNode->GetRCurr()*tilde_Rnhr;

	Vec3 XDelta = RchT*(pNode->GetXCurr() + fn) - tilde_Xc - XDrv.Get();

	Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
	Mat3x3 RDelta = RchrT*Rnhr*R0T;
	ThetaDelta = RotManip::VecRot(RDelta);

	Vec3 FTmp(Rch*F);
	Vec3 MTmp(Rchr*M);

	/* Equilibrium, node 2 */
	WorkVec.Sub(1, FTmp);
	WorkVec.Sub(3 + 1, MTmp + fn.Cross(FTmp));

	/* Constraint equations are divided by dCoef */
	if (dCoef != 0.) {

		/* Position constraint:  */
		for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
			WorkVec.PutCoef(6 + 1 + iCnt, -(XDelta(iPosIncid[iCnt])/dCoef));
		}

		/* Rotation constraints: */
		for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
			WorkVec.PutCoef(6 + 1 + nPosConstraints + iCnt, -(ThetaDelta(iRotIncid[iCnt])/dCoef));
		}
	}

	return WorkVec;
}

#if 0
/* Inverse Dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
TotalPinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	/*
	 * identical to regular AssJac's lower-left block
	 */
	DEBUGCOUT("Entering TotalPinJoint::AssJac()" << std::endl);

	if (iGetNumDof() == 12) {
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* original - nodes, nodes */
	WM.ResizeReset(iNumRows - 12, 12);

	/* Recupera gli indici delle varie incognite */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WM.PutRowIndex(iCnt, iFirstReactionIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	Mat3x3 b1Cross_R1(Mat3x3(b1)*R1); // = [ b1 x ] * R1
	Mat3x3 b2Cross_R1(Mat3x3(b2)*R1); // = [ b2 x ] * R1

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vR1(R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb1Cross_R1(b1Cross_R1.GetVec(iPosIncid[iCnt]));
		Vec3 vb2Cross_R1(b2Cross_R1.GetVec(iPosIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(1 + iCnt, 1, vR1);
      		WM.SubT(1 + iCnt, 3 + 1, vb1Cross_R1);

		/* Constraint, node 2 */
      		WM.AddT(1 + iCnt, 6 + 1, vR1);
      		WM.AddT(1 + iCnt, 9 + 1, vb2Cross_R1);
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vR1(R1r.GetVec(iRotIncid[iCnt]));

		/* Constraint, node 1 */
      		WM.SubT(1 + nPosConstraints + iCnt, 3 + 1, vR1);

		/* Constraint, node 2 */
      		WM.AddT(1 + nPosConstraints +  iCnt, 9 + 1, vR1);
	}

	return WorkMat;
}

/* Inverse Dynamics residual assembly */
SubVectorHandler&
TotalPinJoint::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr,
	const VectorHandler& /* XPrimePrimeCurr */,
	int iOrder)
{
	DEBUGCOUT("Entering TotalPinJoint::AssRes(" << iOrder<< ")" << std::endl);

	if (iGetNumDof() == 0) {
		WorkVec.ResizeReset(0);
		return WorkVec;
	}

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);

	/* original - node equations (6 * 2) */
	WorkVec.ResizeReset(iNumRows - 12);

	/* Indici */
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Indici del vincolo */
	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstReactionIndex + iCnt);
	}

	switch (iOrder) {
	case 0:	// Position - Orientation
		/*
		 * identical to regular AssRes' lower block
		 */
		{ /* need brackets to create a "block" */	
			Vec3 b2(pNode2->GetRCurr()*f2);
			Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
		
			Mat3x3 R1 = pNode1->GetRCurr()*R1h;
			Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
			Mat3x3 R2r = pNode2->GetRCurr()*R2hr;
		
			Vec3 XDelta = R1.Transpose()*b1 - tilde_f1 - XDrv.Get();
		
			Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
			Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
			Vec3 ThetaDelta = RotManip::VecRot(RDelta);
		
			/* Position constraint  */
			for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
				WorkVec.DecCoef(1 + iCnt, XDelta(iPosIncid[iCnt]));
			}
	
			/* Rotation constraint  */
			for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
				WorkVec.DecCoef(1 + nPosConstraints + iCnt, ThetaDelta(iRotIncid[iCnt]));
			}
		}// end case 0:
		break;

	case 1:	// Velocity
		/*
		 * first derivative of regular AssRes' lower block
		 */
		{ /* need brackets to create a "block" */
			Vec3 Tmp = XPDrv.Get(); 	
			
			/* Position constraint derivative  */
			for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
				WorkVec.PutCoef(1 + iCnt, Tmp(iPosIncid[iCnt]));
			}

			Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
			Mat3x3 R2r = pNode2->GetRCurr()*R2hr;
		
			Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
			Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
	
			/*This name is only for clarity...*/
			Vec3 WDelta = RDelta * OmegaDrv.Get();
	
			/* Rotation constraint derivative */
			for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
				WorkVec.PutCoef(1 + nPosConstraints + iCnt, WDelta(iRotIncid[iCnt]));
			}
		} // end case 1:	
		break;

	case 2:	// Acceleration
		/*
		 * second derivative of regular AssRes' lower block
		 */
		{ /* need brackets to create a "block" */
			Vec3 b2(pNode2->GetRCurr()*f2);
			Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
		
			Mat3x3 R1 = pNode1->GetRCurr()*R1h;
			
			Vec3 b1Prime(pNode2->GetVCurr() + pNode2->GetWCurr().Cross(b2) - pNode1->GetVCurr());  
			Vec3 Tmp = 	- b1.Cross(pNode1->GetWCurr())
					- b1Prime
					- pNode2->GetVCurr()
					+ b2.Cross(pNode2->GetWCurr())
					+ pNode1->GetVCurr();
			
			Vec3 Tmp2 = Mat3x3(pNode1->GetWCurr())*Tmp;
	
			Tmp2 -= Mat3x3(pNode2->GetWCurr(), pNode2->GetWCurr())*b2;
			Tmp = R1.Transpose()*Tmp2 - XPPDrv.Get();
			
			/* Position constraint second derivative  */
			for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
				WorkVec.DecCoef(1 + iCnt, Tmp(iPosIncid[iCnt]));
			}
	
			Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
			Mat3x3 R2r = pNode2->GetRCurr()*R2hr;
			Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
			Mat3x3 RDelta = R1r.Transpose()*R2r*R0T;
		
			Tmp = R0T * OmegaDrv.Get();
			Tmp2 = R2r * Tmp;
			Tmp = Mat3x3(pNode2->GetWCurr() - pNode1->GetWCurr()) * Tmp;
			Tmp += Mat3x3(pNode1->GetWCurr())*pNode2->GetWCurr();
			Tmp2 = R1r.Transpose()*Tmp;
			Tmp2 += RDelta * OmegaPDrv.Get();
			
			/* Rotation constraint second derivative */
			for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
				WorkVec.PutCoef(1 + nPosConstraints + iCnt, Tmp2(iRotIncid[iCnt]));
			}
		} // end case 2:
		break;
	
	default:
		/*
		 * higher-order derivatives make no sense right now
		 */

		ASSERT(0);
		throw ErrGeneric();
	}

	return WorkVec;
}

/* Inverse Dynamics update */
void 
TotalPinJoint::Update(const VectorHandler& XCurr, int iOrder)
{
	integer iFirstReactionIndex = iGetFirstIndex();
	
	/* Get constraint reactions */

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
			F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
			M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
	}
};
#endif

DofOrder::Order
TotalPinJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(),
		"INDEX ERROR in TotalPinJoint::GetEqType");

	return DofOrder::ALGEBRAIC;
}

/* Output (da mettere a punto) */
void
TotalPinJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 RnhrTmp(pNode->GetRCurr()*tilde_Rnhr);
		Mat3x3 RTmp(RchrT*RnhrTmp);

		Joint::Output(OH.Joints(), "TotalPinJoint", GetLabel(),
			F, M, Rch*F, Rchr*M)
			<< " " << MatR2EulerAngles(RTmp)*dRaDegr
			<< " " << RchrT*pNode->GetWCurr()
			<< " " << ThetaDeltaPrev << std::endl;
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
TotalPinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{

	/* Per ora usa la matrice piena; eventualmente si puo'
	 * passare a quella sparsa quando si ottimizza */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Equazioni: vedi joints.dvi */

	/*	 equazioni ed incognite
	 * F1					  Delta_x1	1	   
	 * M1					  Delta_g1	3 + 1 
	 * FP1  				  Delta_xP1	6 + 1
	 * MP1  				  Delta_w1	9 + 1 
	 * F2					  Delta_x2	12 + 1
	 * M2					  Delta_g2	15 + 1
	 * FP2  				  Delta_xP2	18 + 1  
	 * MP2  				  Delta_w2	21 + 1 
	 * vincolo spostamento  		  Delta_F	24 + 1 
	 * vincolo rotazione			  Delta_M	24 + nPosConstraints  
	 * derivata vincolo spostamento 	  Delta_FP	24 + nConstraints  
	 * derivata vincolo rotazione		  Delta_MP	24 + nConstraints + nPosConstraints  
	 */
	
	
	/* Indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;


	/* Setta gli indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
	}

	/* Setta gli indici delle reazioni */
	for (unsigned int iCnt = 1; iCnt <=  nConstraints; iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutRowIndex(12 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
		WM.PutColIndex(12 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}
	
	/* Recupera i dati che servono */
	Vec3 fn(pNode->GetRCurr()*tilde_fn);
	
	Vec3 Omega(pNode->GetWCurr());
	
	Vec3 bPrime(pNode->GetVCurr() + Omega.Cross(fn));
	
	/* F ed M sono gia' state aggiornate da InitialAssRes;
	 * Recupero FPrime e MPrime*/
	Vec3 MPrime(0.);
	Vec3 FPrime(0.);

	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		FPrime(iPosIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + iCnt);
	}
	
	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		MPrime(iRotIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + nPosConstraints + iCnt);
	}

	Vec3 FTmp(Rch*F);
	Vec3 FPrimeTmp(Rch*FPrime);

	WM.Add(3 + 1, 3 + 1, Mat3x3(FTmp, fn));
	WM.Add(9 + 1, 3 + 1, Mat3x3(FPrimeTmp, fn));

	/* Constraints: Add only active rows/columns*/	
	
	/* Positions contribution:*/

	Mat3x3 fnCross_Rch(Mat3x3(fn)*Rch); // = [ fn x ] * Rc

	for (unsigned iCnt = 0 ; iCnt < nPosConstraints; iCnt++) {
		Vec3 vRch(Rch.GetVec(iPosIncid[iCnt]));
		Vec3 vfnCross_Rch(fnCross_Rch.GetVec(iPosIncid[iCnt]));

		/* Equilibrium, node */
      		WM.Add(1, 12 + 1 + iCnt, vRch);			// * Delta_F
      		WM.Add(3 + 1, 12 + 1 + iCnt, vfnCross_Rch);	// * Delta_F

		/* Constraint, node */
      		WM.AddT(12 + 1 + iCnt, 1, vRch);		// * Delta_x2
      		WM.AddT(12 + 1 + iCnt, 3 + 1, vfnCross_Rch);	// * Delta_g2
	
		/* d/dt(Equilibrium), node */
      		WM.Add(6 + 1, 12 + 1 + nConstraints + iCnt, vRch);		// * Delta_FP
      		WM.Add(9 + 1, 12 + 1 + nConstraints + iCnt, vfnCross_Rch);	// * Delta_FP

		/* d/dt(Constraint), node */
      		WM.AddT(12 + 1 + nConstraints +  iCnt, 6 + 1, vRch);		// * Delta_xP2
      		WM.AddT(12 + 1 + nConstraints +  iCnt, 9 + 1, vfnCross_Rch);	// * Delta_W2
	}

	for (unsigned iCnt = 0 ; iCnt < nRotConstraints; iCnt++) {
		Vec3 vRchr(Rchr.GetVec(iRotIncid[iCnt]));

		/* Equilibrium, node */
      		WM.Add(3 + 1, 12 + 1 + nPosConstraints + iCnt, vRchr);	// * Delta_M

		/* Constraint, node */
      		WM.AddT(12 + 1 + nPosConstraints +  iCnt, 3 + 1, vRchr);	// * Delta_g2

		/* d/dt(Equilibrium), node */
      		WM.Add(9 + 1, 12 + 1 + nConstraints + nPosConstraints + iCnt, vRchr);	// * Delta_MP

		/* d/dt(Constraint), node */
      		WM.AddT(12 + 1 + nConstraints + nPosConstraints + iCnt, 9 + 1, vRchr);	// * Delta_W2
	}

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
TotalPinJoint::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{

	DEBUGCOUT("Entering TotalPinJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iNodeFirstVelIndex = iNodeFirstPosIndex + 6;
	integer iFirstReactionIndex = iGetFirstIndex();
	integer iReactionPrimeIndex = iFirstReactionIndex + nConstraints;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNodeFirstVelIndex + iCnt);
	}

	for (unsigned int iCnt = 1; iCnt <= nConstraints; iCnt++)	{
		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WorkVec.PutRowIndex(12 + nConstraints + iCnt, iReactionPrimeIndex + iCnt);
	}

	/* Recupera i dati */
	/* Recupera i dati che servono */
	Vec3 fn(pNode->GetRCurr()*tilde_fn);
	Mat3x3 Rnhr = pNode->GetRCurr()*tilde_Rnhr;
	Vec3 Omega(pNode->GetWCurr());
	
	Vec3 FPrime(0.);
	Vec3 MPrime(0.);

	/* Aggiorna F ed M, che restano anche per InitialAssJac */
	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		F(iPosIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + iCnt);
		FPrime(iPosIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + iCnt);
	}

	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		M(iRotIncid[iCnt]) = XCurr(iFirstReactionIndex + 1 + nPosConstraints + iCnt);
		MPrime(iRotIncid[iCnt]) = XCurr(iReactionPrimeIndex + 1 + nPosConstraints + iCnt);
	}

	Vec3 FTmp(Rch*F);
	Vec3 MTmp(Rchr*M);
	Vec3 FPrimeTmp(Rch*FPrime);
	Vec3 MPrimeTmp(Rchr*MPrime);

	/* Equilibrium, node 2 */
	WorkVec.Sub(1, FTmp);
	WorkVec.Sub(3 + 1, fn.Cross(FTmp) + MTmp);

	/* d/dt( Equilibrium ) , node 2 */
	WorkVec.Sub(6 + 1, FPrimeTmp);
	WorkVec.Sub(9 + 1, fn.Cross(FPrimeTmp) + MPrimeTmp);

	/* Constraint Equations */
	Vec3 XDelta = RchT*(pNode->GetXCurr() + fn) - tilde_Xc - XDrv.Get();
	Vec3 XDeltaPrime = RchT*(pNode->GetVCurr() + Omega.Cross(fn));
	
	if (XDrv.bIsDifferentiable())	{
		XDeltaPrime -= XDrv.GetP();
	}
	
	Mat3x3 R0T = RotManip::Rot(-ThetaDrv.Get());	// -Theta0 to get R0 transposed
	Mat3x3 RDelta = RchrT*Rnhr*R0T;
	ThetaDelta = RotManip::VecRot(RDelta);
	Vec3 ThetaDeltaPrime = RchrT*Omega;

	if (ThetaDrv.bIsDifferentiable())	{
		ThetaDeltaPrime -= RDelta*ThetaDrv.GetP();
	}

	/* Position constraint:  */
	for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
		WorkVec.PutCoef(12 + 1 + iCnt, -(XDelta(iPosIncid[iCnt])));
		WorkVec.PutCoef(12 + 1 + nConstraints + iCnt, -(XDeltaPrime(iPosIncid[iCnt])));
	}

	/* Rotation constraints: */
	for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
		WorkVec.PutCoef(12 + 1 + nPosConstraints + iCnt, -(ThetaDelta(iRotIncid[iCnt])));
		WorkVec.PutCoef(12 + 1 + nPosConstraints + nConstraints +  iCnt, -(ThetaDeltaPrime(iRotIncid[iCnt])));
	}

	return WorkVec;
}

unsigned int
TotalPinJoint::iGetNumPrivData(void) const
{
	return 18;
}

unsigned int
TotalPinJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned int off = 0;

	switch (s[0]) {
	case 'p':
		/* relative position */
		break;

	case 'r':
		/* relative orientation */
		off += 3;
		break;

	case 'F':
		/* force */
		off += 6;
		break;

	case 'M':
		/* moment */
		off += 9;
		break;

	case 'd':
		/* imposed relative position */
		off += 12;
		break;

	case 't':
		/* imposed relative orientation */
		off += 15;
		break;

	default:
		return 0;
	}

	switch (s[1]) {
	case 'x':
		return off + 1;

	case 'y':
		return off + 2;

	case 'z':
		return off + 3;
	}

	return 0;
}

doublereal
TotalPinJoint::dGetPrivData(unsigned int i) const
{
	switch (i) {
	case 1:
	case 2:
	case 3: {
		Vec3 x(RchT*(pNode->GetXCurr() + pNode->GetRCurr()*tilde_fn) - tilde_Xc);
			return x(i);
		}

	case 4:
	case 5:
	case 6: {
		Vec3 Theta(Unwrap(ThetaDeltaPrev, ThetaDelta) + ThetaDrv.Get());
		return Theta(i - 3);
		}

	case 7:
	case 8:
	case 9:
		return F(i - 6);

	case 10:
	case 11:
	case 12:
		return M(i - 9);

	case 13:
	case 14:
	case 15:
		if (!bPosActive[i - 13]) {
			return 0.;
		}
		return XDrv.Get()(i - 12);

	case 16:
	case 17:
	case 18:
		if (!bRotActive[i - 16]) {
			return 0.;
		}
		return ThetaDrv.Get()(i - 15);

	default:
		ASSERT(0);
	}

	return 0.;
}

/* TotalPinJoint - end */

/* TotalForce - begin */

TotalForce::TotalForce(unsigned int uL, const DriveCaller *pDC,
	bool bForce[3],
	TplDriveCaller<Vec3> *const pDCForce,
	bool bCouple[3],
	TplDriveCaller<Vec3> *const pDCCouple,
	const StructNode *pN1,
	const Vec3& f1Tmp, const Mat3x3& R1hTmp, const Mat3x3& R1hrTmp,
	const StructNode *pN2,
	const Vec3& f2Tmp, const Mat3x3& R2hTmp, const Mat3x3& R2hrTmp,
	flag fOut)
: Elem(uL, fOut),
Force(uL, pDC, fOut),
pNode1(pN1), pNode2(pN2),
f1(f1Tmp), R1h(R1hTmp), R1hr(R1hrTmp),
f2(f2Tmp), R2h(R2hTmp), R2hr(R2hrTmp),
FDrv(pDCForce), MDrv(pDCCouple),
nForces(0), nCouples(0),
M(0.), F(0.)
{
	/* Active-Inactive components: */ 

	for (unsigned int i = 0; i < 3; i++) {
		bForceActive[i] = bForce[i];
		bCoupleActive[i] = bCouple[i];
	
		if (bForceActive[i]) {
			iForceIncid[nForces] = i + 1;
			nForces++;
		}
		if (bCoupleActive[i]) {
			iCoupleIncid[nCouples] = i + 1;
			nCouples++;
		}
	}
}

std::ostream&
TotalForce::Restart(std::ostream& out) const
{
	/*   Force::Restart(out) << ", total, "
	<< pNode1->GetLabel()
	<< ", reference, node, ",
	Dir.Write(out, ", ")
	<< ", reference, node, ",
	f1.Write(out, ", ") << ", "
	<< pNode2->GetLabel()
	<< ", reference, node, ",
	f2.Write(out, ", ") << ", ";
	*/
	return pGetDriveCaller()->Restart(out) << ';'      
	<< std::endl;
}

/* Assemblaggio jacobiano */
VariableSubMatrixHandler&
TotalForce::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering TotalForce::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	/* Moltiplica il momento e la forza per il coefficiente del metodo */
	Vec3 FTmp(R1*(F*dCoef));
	Vec3 MTmp(R1r*(M*dCoef));

	/* Equilibrium: ((Phi/q)^T*F)/q */

	Mat3x3 Tmp;

	/* [ F x ] */
	Tmp = Mat3x3(FTmp);

	/* Lines 1->3: */
	WM.Add(1, 3 + 1, Tmp);

	/* Lines 4->6: */
	WM.Sub(3 + 1, 1, Tmp);

	WM.Add(3 + 1, 6 + 1, Tmp);

	/* Lines 7->9: */
	WM.Sub(6 + 1, 3 + 1, Tmp);

	/* [ F x ] [ b2 x ] */
	Tmp = Mat3x3(FTmp, b2);

	/* Lines 4->6: */
	WM.Sub(3 + 1, 9 + 1, Tmp);

	/* Lines 10->12: */
	WM.Add(9 + 1, 9 + 1, Tmp);

	/* [ b1 x ] [ F x ] + [ M x ] */

	/* Lines 4->6: */
	WM.Add(3 + 1, 3 + 1, Mat3x3(b1, FTmp) + Mat3x3(MTmp));

	/* [ b2 x ] [ F x ] + [ M x ] */

	/* Lines 10->12: */
	WM.Sub(9 + 1, 3 + 1, Mat3x3(b2, FTmp) + Mat3x3(MTmp));

	return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler&
TotalForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering TotalJoint::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex + iCnt);
	}

	/* Get Forces */
	
	Vec3 FTmp(FDrv.Get());
	Vec3 MTmp(MDrv.Get());

	for (unsigned iCnt = 0; iCnt < nForces; iCnt++) {
		F(iForceIncid[iCnt]) = FTmp(iCnt + 1);
	}

	for (unsigned iCnt = 0; iCnt < nCouples; iCnt++) {
		M(iCoupleIncid[iCnt]) = MTmp(iCnt + 1);
	}

	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());

	Mat3x3 R1 = pNode1->GetRCurr()*R1h;
	Mat3x3 R1r = pNode1->GetRCurr()*R1hr;
	Mat3x3 R2r = pNode2->GetRCurr()*R2hr;

	FTmp = R1*F;
	MTmp = R1r*M;

	/* Equilibrium, node 1 */
	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, MTmp + b1.Cross(FTmp));

	/* Equilibrium, node 2 */
	WorkVec.Sub(6 + 1, FTmp);
	WorkVec.Sub(9 + 1, MTmp + b2.Cross(FTmp));

	return WorkVec;
}

VariableSubMatrixHandler&
TotalForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{

	/* Per ora usa la matrice piena; eventualmente si puo'
	 * passare a quella sparsa quando si ottimizza */
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;


	/* Setta gli indici dei nodi */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	
	Vec3 FTmp(R1 * F);
	Vec3 MTmp(R1r * M);

	Mat3x3 Tmp;

	/* [ F x ] */
	Tmp = Mat3x3(FTmp);

	/* Force, Node 1, Lines 1->3: */
	WM.Add(1, 3 + 1, Tmp);	// * Delta_g1

	/* Moment, Node 1, Lines 4->6: */
	WM.Sub(3 + 1, 1, Tmp);	// * Delta_x1

	WM.Add(3 + 1, 12 + 1, Tmp);	// * Delta_x2

	/* Force, Node 2, Lines 13->15: */
	WM.Sub(12 + 1, 3 + 1, Tmp);	// * Delta_g1
	
	/* [ F x ] [ b2 x ] */
	Tmp = Mat3x3(FTmp, b2);

	/* Moment, Node1, Lines 4->6: */
	WM.Sub(3 + 1, 15 + 1, Tmp);	// * Delta_g2

	/* Moment, Node2, Lines 16->18: */
	WM.Add(15 + 1, 15 + 1, Tmp);	// * Delta_g2
	
	/* [ b1 x ] [ F x ] + [ M x ] */

	/* Moment, Node1, Lines 4->6: */
	WM.Add(3 + 1, 3 + 1, Mat3x3(b1, FTmp) + Mat3x3(MTmp));	// *Delta_g1

	/* [ b2 x ] [ F x ] + [ M x ] */

	/* Moment, Node2, Lines 16->18: */
	WM.Sub(15 + 1, 3 + 1, Mat3x3(b2, FTmp) + Mat3x3(MTmp));	// * Delta_g1
	
	/* d/dt(Moment), Node2, Lines 22->24: */
	WM.Sub(21 + 1, 9 + 1, Mat3x3(b2, FTmp) + Mat3x3(MTmp));	// * Delta_W1
	
	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
TotalForce::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{

	DEBUGCOUT("Entering TotalJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	/* Recupera i dati */
	/* Recupera i dati che servono */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R1r(pNode1->GetRRef()*R1hr);
	Mat3x3 R2r(pNode2->GetRCurr()*R2hr);
	
	Vec3 b2(pNode2->GetRCurr()*f2);
	Vec3 b1(pNode2->GetXCurr() + b2 - pNode1->GetXCurr());
	
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	
	Mat3x3 Omega2Cross(Omega2);
	Vec3 Omega2Crossb2(Omega2Cross*b2);
	Vec3 b1Prime(pNode2->GetVCurr() + Omega2.Cross(b2) - pNode1->GetVCurr());
	
	/* Aggiorna F ed M, che restano anche per InitialAssJac */

	Vec3 FTmp(FDrv.Get());
	Vec3 MTmp(MDrv.Get());

	for (unsigned iCnt = 0; iCnt < nForces; iCnt++) {
		F(iForceIncid[iCnt]) = FTmp(iCnt + 1);
	}

	for (unsigned iCnt = 0; iCnt < nCouples; iCnt++) {
		M(iCoupleIncid[iCnt]) = MTmp(iCnt + 1);
	}

	FTmp = R1 * F;
	MTmp = R1r * M;

	/* Equilibrium, node 1 */
	WorkVec.Add(1, FTmp);
	WorkVec.Add(3 + 1, b1.Cross(FTmp) + MTmp);

	/* Equilibrium, node 2 */
	WorkVec.Sub(12 + 1, FTmp);
	WorkVec.Sub(15 + 1, b2.Cross(FTmp) + MTmp);

return WorkVec;
}


