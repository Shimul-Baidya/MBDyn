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

/* inertia element */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <set>
#include <limits>

#include "inertia.h"
#include "dataman.h"
#include "Rot.hh"
#include <sstream> // not sure why this is now needed...

/* CenterOfMass - begin */

std::ostream&
CenterOfMass::Output_int(std::ostream& out) const
{
     return out
          << "    mass:        " << dMass << '\n'
          << "    J:           " << J << '\n'
          << "    Xcg:         " << X_cm << '\n'
          << "    Jcg:         " << J_cm << '\n'
          << "    Vcg:         " << V_cm << '\n'
          << "    Wcg:         " << Omega_cm << '\n';
}

void
CenterOfMass::Collect_int(void)
{
     dMass = 0.;
     S = Zero3;
     J = Zero3x3;

     B = Zero3;
     G_cm = Zero3;

     for (const ElemGravityOwner* pElem: elements)
     {
          dMass += pElem->dGetM();
          S += pElem->GetS();
          J += pElem->GetJ();

          B += pElem->GetB();
          G_cm += pElem->GetG();
     }

     J_cm = J;

     if (!dMass) {
          X_cm = Zero3;
          V_cm = Zero3;
          Omega_cm = Zero3;
          J_cm = J_cm.Symm();

     } else {
          X_cm = S/dMass;
          V_cm = B/dMass;
          G_cm -= X_cm.Cross(B);

          /*
           * FIXME: should also rotate it in the principal
           * reference frame, and log the angles
           */
          J_cm += Mat3x3(MatCrossCross, S, X_cm);
          J_cm = J_cm.Symm();

          ASSERT(J_cm.IsSymmetric()); // NOTE: should be a run time test
          Omega_cm = J_cm.LDLSolve(G_cm);
     }
}

/* Costruttore definitivo (da mettere a punto) */
CenterOfMass::CenterOfMass(std::set<const ElemGravityOwner *>&& elements) :
     elements(std::move(elements)),
     dMass(0.), S(Zero3), J(Zero3x3), X_cm(Zero3), Omega_cm(Zero3),
     J_cm(Zero3x3), B(Zero3), G_cm(Zero3)
{
     NO_OP;
}

CenterOfMass::~CenterOfMass(void)
{
     NO_OP;
}

/* CenterOfMass - end */

/* Inertia - begin */

/* momento statico */
Vec3
Inertia::GetS_int(void) const
{
     return S;
}

/* momento d'inerzia */
Mat3x3
Inertia::GetJ_int(void) const
{
     return J;
}

std::ostream&
Inertia::Output_int(std::ostream& out) const
{
     out
          << "inertia: " << GetLabel()
          << ( GetName().empty() ? "" : ( std::string(" \"") + GetName() + "\"" ) )
          << '\n';
     Vec3 DX(X_cm - X0);
     Mat3x3 JX(J_cm - Mat3x3(MatCrossCross, DX, DX*dMass));
     CenterOfMass::Output_int(out)
          << "    Xcg-X:       " << DX << '\n'
          << "    R^T*(Xcg-X): " << R0.MulTV(DX) << '\n'
          << "    J(X):        " << JX << '\n'
          << "    R^T*J(X)*R:  " << R0.MulTM(JX)*R0 << '\n'
          << "    Rp:          " << R_princ << '\n'
          << "    Thetap:      " << RotManip::VecRot(R_princ) << '\n'
          << "    Jp:          " << J_princ << '\n'
          << "    beta:        " << B << '\n'
          << "    gamma:       " << G_cm << '\n';

     //printf("\n\nmassa %1.16e\n",dMass);
     //Vec3 tmp = R0.MulTV(X_cm-X0);
     //printf("xcg %1.16e %1.16e %1.16e\n", tmp(1), tmp(2), tmp(3) );
     //printf("RR %1.16e %1.16e %1.16e\n", R_princ(1,1), R_princ(1,2), R_princ(1,3) );
     //printf("RR %1.16e %1.16e %1.16e\n", R_princ(2,1), R_princ(2,2), R_princ(2,3) );
     //printf("RR %1.16e %1.16e %1.16e\n", R_princ(3,1), R_princ(3,2), R_princ(3,3) );
     //printf("JJ %1.16e %1.16e %1.16e\n\n\n", J_princ(1), J_princ(2), J_princ(3) );
     return out;
}

/* Costruttore definitivo (da mettere a punto) */
Inertia::Inertia(unsigned int uL, const std::string& sN, std::set<const ElemGravityOwner *>&& elements,
                 const Vec3& x0, const Mat3x3& r0, flag fOut)
     : Elem(uL, fOut),
       ElemGravityOwner(uL, fOut),
       InitialAssemblyElem(uL, fOut),
       CenterOfMass(std::move(elements)),
       X0(x0), R0(r0), J0(Zero3x3), R_princ(Zero3x3), J_princ(Zero3)
{
     this->PutName(sN);

     Collect(CollectType::FIRST_TIME);

     R_princ_Prev = R_princ;
}

Inertia::~Inertia(void)
{
     NO_OP;
}

/* massa totale */
doublereal
Inertia::dGetM(void) const
{
     return dMass;
}

/* Tipo dell'elemento (usato solo per debug ecc.) */
Elem::Type
Inertia::GetElemType(void) const
{
     return Elem::INERTIA;
}

/* Numero gdl durante l'assemblaggio iniziale */
unsigned int
Inertia::iGetInitialNumDof(void) const
{
     return 0;
}

void
Inertia::Collect(const CollectType eCollectType)
{
     CenterOfMass::Collect_int();

     if (!dMass) {
          silent_cerr("Inertia(" << GetLabel() << "): "
                      "mass is null\n");

          R_princ = Eye3;
          J_princ = Zero3;
     } else {
          bool bStatus;

          if (eCollectType == CollectType::FIRST_TIME) {
               DEBUGCERR("Initial update of principal axes:\n");

               bStatus = J_cm.PrincipalAxes(J_princ, R_princ);
          } else {
               DEBUGCERR("Incremental update of principal axes:\n");

               ASSERT(eCollectType == CollectType::INCREMENTAL);

               bStatus = RotationIncrement();

               if (!bStatus) {
                    DEBUGCERR("Incremental update of principal axes failed!\n"
                              "Perform an initial update of principal axes instead:\n");
                    // FIXME: This could again cause discontinuities of the principal axes
                    bStatus = J_cm.PrincipalAxes(J_princ, R_princ);
               }
          }

          if (!bStatus) {
               pedantic_cerr("Inertia(" << GetLabel() << "): warning: principal axes of inertia tensor could not be computed!\n");
          }
     }

     const Vec3 DX = X_cm - X0;

     J0 = R0.MulTM((J_cm - Mat3x3(MatCrossCross, DX, DX * dMass)) * R0);
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
Inertia::Restart(std::ostream& out) const
{
     return out;
}

void
Inertia::Output(OutputHandler& OH) const
{
     if (bToBeOutput()) {
          const Vec3 DX(X_cm - X0);
          const Vec3 dx = R0.MulTV(DX);
          const Vec3 Phip = RotManip::VecRot(R_princ);
#ifdef USE_NETCDF
          if (OH.UseNetCDF(OutputHandler::INERTIA_ELEMENTS)) {
               OH.WriteNcVar(Var_B, B);
               OH.WriteNcVar(Var_G_cm, G_cm);
               OH.WriteNcVar(Var_dMass, dMass);
               OH.WriteNcVar(Var_X_cm, X_cm);
               OH.WriteNcVar(Var_V_cm, V_cm);
               OH.WriteNcVar(Var_Omega_cm, Omega_cm);

               OH.WriteNcVar(Var_DX, DX);
               OH.WriteNcVar(Var_dx, dx);
               OH.WriteNcVar(Var_Jp, J_princ);
               OH.WriteNcVar(Var_Phip, Phip);
               OH.WriteNcVar(Var_J_cm, J_cm);
               OH.WriteNcVar(Var_J0, J0);
          }
#endif // USE_NETCDF
          if (OH.UseText(OutputHandler::INERTIA_ELEMENTS)) {
               OH.InertiaElements() << std::setw(8) << GetLabel() // [1]
                                    << " " << B                   // [2-4]
                                    << " " << G_cm                // [5-7] With respect to center of mass instead of position of node
                                    << " " << dMass               // [8]
                                    << " " << X_cm                // [9-11]
                                    << " " << V_cm                // [12-14]
                                    << " " << Omega_cm            // [15-17]
                                    << " " << DX                  // [18-20]
                                    << " " << dx                  // [21-23]
                                    << " " << J_princ             // [24-26]
                                    << " " << Phip                // [27-29]
                                    << " " << J_cm                // [30-38]
                                    << " " << J0                  // [39-47]
                                    << "\n";
          }
     }
}

void
Inertia::OutputPrepare_int(OutputHandler &OH)
{
     if (bToBeOutput()) {
#ifdef USE_NETCDF
          ASSERT(OH.IsOpen(OutputHandler::NETCDF));

          std::ostringstream os;
          os << "elem.inertia." << GetLabel();
          m_sOutputNameBase = os.str();
          (void)OH.CreateVar(m_sOutputNameBase, "inertia");
#endif // USE_NETCDF
     }
}

void
Inertia::OutputPrepare(OutputHandler &OH)
{
     if (bToBeOutput()) {
#ifdef USE_NETCDF
          if (OH.UseNetCDF(OutputHandler::INERTIA_ELEMENTS))
          {
               OutputPrepare_int(OH);

               Var_B = OH.CreateVar<Vec3>(m_sOutputNameBase + ".B",
                                          OutputHandler::Dimensions::Momentum,
                                          "total momentum (x, y, z) w.r.t. global frame");

               Var_G_cm = OH.CreateVar<Vec3>(m_sOutputNameBase + ".G_cm",
                                             OutputHandler::Dimensions::MomentaMoment,
                                             "total momenta moment (x, y, z) at center of mass w.r.t. global frame");

               Var_dMass = OH.CreateVar<doublereal>(m_sOutputNameBase + "." "M",
                                                    OutputHandler::Dimensions::Mass,
                                                    "total mass");
               Var_X_cm = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "X_cm",
                                             OutputHandler::Dimensions::Length,
                                             "position of center of mass (x, y, z) w.r.t. global frame");
               Var_V_cm = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "V_cm",
                                             OutputHandler::Dimensions::Velocity,
                                             "velocity of center of mass (x, y, z) w.r.t. global frame");
               Var_Omega_cm = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "Omega_cm",
                                                 OutputHandler::Dimensions::AngularVelocity,
                                                 "angular velocity around center of mass (x, y, z) w.r.t. global frame");

               Var_DX = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "DX",
                                           OutputHandler::Dimensions::Length,
                                           "relative center of mass position, global frame (x, y, z)");
               Var_dx = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "dx",
                                           OutputHandler::Dimensions::Length,
                                           "relative center of mass position, local frame (x, y, z)");
               Var_Jp = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "Jp",
                                           OutputHandler::Dimensions::MomentOfInertia,
                                           "global inertia matrix, w.r.t. principal axes");
               Var_Phip = OH.CreateVar<Vec3>(m_sOutputNameBase + "." "Phip",
                                             OutputHandler::Dimensions::Dimensionless,
                                             "orientation vector of principal axes, global frame");

               Var_J_cm = OH.CreateVar<Mat3x3>(m_sOutputNameBase + ".J_cm",
                                               OutputHandler::Dimensions::MomentOfInertia,
                                               "matrix of inertia at center of mass w.r.t. global frame");

               Var_J0 = OH.CreateVar<Mat3x3>(m_sOutputNameBase + ".J0",
                                             OutputHandler::Dimensions::MomentOfInertia,
                                             "total inertia matrix with respect to reference position and reference orientation");
          }
#endif // USE_NETCDF
     }
}

void
Inertia::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 0;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
Inertia::AssJac(VariableSubMatrixHandler& WorkMat,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr)
{
     WorkMat.SetNullMatrix();
     return WorkMat;
}

SubVectorHandler&
Inertia::AssRes(SubVectorHandler& WorkVec,
                doublereal dCoef,
                const VectorHandler& XCurr,
                const VectorHandler& XPrimeCurr)
{
     WorkVec.Resize(0);

     Collect(CollectType::INCREMENTAL);

     return WorkVec;
}

bool
Inertia::bInverseDynamics(void) const
{
     return true;
}

/* Inverse Dynamics residual assembly */
SubVectorHandler&
Inertia::AssRes(SubVectorHandler& WorkVec,
                const VectorHandler& XCurr,
                const VectorHandler&  XPrimeCurr,
                const VectorHandler&  /* XPrimePrimeCurr */,
                InverseDynamics::Order iOrder)
{
     WorkVec.Resize(0);

     Collect(CollectType::INCREMENTAL);

     return WorkVec;
}

/* inverse dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
Inertia::AssJac(VariableSubMatrixHandler& WorkMat,
                const VectorHandler& /* XCurr */)
{
     WorkMat.SetNullMatrix();
     return WorkMat;
}

/* Dimensione del workspace durante l'assemblaggio iniziale.
 * Occorre tener conto del numero di dof che l'elemento definisce
 * in questa fase e dei dof dei nodi che vengono utilizzati.
 * Sono considerati dof indipendenti la posizione e la velocita'
 * dei nodi */
void
Inertia::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 0;
     *piNumCols = 0;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
Inertia::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                       const VectorHandler& XCurr)
{
     WorkMat.SetNullMatrix();
     return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
Inertia::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
     WorkVec.Resize(0);

     Collect(CollectType::INCREMENTAL);

     return WorkVec;
}

/* Usata per inizializzare la quantita' di moto */
void
Inertia::SetValue(DataManager *pDM,
                  VectorHandler& X, VectorHandler& XP,
                  SimulationEntity::Hints *ph)
{
     R_princ_Prev = R_princ;
}

unsigned int
Inertia::iGetNumPrivData(void) const
{
     return
          +3		// X[1-3]
          +3		// Phi[1-3]
          +3		// V[1-3]
          +3		// Omega[1-3]
          +3		// JP[1-3] principal inertia moments
          +3*3		// J[1-3,1-3] inertia tensor
          +1		// mass
          +3            // beta[1-3]
          +3            // gamma[1-3]
          +3*3          // Jcg[1-3,1-3] inertia tensor with respect to center of mass
          ;
}

unsigned int
Inertia::iGetPrivDataIdx(const char *s) const
{
        unsigned int idx = 0;
        bool bIsMatrix = false;

        switch (s[0]) {
        case 'X':
                break;

        case 'P':
                if (strncasecmp(s, "Phi", STRLENOF("Phi")) != 0) {
                        return 0;
                }
                s += STRLENOF("Phi") - 1;
                idx = 3;
                break;

        case 'V':
                idx = 6;
                break;

        case 'O':
                if (strncasecmp(s, "Omega", STRLENOF("Omega")) != 0) {
                        return 0;
                }
                s += STRLENOF("Omega") - 1;
                idx = 9;
                break;

        case 'J':
                if (strncasecmp(s, "JP", STRLENOF("JP")) == 0) {
                        s += STRLENOF("JP") - 1;
                        idx = 12;
                        break;
                }

                if (strncasecmp(s, "Jcg", STRLENOF("Jcg")) == 0) {
                        s += STRLENOF("Jcg") - 1;
                        idx = 31;
                        bIsMatrix = true;
                        break;
                }

                idx = 15;
                bIsMatrix = true;
                break;

        case 'm':
                if (s[1] != '\0') {
                        return 0;
                }
                return 25;

        case 'b':
                if (strncasecmp(s, "beta", STRLENOF("beta")) != 0) {
                        return 0;
                }
                s += STRLENOF("beta") - 1;
                idx = 25;
                break;

        case 'g':
                if (strncasecmp(s, "gamma", STRLENOF("gamma")) != 0) {
                        return 0;
                }
                s += STRLENOF("gamma") - 1;
                idx = 28;
                break;

        default:
                return 0;
        }

        s++;
        if (s[0] != '[') {
                return 0;
        }

        s++;
        switch (s[0]) {
        case '1':
        case '2':
        case '3':
                idx += s[0] - '0';
                s++;
                if (bIsMatrix) {
                        if (s[0] != ',') {
                                return 0;
                        }
                        s++;
                        switch (s[0]) {
                        case '1':
                        case '2':
                        case '3':
                                break;

                        default:
                                return 0;
                        }
                        idx += 3*(s[0] - '1');
                        s++;
                }
                break;

        default:
                return 0;
        }

        if (strcmp(s, "]") != 0) {
                return 0;
        }

        return idx;
}

doublereal
Inertia::dGetPrivData(unsigned int i) const
{
        switch (i) {
        case 1:
        case 2:
        case 3:
                // center of mass position
                return X_cm(i);
        case 4:
        case 5:
        case 6:
                // principal inertia axes Euler Rodriguez's parameters
                return RotManip::VecRot(R_princ)(i - 3);
        case 7:
        case 8:
        case 9:
                // center of mass velocity
                return V_cm(i - 6);
        case 10:
        case 11:
        case 12:
                // angular velocity
                return Omega_cm(i - 9);
        case 13:
        case 14:
        case 15:
                // principal inertia moments
                return J_princ(i - 12);
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
        case 21:
        case 22:
        case 23:
        case 24: {
                // inertia tensor with respect to X0 and R0
                int ir = (i - 15 - 1)%3 + 1;
                int ic = (i - 15 - 1)/3 + 1;

                return J0(ir, ic);
        }
        case 25:
                // mass
                return dMass;
        case 26:
        case 27:
        case 28:
                // momentum
                return B(i - 25);
        case 29:
        case 30:
        case 31:
                // momenta moment with respect to center of mass
                return G_cm(i - 28);
        case 32:
        case 33:
        case 34:
        case 35:
        case 36:
        case 37:
        case 38:
        case 39:
        case 40: {
                // global inertia tensor with respect to center of mass
                int ir = (i - 31 - 1)%3 + 1;
                int ic = (i - 31 - 1)/3 + 1;

                return J_cm(ir, ic);
        }
        default:
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

Vec3 Inertia::RotationParamRes(const Vec3& g) const
{
     const Mat3x3 RDelta(CGR_Rot::MatR, g);
     const Mat3x3 RCurr = RDelta * R_princ_Prev;
     const Mat3x3 JBar = RCurr.MulTM(J_cm * RCurr);

     return Vec3(JBar(1, 2) + JBar(2, 1), JBar(1, 3) + JBar(3, 1), JBar(2, 3) + JBar(3, 2)) / (2. * J_cm.Trace());
}

Mat3x3 Inertia::RotationParamJac(const Vec3& g, const Vec3& F0) const
{
     Mat3x3 Jac;

     const doublereal dg = (g.Norm() + 1.) * sqrt(std::numeric_limits<doublereal>::epsilon());

     for (integer j = 1; j <= 3; ++j) {
          Vec3 gDelta = g;

          gDelta(j) += dg;

          const Vec3 F = RotationParamRes(gDelta);

          for (integer i = 1; i <= 3; ++i) {
               Jac(i, j) = (F(i) - F0(i)) / dg;
          }
     }

     return Jac;
}

bool Inertia::RotationIncrement()
{
     DEBUGCOUTFNAME("Inertia::RotationIncrement");

     const doublereal dSingularTol = fabs(std::numeric_limits<doublereal>::epsilon() * J_cm.Trace());
     const doublereal fTol = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.8);
     const doublereal SolTol = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.8);
     constexpr integer iMaxIter = 50;

     bool bConverged = false;

     integer iIter;
     Vec3 gCurr(Zero3);

     for (iIter = 0; iIter < iMaxIter; ++iIter) {
          const Vec3 F = RotationParamRes(gCurr);
          const Mat3x3 Jac = RotationParamJac(gCurr, F);
          const doublereal detJac = Jac.dDet();

          if (fabs(detJac) < dSingularTol) {
               return false;
          }

          const Vec3 deltag = Jac.Solve(detJac, F);

          const doublereal fRes = F.Norm();
          const doublereal fSol = deltag.Norm();

          gCurr -= deltag;

          DEBUGCERR("Iteration(" << iIter << "): fRes=" << fRes << "\n\tfSol=" << fSol << "\n\tcond=" << Jac.dDet() * Jac.Inv().dDet() << "\n");

          if (fRes < fTol && fSol < SolTol) {
               bConverged = true;
               break;
          }
     }

     R_princ = Mat3x3(CGR_Rot::MatR, gCurr) * R_princ_Prev;

     const Mat3x3 Jtmp = R_princ.MulTM(J_cm * R_princ);

     ASSERT(!bConverged || Jtmp.IsDiag(fTol * J_cm.Trace()));

     for (integer i = 1; i <= 3; ++i) {
          J_princ(i) = Jtmp(i, i);
     }

     DEBUGCERR("Incremental update of principal axes " << (bConverged ? "converged" : "did not converge") <<  " after " << iIter << " iterations\n");

     return bConverged;
}

void
Inertia::AfterConvergence(const VectorHandler& X,
                          const VectorHandler& XP)
{
     R_princ_Prev = R_princ;
}

void
Inertia::AfterConvergence(const VectorHandler& X,
                          const VectorHandler& XP,
                          const VectorHandler& XPP)
{
     Inertia::AfterConvergence(X, XP);
}
