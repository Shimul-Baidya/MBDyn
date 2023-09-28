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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2023(-2023) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <algorithm>
#include "dataman.h"
#include "fdjac.h"
#include "fullmh.h"

constexpr integer FiniteDifferenceOperator<2>::N;
constexpr std::array<doublereal, 2> FiniteDifferenceOperator<2>::pertFD;
constexpr std::array<integer, 2> FiniteDifferenceOperator<2>::idxFD;
constexpr std::array<doublereal, 2> FiniteDifferenceOperator<2>::coefFD;

constexpr integer FiniteDifferenceOperator<3>::N;
constexpr std::array<doublereal, 3> FiniteDifferenceOperator<3>::pertFD;
constexpr std::array<integer, 3> FiniteDifferenceOperator<3>::idxFD;
constexpr std::array<doublereal, 3 - 1> FiniteDifferenceOperator<3>::coefFD;

constexpr integer FiniteDifferenceOperator<5>::N;
constexpr std::array<doublereal, 5> FiniteDifferenceOperator<5>::pertFD;
constexpr std::array<integer, 5> FiniteDifferenceOperator<5>::idxFD;
constexpr std::array<doublereal, 5 - 1> FiniteDifferenceOperator<5>::coefFD;

constexpr integer FiniteDifferenceOperator<7>::N;
constexpr std::array<doublereal, 7> FiniteDifferenceOperator<7>::pertFD;
constexpr std::array<integer, 7> FiniteDifferenceOperator<7>::idxFD;
constexpr std::array<doublereal, 7 - 1> FiniteDifferenceOperator<7>::coefFD;

FiniteDifferenceJacobianBase::JacobianStat::JacobianStat()
     :dTimeSample(-std::numeric_limits<doublereal>::max()),
      dMaxDiff(-std::numeric_limits<doublereal>::max()),
      iRowMaxDiff(-1),
      iColMaxDiff(-1) {
}

FiniteDifferenceJacobianBase::FiniteDifferenceJacobianBase(DataManager* const pDM, FiniteDifferenceJacobianParam&& oParam)
     :FiniteDifferenceJacobianParam(std::move(oParam)),
      pDM(pDM),
      dTimePrev(-std::numeric_limits<doublereal>::max()),
      iJacobians(0),
      iIterations(0)
 {
}

FiniteDifferenceJacobianBase::~FiniteDifferenceJacobianBase()
{
     if (iJacobians > 0 && (uOutputFlags & FDJAC_OUTPUT_STAT_END)) {
          silent_cerr("Finite difference Jacobian matrix check #" << iJacobians << " Time=" << oJacStatMax.dTimeSample << ":\n");
          silent_cerr("maximum difference at Jac(" << oJacStatMax.iRowMaxDiff << "," << oJacStatMax.iColMaxDiff << "): " << oJacStatMax.dMaxDiff << "\n");
     }
}

void FiniteDifferenceJacobianBase::JacobianCheck(const NonlinearProblem* const pNLP, const MatrixHandler* const pJac)
{
     // Finite difference check of Jacobian matrix
     // NOTE: might not be safe!

     JacobianStat oJacStatCurr;

     oJacStatCurr.dTimeSample = pDM->dGetTime();

     if (oJacStatCurr.dTimeSample != dTimePrev) {
          iIterations = 0;
     }

     ++iIterations;

     if (pFDJacMeterStep->dGet() && pFDJacMeterIter->dGet(iIterations)) {
          ++iJacobians;

          Attach(pJac);

          if (pFDJac) {
               pFDJac->Reset();
          }

          inc.Reset();

          JacobianCheckImpl(pNLP, pJac, oJacStatCurr);
     }

     dTimePrev = oJacStatCurr.dTimeSample;
}

void FiniteDifferenceJacobianBase::Output(const MatrixHandler* const pJac, const JacobianStat& oJacStatCurr)
{
     if (oJacStatCurr.dMaxDiff >= oJacStatMax.dMaxDiff) {
          oJacStatMax = oJacStatCurr;
     }

     if (uOutputFlags & (FDJAC_OUTPUT_STAT_PER_ITER | FDJAC_OUTPUT_MAT_PER_ITER)) {
          silent_cerr("Finite difference Jacobian check #" << iJacobians << " Time=" << oJacStatCurr.dTimeSample << " Iteration(" << iIterations << "):\n");
     }

     if (uOutputFlags & FDJAC_OUTPUT_STAT_PER_ITER) {
          silent_cerr("maximum difference at Jac(" << oJacStatCurr.iRowMaxDiff << "," << oJacStatCurr.iColMaxDiff << "): " << oJacStatCurr.dMaxDiff << "\n");
     }

     if (pFDJac) {
          silent_cerr("\nxxxxxxxxxxxxxxx\n\n");

          if (silent_err) {
               pJac->Print(std::cerr, MatrixHandler::MAT_PRINT_TRIPLET);
          }

          silent_cerr("\n\n---------------\n\n");

          if (silent_err) {
               pFDJac->Print(std::cerr, MatrixHandler::MAT_PRINT_TRIPLET);
          }

          silent_cerr("\n\n===============\n\n");
     }
}

void FiniteDifferenceJacobianBase::Attach(const MatrixHandler* pJac)
{
     if (inc.iGetSize() != pJac->iGetNumCols()) {
          if (uOutputFlags & FDJAC_OUTPUT_MAT_PER_ITER) {
               pFDJac.reset(new FullMatrixHandler(pJac->iGetNumRows(), pJac->iGetNumCols()));
          }

          inc.Resize(pJac->iGetNumCols());
     }

     ASSERT(!pFDJac || pFDJac->iGetNumRows() == pJac->iGetNumRows());
     ASSERT(!pFDJac || pFDJac->iGetNumCols() == pJac->iGetNumCols());
     ASSERT(inc.iGetSize() == pJac->iGetNumCols());
}

template <integer N>
FiniteDifferenceJacobian<N>::FiniteDifferenceJacobian(DataManager* const pDM, FiniteDifferenceJacobianParam&& oParam)
     :FiniteDifferenceJacobianBase(pDM, std::move(oParam))
{
}

template <integer N>
FiniteDifferenceJacobian<N>::~FiniteDifferenceJacobian()
{
}

template <integer N>
void FiniteDifferenceJacobian<N>::Attach(const MatrixHandler* pJac)
{
     FiniteDifferenceJacobianBase::Attach(pJac);

     if (incsol[0].iGetSize() != pJac->iGetNumCols()) {
          for (size_t k = 0; k < incsol.size(); ++k) {
               incsol[k].Resize(pDM->iGetNumDofs());
          }
     }

     ASSERT(incsol[0].iGetSize() == pJac->iGetNumRows());
}

template <integer N>
void FiniteDifferenceJacobian<N>::JacobianCheckImpl(const NonlinearProblem* const pNLP, const MatrixHandler* const pJac, JacobianStat& oJacStatCurr)
{
     // Finite difference check of Jacobian matrix
     // NOTE: might not be safe!

     for (integer j = 1; j <= pJac->iGetNumCols(); ++j) {
          const doublereal Xj = pDM->GetDofType(j) == DofOrder::DIFFERENTIAL ? pDM->GetpXPCurr()->dGetCoef(j) : pDM->GetpXCurr()->dGetCoef(j);
          const doublereal h = dFDJacCoef * (1. + fabs(Xj));

          inc.PutCoef(j, pertFD[0] * h);

          ASSERT(incsol.size() >= idxFD.size());
          static_assert(pertFD.size() >= idxFD.size(), "size mismatch");

          for (size_t k = 0; k < idxFD.size(); ++k) {
               pNLP->Update(&inc);

               ASSERT(incsol[idxFD[k]].iGetSize() == pJac->iGetNumRows());

               incsol[idxFD[k]].Reset();
               pNLP->Residual(&incsol[idxFD[k]]);

               inc.PutCoef(j, (k + 1 < pertFD.size()) ? (pertFD[k + 1] - pertFD[k]) * h : 0.);
          }

          doublereal normColj = 0.;

          for (integer i = 1; i <= pJac->iGetNumRows(); ++i) {
               doublereal Jacij = 0.;

               static_assert(idxFD.size() >= coefFD.size(), "size mismatch");

               for (size_t k = 0; k < coefFD.size(); ++k) {
                    Jacij -= incsol[idxFD[k]](i) * coefFD[k];
               }

               Jacij /=  h;

               normColj += Jacij * Jacij;

               incsol[0](i) = Jacij;

               if (pFDJac) {
                    pFDJac->PutCoef(i, j, Jacij);
               }
          }

          normColj = sqrt(normColj);

          for (integer i = 1; i <= pJac->iGetNumRows(); ++i) {
               const doublereal dCurrDiff = fabs(incsol[0](i) - pJac->dGetCoef(i, j)) / normColj;

               if (dCurrDiff > oJacStatCurr.dMaxDiff) {
                    oJacStatCurr.iRowMaxDiff = i;
                    oJacStatCurr.iColMaxDiff = j;
                    oJacStatCurr.dMaxDiff = dCurrDiff;
               }
          }
     }

     Output(pJac, oJacStatCurr);
}

template
class FiniteDifferenceJacobian<2>;

template
class FiniteDifferenceJacobian<3>;

template
class FiniteDifferenceJacobian<5>;

template
class FiniteDifferenceJacobian<7>;

AdForwardModeJacobian::AdForwardModeJacobian(DataManager* pDM, FiniteDifferenceJacobianParam&& oParam)
     :FiniteDifferenceJacobianBase(pDM, std::move(oParam)) {
}

AdForwardModeJacobian::~AdForwardModeJacobian()
{

}

void AdForwardModeJacobian::JacobianCheckImpl(const NonlinearProblem* pNLP, const MatrixHandler* pJac, JacobianStat& oJacStatCurr)
{
     for (integer j = 1; j <= pJac->iGetNumCols(); ++j) {
          inc.PutCoef(j, 1.);
          pNLP->Jacobian(&JY, &inc);

          doublereal normColj = 0.;

          for (integer i = 1; i <= pJac->iGetNumRows(); ++i) {
               const doublereal Jacij = JY.dGetCoef(i);

               normColj += Jacij * Jacij;

               if (pFDJac) {
                    pFDJac->PutCoef(i, j, Jacij);
               }
          }

          normColj = sqrt(normColj);

          for (integer i = 1; i <= pJac->iGetNumRows(); ++i) {
               const doublereal dCurrDiff = fabs(JY.dGetCoef(i) - pJac->dGetCoef(i, j)) / normColj;

               if (dCurrDiff > oJacStatCurr.dMaxDiff) {
                    oJacStatCurr.iRowMaxDiff = i;
                    oJacStatCurr.iColMaxDiff = j;
                    oJacStatCurr.dMaxDiff = dCurrDiff;
               }
          }

          inc.PutCoef(j, 0.);
     }

     Output(pJac, oJacStatCurr);
}

void AdForwardModeJacobian::Attach(const MatrixHandler* pJac)
{
     FiniteDifferenceJacobianBase::Attach(pJac);

     if (JY.iGetSize() != pJac->iGetNumRows()) {
          JY.Resize(pJac->iGetNumRows());
     }
}
