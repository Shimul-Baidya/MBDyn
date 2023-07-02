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

FiniteDifferenceJacobianBase::FiniteDifferenceJacobianBase(DataManager* const pDM, FiniteDifferenceJacobianParam&& oParam)
     :FiniteDifferenceJacobianParam(std::move(oParam)),
      pDM(pDM),
      dTimePrev(-std::numeric_limits<doublereal>::max()),
      dTimeMaxDiff(-std::numeric_limits<doublereal>::max()),
      dMaxDiffAll(-std::numeric_limits<doublereal>::max()),
      iRowMaxDiffAll(-1),
      iColMaxDiffAll(-1),
      iJacobians(0),
      iIterations(0)
 {
}

FiniteDifferenceJacobianBase::~FiniteDifferenceJacobianBase()
{
     if (iJacobians > 0 && (uOutputFlags & FDJAC_OUTPUT_STAT_END)) {
          silent_cerr("Finite difference Jacobian matrix check #" << iJacobians << " Time=" << dTimeMaxDiff << ":\n");
          silent_cerr("maximum difference at Jac(" << iRowMaxDiffAll << "," << iColMaxDiffAll << "): " << dMaxDiffAll << "\n");
     }
}

void FiniteDifferenceJacobianBase::Output(const MatrixHandler* const pJac, const doublereal dTimeCurr, const doublereal dMaxDiff, const integer iRowMaxDiff, const integer iColMaxDiff)
{
     if (dMaxDiff >= dMaxDiffAll) {
          dTimeMaxDiff = dTimeCurr;
          dMaxDiffAll = dMaxDiff;
          iRowMaxDiffAll = iRowMaxDiff;
          iColMaxDiffAll = iColMaxDiff;
     }

     if (uOutputFlags & (FDJAC_OUTPUT_STAT_PER_ITER | FDJAC_OUTPUT_MAT_PER_ITER)) {
          silent_cerr("Finite difference Jacobian check #" << iJacobians << " Time=" << dTimeCurr << " Iteration(" << iIterations << "):\n");
     }

     if (uOutputFlags & FDJAC_OUTPUT_STAT_PER_ITER) {
          silent_cerr("maximum difference at Jac(" << iRowMaxDiff << "," << iColMaxDiff << "): " << dMaxDiff << "\n");
     }

     if (pFDJac) {
          silent_cerr("\nxxxxxxxxxxxxxxx\n\n");
          silent_cerr(*pJac << "\n");
          silent_cerr("\n---------------\n\n");
          silent_cerr(*pFDJac << "\n");
          silent_cerr("\n===============\n\n");
     }
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
     if (inc.iGetSize() != pJac->iGetNumCols()) {
          if (uOutputFlags & FDJAC_OUTPUT_MAT_PER_ITER) {
               pFDJac.reset(new FullMatrixHandler(pJac->iGetNumRows(), pJac->iGetNumCols()));
          }

          inc.Resize(pJac->iGetNumCols());

          for (integer k = 0; k < N; ++k) {
               incsol[k].Resize(pDM->iGetNumDofs());
          }
     }

     ASSERT(!pFDJac || pFDJac->iGetNumRows() == pJac->iGetNumRows());
     ASSERT(!pFDJac || pFDJac->iGetNumCols() == pJac->iGetNumCols());
     ASSERT(inc.iGetSize() == pJac->iGetNumCols());
     ASSERT(incsol[0].iGetSize() == pJac->iGetNumRows());
}

template <integer N>
void FiniteDifferenceJacobian<N>::JacobianCheck(const NonlinearProblem* const pNLP, const MatrixHandler* const pJac)
{
     // Finite difference check of Jacobian matrix
     // NOTE: might not be safe!

     const doublereal dTimeCurr = pDM->dGetTime();

     if (dTimeCurr != dTimePrev) {
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

          doublereal dMaxDiff = -std::numeric_limits<doublereal>::max();
          integer iRowMaxDiff = -1;
          integer iColMaxDiff = -1;

          for (integer j = 1; j <= pJac->iGetNumCols(); ++j) {
               inc.PutCoef(j, pertFD[0] * dFDJacCoef);

               ASSERT(incsol.size() >= idxFD.size());
               static_assert(pertFD.size() >= idxFD.size());

               for (size_t k = 0; k < idxFD.size(); ++k) {
                    pNLP->Update(&inc);

                    ASSERT(incsol[idxFD[k]].iGetSize() == pJac->iGetNumRows());

                    incsol[idxFD[k]].Reset();
                    pNLP->Residual(&incsol[idxFD[k]]);

                    inc.PutCoef(j, (k + 1 < pertFD.size()) ? (pertFD[k + 1] - pertFD[k]) * dFDJacCoef : 0.);
               }

               doublereal normColj = 0.;

               for (integer i = 1; i <= pJac->iGetNumRows(); ++i) {
                    doublereal Jacij = 0.;

                    static_assert(idxFD.size() >= coefFD.size());

                    for (size_t k = 0; k < coefFD.size(); ++k) {
                         Jacij -= incsol[idxFD[k]](i) * coefFD[k];
                    }

                    Jacij /=  dFDJacCoef;

                    normColj += Jacij * Jacij;

                    incsol[0](i) = Jacij;

                    if (pFDJac) {
                         pFDJac->PutCoef(i, j, Jacij);
                    }
               }

               normColj = sqrt(normColj);

               for (integer i = 1; i <= pJac->iGetNumRows(); ++i) {
                    const doublereal dCurrDiff = fabs(incsol[0](i) - pJac->dGetCoef(i, j)) / normColj;

                    if (dCurrDiff > dMaxDiff) {
                         iRowMaxDiff = i;
                         iColMaxDiff = j;
                         dMaxDiff = dCurrDiff;
                    }
               }
          }

          Output(pJac, dTimeCurr, dMaxDiff, iRowMaxDiff, iColMaxDiff);
     }
     dTimePrev = dTimeCurr;
}

template
class FiniteDifferenceJacobian<2>;

template
class FiniteDifferenceJacobian<3>;

template
class FiniteDifferenceJacobian<5>;

template
class FiniteDifferenceJacobian<7>;
