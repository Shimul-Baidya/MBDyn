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

#ifndef FDJAC_H__INCLUDED
#define FDJAC_H__INCLUDED

#include <array>
#include <memory>

#include "vh.h"
#include "mh.h"
#include "datamanforward.h"
#include "drive.h"

template <integer N>
struct FiniteDifferenceOperator;

template <>
struct FiniteDifferenceOperator<2> {
     static constexpr integer N = 2;
     static constexpr std::array<doublereal, N> pertFD{1., 0.};
     static constexpr std::array<integer, N> idxFD{1, 0};
     static constexpr std::array<doublereal, N> coefFD{1., -1.};
};

template <>
struct FiniteDifferenceOperator<3> {
     static constexpr integer N = 3;
     static constexpr std::array<doublereal, N> pertFD{-1., 1., 0.};
     static constexpr std::array<integer, N> idxFD{1, 2, 0};
     static constexpr std::array<doublereal, N - 1> coefFD{-1./2., 1./2.};
};

template <>
struct FiniteDifferenceOperator<5> {
     static constexpr integer N = 5;
     static constexpr std::array<doublereal, N> pertFD{-2., -1., 1., 2., 0.};
     static constexpr std::array<integer, N> idxFD{1, 2, 3, 4, 0};
     static constexpr std::array<doublereal, N - 1> coefFD{1. / 12., -8. / 12., 8. / 12., -1. / 12.};
};

template <>
struct FiniteDifferenceOperator<7> {
     static constexpr integer N = 7;
     static constexpr std::array<doublereal, N> pertFD{-3., -2., -1., 1., 2., 3., 0.};
     static constexpr std::array<integer, N> idxFD{1, 2, 3, 4, 5, 6, 0};
     static constexpr std::array<doublereal, N - 1> coefFD{-1. / 60., 9. / 60., -45. / 60., 45. / 60., -9./60., 1./60.};
};

struct FiniteDifferenceJacobianParam {
     enum OutputFlags: unsigned {
          FDJAC_OUTPUT_NONE = 0x0u,
          FDJAC_OUTPUT_STAT_PER_ITER = 0x1u,
          FDJAC_OUTPUT_MAT_PER_ITER = 0x2u,
          FDJAC_OUTPUT_STAT_END = 0x4,
          FDJAC_OUTPUT_ALL = 0xFFFFFFFFu,
     };

     FiniteDifferenceJacobianParam()
          :dFDJacCoef(1e-4),
           uOutputFlags(FDJAC_OUTPUT_ALL) {
     }

     FiniteDifferenceJacobianParam(FiniteDifferenceJacobianParam&& oParam)
          :pFDJacMeterStep(std::move(oParam.pFDJacMeterStep)),
           pFDJacMeterIter(std::move(oParam.pFDJacMeterIter)),
           dFDJacCoef(oParam.dFDJacCoef),
           uOutputFlags(oParam.uOutputFlags) {
     }

     std::unique_ptr<DriveCaller> pFDJacMeterStep;
     std::unique_ptr<DriveCaller> pFDJacMeterIter;
     doublereal dFDJacCoef;
     unsigned uOutputFlags;
};

class FiniteDifferenceJacobianBase: public FiniteDifferenceJacobianParam {
public:
     FiniteDifferenceJacobianBase(DataManager* pDM, FiniteDifferenceJacobianParam&& oParam);
     virtual ~FiniteDifferenceJacobianBase();

     void JacobianCheck(const NonlinearProblem* pNLP, const MatrixHandler* pJac);

protected:
     struct JacobianStat {
          JacobianStat();
          doublereal dTimeSample;
          doublereal dMaxDiff;
          integer iRowMaxDiff;
          integer iColMaxDiff;
     };

     virtual void JacobianCheckImpl(const NonlinearProblem* pNLP, const MatrixHandler* pJac, JacobianStat& oJacStatCurr)=0;

     virtual void Attach(const MatrixHandler* pJac);

     void
     Output(const MatrixHandler* pJac,
            const JacobianStat& oJacStatCurr);

     DataManager* const pDM;
     std::unique_ptr<MatrixHandler> pFDJac;
     MyVectorHandler inc;

     doublereal dTimePrev;

     JacobianStat oJacStatMax;

     integer iJacobians;
     integer iIterations;
};

// Allow us to validate the Jacobian
template <integer N>
class FiniteDifferenceJacobian: public FiniteDifferenceJacobianBase, private FiniteDifferenceOperator<N> {
public:
     FiniteDifferenceJacobian(DataManager* pDM, FiniteDifferenceJacobianParam&& oParam);
     virtual ~FiniteDifferenceJacobian();

private:
     virtual void JacobianCheckImpl(const NonlinearProblem* pNLP, const MatrixHandler* pJac, JacobianStat& oJacStatCurr) override final;

     virtual void Attach(const MatrixHandler* pJac) override final;

     using FiniteDifferenceOperator<N>::pertFD;
     using FiniteDifferenceOperator<N>::idxFD;
     using FiniteDifferenceOperator<N>::coefFD;

     std::array<MyVectorHandler, N> incsol;
};

// Allow us to validate if sp_grad::SpGradient and sp_grad:GpGradProd are consistent
class AdForwardModeJacobian: public FiniteDifferenceJacobianBase {
public:
     AdForwardModeJacobian(DataManager* pDM, FiniteDifferenceJacobianParam&& oParam);
     virtual ~AdForwardModeJacobian();

private:
     virtual void JacobianCheckImpl(const NonlinearProblem* pNLP, const MatrixHandler* pJac, JacobianStat& oJacStatCurr) override final;

     virtual void Attach(const MatrixHandler* pJac) override final;

     MyVectorHandler JY;
};

#endif
