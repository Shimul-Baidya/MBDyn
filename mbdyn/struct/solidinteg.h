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

#ifndef __SOLID_INTEG_H__INCLUDED___
#define __SOLID_INTEG_H__INCLUDED___

#include <cmath>

#include "sp_matrix_base.h"
#include "solid.h"

class Gauss2 {
public:
     static constexpr sp_grad::index_type iGaussOrder = 2;
     static constexpr doublereal ri[] = {0.577350269189626, -0.577350269189626};
     static constexpr doublereal alphai[] = {1.0, 1.0};
     static constexpr doublereal ri_lumped[] = {1., -1.};
     static constexpr doublereal alphai_lumped[] = {1.0, 1.0};
};

class Gauss3 {
public:
     static constexpr sp_grad::index_type iGaussOrder = 3;

     static constexpr doublereal ri[] = {0.774596669241483, 0., -0.774596669241483};
     static constexpr doublereal alphai[] = {0.555555555555556, 0.888888888888889, 0.555555555555556};
     static constexpr doublereal ri_lumped[] = {1., 0., -1.};
     static constexpr doublereal alphai_lumped[] = {2./3., 2./3., 2./3.};
};

class Gauss2x2: private Gauss2 {
public:
     static constexpr sp_grad::index_type iNumEvalPoints = iGaussOrder * iGaussOrder;

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);

private:
     static constexpr sp_grad::index_type ridx[] = {0, 0, 1, 1};
     static constexpr sp_grad::index_type sidx[] = {0, 1, 0, 1};
};

class Gauss3x3: private Gauss3 {
public:
     static constexpr sp_grad::index_type iNumEvalPoints = iGaussOrder * iGaussOrder;

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);

private:
     static constexpr sp_grad::index_type ridx[] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
     static constexpr sp_grad::index_type sidx[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
};

class CollocTria6h {
public:
     static constexpr sp_grad::index_type iNumEvalPoints = 7;

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);

private:
     static constexpr double A = 0.470142064105115;
     static constexpr double B = 0.101286507323456;
     static constexpr double P1 = 0.066197076394253;
     static constexpr double P2 = 0.062969590272413;
     static constexpr doublereal zeta[7] =  {1./3., A, 1. - 2. * A, A, B, 1. - 2. * B, B};
     static constexpr doublereal eta[7] = {1./3., A, A, 1. - 2. * A, B, B, 1. - 2. * B};
     static constexpr doublereal w[7] = {9./80., P1, P1, P1, P2, P2, P2};
};

class Gauss2x2x2: private Gauss2 {
public:
     static constexpr sp_grad::index_type iNumEvalPointsStiffness = std::pow(iGaussOrder, 3);
     static constexpr sp_grad::index_type iNumEvalPointsMass = iNumEvalPointsStiffness;
     static constexpr sp_grad::index_type iNumEvalPointsMassLumped = iNumEvalPointsStiffness;

     static inline void
     GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightStiffness(sp_grad::index_type i);

     static inline void
     GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r) {
          GetPositionStiffness(i, r);
     }

     static inline doublereal
     dGetWeightMass(sp_grad::index_type i) {
          return dGetWeightStiffness(i);
     }

     static inline void
     GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightMassLumped(sp_grad::index_type i);

private:
     static constexpr sp_grad::index_type ridx[] = {0, 0, 0, 0, 1, 1, 1, 1};
     static constexpr sp_grad::index_type sidx[] = {0, 0, 1, 1, 0, 0, 1, 1};
     static constexpr sp_grad::index_type tidx[] = {0, 1, 0, 1, 0, 1, 0, 1};
};

class Gauss3x3x3: private Gauss3 {
public:
     static constexpr sp_grad::index_type iNumEvalPointsStiffness = std::pow(iGaussOrder, 3);
     static constexpr sp_grad::index_type iNumEvalPointsMass = iNumEvalPointsStiffness;
     static constexpr sp_grad::index_type iNumEvalPointsMassLumped = iNumEvalPointsStiffness;

     static inline void
     GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightStiffness(sp_grad::index_type i);

     static inline void
     GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r) {
          GetPositionStiffness(i, r);
     }

     static inline doublereal
     dGetWeightMass(sp_grad::index_type i) {
          return dGetWeightStiffness(i);
     }

     static inline void
     GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightMassLumped(sp_grad::index_type i);

private:
     static constexpr sp_grad::index_type ridx[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2};
     static constexpr sp_grad::index_type sidx[] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2};
     static constexpr sp_grad::index_type tidx[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
};

class GaussH20r: public Gauss2x2x2, public Gauss3x3x3 {
public:
     using Gauss2x2x2::iNumEvalPointsStiffness;
     using Gauss3x3x3::iNumEvalPointsMass;
     using Gauss3x3x3::iNumEvalPointsMassLumped;

     using Gauss2x2x2::GetPositionStiffness;
     using Gauss2x2x2::dGetWeightStiffness;

     using Gauss3x3x3::GetPositionMass;
     using Gauss3x3x3::dGetWeightMass;

     using Gauss3x3x3::GetPositionMassLumped;
     using Gauss3x3x3::dGetWeightMassLumped;
};

class CollocPenta15 {
     static constexpr sp_grad::index_type M = 7;
     static constexpr sp_grad::index_type N = 3;
     static constexpr sp_grad::index_type M_lumped = 6;
public:
     static constexpr sp_grad::index_type iNumEvalPointsStiffness = M * N;
     static constexpr sp_grad::index_type iNumEvalPointsMass = iNumEvalPointsStiffness;
     static constexpr sp_grad::index_type iNumEvalPointsMassLumped = M_lumped * N;

     static inline void
     GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightStiffness(sp_grad::index_type i);

     static inline void
     GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r) {
          GetPositionStiffness(i, r);
     }

     static inline doublereal
     dGetWeightMass(sp_grad::index_type i) {
          return dGetWeightStiffness(i);
     }

     static inline void
     GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightMassLumped(sp_grad::index_type i);

private:
     static constexpr doublereal r1 = 0.1012865073235;
     static constexpr doublereal r2 = 0.7974269853531;
     static constexpr doublereal r3 = r1;
     static constexpr doublereal r4 = 0.4701420641051;
     static constexpr doublereal r5 = r4;
     static constexpr doublereal r6 = 0.0597158717898;
     static constexpr doublereal r7 = 0.3333333333333;
     static constexpr doublereal s1 = r1;
     static constexpr doublereal s2 = r1;
     static constexpr doublereal s3 = r2;
     static constexpr doublereal s4 = r6;
     static constexpr doublereal s5 = r4;
     static constexpr doublereal s6 = r4;
     static constexpr doublereal s7 = r7;
     static constexpr doublereal w1 = 0.1259391805448;
     static constexpr doublereal w2 = w1;
     static constexpr doublereal w3 = w1;
     static constexpr doublereal w4 = 0.1323941527885;
     static constexpr doublereal w5 = w4;
     static constexpr doublereal w6 = w4;
     static constexpr doublereal w7 = 0.2250000000001;
     static constexpr doublereal ri[] = {r1, r2, r3, r4, r5, r6, r7};
     static constexpr doublereal si[] = {s1, s2, s3, s4, s5, s6, s7};
     static constexpr doublereal wi[] = {w1, w2, w3, w4, w5, w6, w7};
     static constexpr doublereal ti[] = {0.774596669241483, 0., -0.774596669241483};
     static constexpr doublereal alphai[] = {0.555555555555556, 0.888888888888889, 0.555555555555556};

     static constexpr doublereal ri_lumped[] = {0, 1, 0, 0.5,  0.5,    0};
     static constexpr doublereal si_lumped[] = {0, 0, 1,   0,  0.5,  0.5};
     static constexpr doublereal wi_lumped[] = {1./6., 1./6., 1./6., 1./6., 1./6., 1./6.};
     static constexpr doublereal ti_lumped[] = {1., 0., -1.};
     static constexpr doublereal alphai_lumped[] = {2./3., 2./3., 2./3.};
};

class CollocTet10h {
public:
     static constexpr sp_grad::index_type iNumEvalPointsStiffness = 5;
     static constexpr sp_grad::index_type iNumEvalPointsMass = 15;
     static constexpr sp_grad::index_type iNumEvalPointsMassLumped = -1;

     static inline void
     GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightStiffness(sp_grad::index_type i);

     static inline void
     GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightMass(sp_grad::index_type i);

     static inline void
     GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeightMassLumped(sp_grad::index_type i);

private:
     static constexpr doublereal a1 = 0.25;
     static constexpr doublereal b1 = 1. / 6.;
     static constexpr doublereal c1 = 0.5;
     static constexpr doublereal d1 = -2. / 15.;
     static constexpr doublereal e1 = 3. / 40.;
     static constexpr sp_grad::index_type N1 = 5;
     static constexpr doublereal r1[] = {a1, b1, b1, b1, c1};
     static constexpr doublereal s1[] = {a1, b1, b1, c1, b1};
     static constexpr doublereal t1[] = {a1, b1, c1, b1, b1};
     static constexpr doublereal w1[] = {d1, e1, e1, e1, e1};

     static constexpr doublereal a2 = 0.25;
     static constexpr doublereal b2_1 = (7. + sqrt(15.)) / 34.;
     static constexpr doublereal b2_2 = (7. - sqrt(15.)) / 34.;
     static constexpr doublereal c2_1 = (13. - 3. * sqrt(15.)) / 34.;
     static constexpr doublereal c2_2 = (13. + 3. * sqrt(15.)) / 34.;
     static constexpr doublereal d2 = (5. - sqrt(15.)) / 20.;
     static constexpr doublereal e2 = (5. + sqrt(15.)) / 20.;
     static constexpr doublereal f2 = 8. / 405.;
     static constexpr doublereal g2 = (2665. - 14. * sqrt(15.)) / 226800.;
     static constexpr doublereal h2 = (2665. + 14. * sqrt(15.)) / 226800.;
     static constexpr doublereal i2 = 5. / 567.;
     static constexpr sp_grad::index_type N2 = 15;

     static constexpr doublereal r2[] = {a2, b2_1, b2_1, c2_1, b2_1, b2_2, b2_2, c2_2, b2_2, d2, e2, d2, e2, d2, e2};
     static constexpr doublereal s2[] = {a2, b2_1, c2_1, b2_1, b2_1, b2_2, c2_2, b2_2, b2_2, e2, d2, d2, e2, e2, d2};
     static constexpr doublereal t2[] = {a2, b2_1, b2_1, b2_1, c2_1, b2_2, b2_2, b2_2, c2_2, d2, d2, e2, d2, e2, e2};
     static constexpr doublereal w2[] = {f2, g2, g2, g2, g2, h2, h2, h2, h2, i2, i2, i2, i2, i2, i2};

     static constexpr doublereal a3 = (5. - sqrt(5.)) / 20.;
     static constexpr doublereal b3 = (5. + 3. * sqrt(5)) / 20.;
     static constexpr doublereal c3 = 1. / 24.;
     static constexpr sp_grad::index_type N3 = 4;

     static constexpr doublereal r3[] = {a3, a3, a3, b3};
     static constexpr doublereal s3[] = {a3, a3, b3, a3};
     static constexpr doublereal t3[] = {a3, b3, a3, a3};
     static constexpr doublereal w3[] = {c3, c3, c3, c3};
};

void
Gauss2x2::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);

     static_assert(sizeof(ri) / sizeof(ri[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints);

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
}

doublereal
Gauss2x2::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);

     static_assert(sizeof(alphai) / sizeof(alphai[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints);

     return alphai[ridx[i]] * alphai[sidx[i]];
}

void
Gauss3x3::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);

     static_assert(sizeof(ri) / sizeof(ri[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints);

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
}

doublereal
Gauss3x3::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);

     static_assert(sizeof(alphai) / sizeof(alphai[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints);

     return alphai[ridx[i]] * alphai[sidx[i]];
}

void
CollocTria6h::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);

     static_assert(sizeof(zeta) / sizeof(zeta[0]) == iNumEvalPoints);
     static_assert(sizeof(eta) / sizeof(eta[0]) == iNumEvalPoints);

     r(1) = zeta[i];
     r(2) = eta[i];
}

doublereal
CollocTria6h::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);

     static_assert(sizeof(w) / sizeof(w[0]) == iNumEvalPoints);

     return w[i];
}

void
Gauss2x2x2::GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(ri) / sizeof(ri[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsStiffness);

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
     r(3) = ri[tidx[i]];
}

doublereal
Gauss2x2x2::dGetWeightStiffness(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(alphai) / sizeof(alphai[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsStiffness);

     return alphai[ridx[i]] * alphai[sidx[i]] * alphai[tidx[i]];
}

void
Gauss2x2x2::GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMassLumped);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(ri_lumped) / sizeof(ri_lumped[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsMassLumped);

     r(1) = ri_lumped[ridx[i]];
     r(2) = ri_lumped[sidx[i]];
     r(3) = ri_lumped[tidx[i]];
}

doublereal
Gauss2x2x2::dGetWeightMassLumped(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMassLumped);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(alphai_lumped) / sizeof(alphai_lumped[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsMassLumped);

     return alphai_lumped[ridx[i]] * alphai_lumped[sidx[i]] * alphai_lumped[tidx[i]];
}

void
Gauss3x3x3::GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(ri) / sizeof(ri[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsStiffness);

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
     r(3) = ri[tidx[i]];
}

doublereal
Gauss3x3x3::dGetWeightStiffness(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(alphai) / sizeof(alphai[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsStiffness);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsStiffness);

     return alphai[ridx[i]] * alphai[sidx[i]] * alphai[tidx[i]];
}

void
Gauss3x3x3::GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMassLumped);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(ri_lumped) / sizeof(ri_lumped[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsMassLumped);

     r(1) = ri_lumped[ridx[i]];
     r(2) = ri_lumped[sidx[i]];
     r(3) = ri_lumped[tidx[i]];
}

doublereal
Gauss3x3x3::dGetWeightMassLumped(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMassLumped);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(alphai_lumped) / sizeof(alphai_lumped[0]) == iGaussOrder);
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPointsMassLumped);
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPointsMassLumped);

     return alphai_lumped[ridx[i]] * alphai_lumped[sidx[i]] * alphai_lumped[tidx[i]];
}

void
CollocPenta15::GetPositionStiffness(sp_grad::index_type idx, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(idx >= 0);
     ASSERT(idx < iNumEvalPointsStiffness);

     using namespace sp_grad;

     const index_type i = idx % M;
     const index_type j = idx / M;

     ASSERT(i >= 0);
     ASSERT(i < M);
     ASSERT(j >= 0);
     ASSERT(j < N);

     static_assert(sizeof(ri) / sizeof(ri[0]) == M);
     static_assert(sizeof(si) / sizeof(si[0]) == M);
     static_assert(sizeof(ti) / sizeof(ti[0]) == N);

     r(1) = ri[i];
     r(2) = si[i];
     r(3) = ti[j];
}

doublereal
CollocPenta15::dGetWeightStiffness(sp_grad::index_type idx)
{
     ASSERT(idx >= 0);
     ASSERT(idx < iNumEvalPointsStiffness);

     using namespace sp_grad;

     const index_type i = idx % M;
     const index_type j = idx / M;

     ASSERT(i >= 0);
     ASSERT(i < M);
     ASSERT(j >= 0);
     ASSERT(j < N);

     static_assert(sizeof(wi) / sizeof(wi[0]) == M);
     static_assert(sizeof(alphai) / sizeof(alphai[0]) == N);

     return 0.5 * wi[i] * alphai[j];
}

void
CollocPenta15::GetPositionMassLumped(sp_grad::index_type idx, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(idx >= 0);
     ASSERT(idx < iNumEvalPointsMassLumped);

     using namespace sp_grad;

     const index_type i = idx % M_lumped;
     const index_type j = idx / M_lumped;

     ASSERT(i >= 0);
     ASSERT(i < M_lumped);
     ASSERT(j >= 0);
     ASSERT(j < N);

     static_assert(sizeof(ri_lumped) / sizeof(ri_lumped[0]) == M_lumped);
     static_assert(sizeof(si_lumped) / sizeof(si_lumped[0]) == M_lumped);
     static_assert(sizeof(ti_lumped) / sizeof(ti_lumped[0]) == N);

     r(1) = ri_lumped[i];
     r(2) = si_lumped[i];
     r(3) = ti_lumped[j];
}

doublereal
CollocPenta15::dGetWeightMassLumped(sp_grad::index_type idx)
{
     ASSERT(idx >= 0);
     ASSERT(idx < iNumEvalPointsMassLumped);

     using namespace sp_grad;

     const index_type i = idx % M_lumped;
     const index_type j = idx / M_lumped;

     ASSERT(i >= 0);
     ASSERT(i < M_lumped);
     ASSERT(j >= 0);
     ASSERT(j < N);

     static_assert(sizeof(wi_lumped) / sizeof(wi_lumped[0]) == M_lumped);
     static_assert(sizeof(alphai_lumped) / sizeof(alphai_lumped[0]) == N);

     return 0.5 * wi_lumped[i] * alphai_lumped[j];
}



void
CollocTet10h::GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     static_assert(sizeof(r1) / sizeof(r1[0]) == N1);
     static_assert(sizeof(s1) / sizeof(s1[0]) == N1);
     static_assert(sizeof(t1) / sizeof(t1[0]) == N1);

     static_assert(iNumEvalPointsStiffness == N1);

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);

     r(1) = r1[i];
     r(2) = s1[i];
     r(3) = t1[i];
}

doublereal
CollocTet10h::dGetWeightStiffness(sp_grad::index_type i)
{
     static_assert(sizeof(w1) / sizeof(w1[0]) == N1);

     static_assert(iNumEvalPointsStiffness == N1);

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);

     return w1[i];
}

void
CollocTet10h::GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     static_assert(sizeof(r2) / sizeof(r2[0]) == N2);
     static_assert(sizeof(s2) / sizeof(s2[0]) == N2);
     static_assert(sizeof(t2) / sizeof(t2[0]) == N2);
     static_assert(iNumEvalPointsMass == N2);

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMass);

     r(1) = r2[i];
     r(2) = s2[i];
     r(3) = t2[i];
}

doublereal
CollocTet10h::dGetWeightMass(sp_grad::index_type i)
{
     static_assert(sizeof(w2) / sizeof(w2[0]) == N2);
     static_assert(iNumEvalPointsMass == N2);

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMass);

     return w2[i];
}

void
CollocTet10h::GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

doublereal
CollocTet10h::dGetWeightMassLumped(sp_grad::index_type i)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

#endif
