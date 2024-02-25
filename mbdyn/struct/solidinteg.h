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
#include <constexpr_math.h>

#include "sp_matrix_base.h"
#include "solid.h"

struct Gauss2_1D {
     static constexpr sp_grad::index_type iGaussOrder = 2;
     static constexpr doublereal ri[] = {constexpr_math::sqrt(1./3.), -constexpr_math::sqrt(1./3.)};
     static constexpr doublereal alphai[] = {1.0, 1.0};
};

struct Gauss3_1D {
     static constexpr sp_grad::index_type iGaussOrder = 3;
     static constexpr doublereal ri[] = {constexpr_math::sqrt(3./5.), 0., -constexpr_math::sqrt(3./5.)};
     static constexpr doublereal alphai[] = {5./9., 8./9., 5./9.};
};

struct Gauss2Lumped_1D {
     static constexpr sp_grad::index_type iGaussOrder = 2;
     static constexpr doublereal ri[] = {1., -1.};
     static constexpr doublereal alphai[] = {1.0, 1.0};
};

struct Gauss3Lumped_1D {
     static constexpr sp_grad::index_type iGaussOrder = 3;
     static constexpr doublereal ri[] = {1., 0., -1.};
     static constexpr doublereal alphai[] = {2./3., 2./3., 2./3.};
};

struct IntegLayout2_2D {
     static constexpr sp_grad::index_type ridx[] = {0, 0, 1, 1};
     static constexpr sp_grad::index_type sidx[] = {0, 1, 0, 1};
};

struct IntegLayout3_2D {
     static constexpr sp_grad::index_type ridx[] = {0, 0, 0, 1, 1, 1, 2, 2, 2};
     static constexpr sp_grad::index_type sidx[] = {0, 1, 2, 0, 1, 2, 0, 1, 2};
};

struct IntegLayout2_3D {
     static constexpr sp_grad::index_type ridx[] = {0, 0, 0, 0, 1, 1, 1, 1};
     static constexpr sp_grad::index_type sidx[] = {0, 0, 1, 1, 0, 0, 1, 1};
     static constexpr sp_grad::index_type tidx[] = {0, 1, 0, 1, 0, 1, 0, 1};
};

struct IntegLayout3_3D {
     static constexpr sp_grad::index_type ridx[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2};
     static constexpr sp_grad::index_type sidx[] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2};
     static constexpr sp_grad::index_type tidx[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
};

template <typename GaussType_1D, typename IntegLayout_2D>
class Gauss_2D: private GaussType_1D, private IntegLayout_2D {
private:
     using GaussType_1D::iGaussOrder;
     using GaussType_1D::ri;
     using GaussType_1D::alphai;
     using IntegLayout_2D::ridx;
     using IntegLayout_2D::sidx;
public:
     static constexpr sp_grad::index_type iNumEvalPoints = iGaussOrder * iGaussOrder;

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);
};

class Gauss2x2: public Gauss_2D<Gauss2_1D, IntegLayout2_2D> {};
class Gauss3x3: public Gauss_2D<Gauss3_1D, IntegLayout3_2D> {};

class Gauss2x2Lumped: public Gauss_2D<Gauss2Lumped_1D, IntegLayout2_2D> {};
class Gauss3x3Lumped: public Gauss_2D<Gauss3Lumped_1D, IntegLayout3_2D> {};

class CollocTria6h {
public:
     static constexpr sp_grad::index_type iNumEvalPoints = 7;

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);

     // https://zhilin.math.ncsu.edu/
     // Chapter 24: Implementation of Iso-P Triangular Elements
     static constexpr doublereal A = (constexpr_math::sqrt(15.) + 6.) / 2.1E+1;
     static constexpr doublereal B = (6. - constexpr_math::sqrt(15.)) / 2.1E+1;
     static constexpr doublereal P1 = (constexpr_math::sqrt(15.) + 155.) / 2.4E+3;
     static constexpr doublereal P2 = (155. - constexpr_math::sqrt(15.)) / 2.4E+3;
     static constexpr doublereal zeta[7] = {1./3., A, 1. - 2. * A, A, B, 1. - 2. * B, B};
     static constexpr doublereal eta[7]  = {1./3., A, A, 1. - 2. * A, B, B, 1. - 2. * B};
     static constexpr doublereal w[7] = {9./80., P1, P1, P1, P2, P2, P2};
};

template <typename GaussType_1D, typename IntegLayout_3D>
class Gauss_3D: private GaussType_1D, private IntegLayout_3D {
     using GaussType_1D::iGaussOrder;
     using GaussType_1D::ri;
     using GaussType_1D::alphai;
     using IntegLayout_3D::ridx;
     using IntegLayout_3D::sidx;
     using IntegLayout_3D::tidx;
public:
     static constexpr sp_grad::index_type iNumEvalPoints = iGaussOrder*iGaussOrder*iGaussOrder;

     static inline void
     GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r);

     static inline doublereal
     dGetWeight(sp_grad::index_type i);
};

template <typename GaussTypeStiffness_1D, typename GaussTypeMass_1D, typename GaussTypeMassLumped_1D, typename IntegLayoutStiffness_3D, typename IntegLayoutMass_3D, typename IntegLayoutMassLumped_3D>
class GaussSolidStruct_3D {
     typedef Gauss_3D<GaussTypeStiffness_1D, IntegLayoutStiffness_3D> Stiffness;
     typedef Gauss_3D<GaussTypeMass_1D, IntegLayoutMass_3D> Mass;
     typedef Gauss_3D<GaussTypeMassLumped_1D, IntegLayoutMassLumped_3D> MassLumped;
public:
     static constexpr sp_grad::index_type iNumEvalPointsStiffness = Stiffness::iNumEvalPoints;
     static constexpr sp_grad::index_type iNumEvalPointsMass = Mass::iNumEvalPoints;
     static constexpr sp_grad::index_type iNumEvalPointsMassLumped = MassLumped::iNumEvalPoints;

     static inline void
     GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r) {
          Stiffness::GetPosition(i, r);
     }

     static inline doublereal
     dGetWeightStiffness(sp_grad::index_type i) {
          return Stiffness::dGetWeight(i);
     }

     static inline void
     GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r) {
          Mass::GetPosition(i, r);
     }

     static inline doublereal
     dGetWeightMass(sp_grad::index_type i) {
          return Mass::dGetWeight(i);
     }

     static inline void
     GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r) {
          return MassLumped::GetPosition(i, r);
     }

     static inline doublereal
     dGetWeightMassLumped(sp_grad::index_type i) {
          return MassLumped::dGetWeight(i);
     }
};

class Gauss2x2x2: public GaussSolidStruct_3D<Gauss2_1D, Gauss2_1D, Gauss2Lumped_1D, IntegLayout2_3D, IntegLayout2_3D, IntegLayout2_3D> {};
class Gauss3x3x3: public GaussSolidStruct_3D<Gauss3_1D, Gauss3_1D, Gauss3Lumped_1D, IntegLayout3_3D, IntegLayout3_3D, IntegLayout3_3D> {};
class GaussH20r: public GaussSolidStruct_3D<Gauss2_1D, Gauss3_1D, Gauss3Lumped_1D, IntegLayout2_3D, IntegLayout3_3D, IntegLayout3_3D> {};

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

     static constexpr doublereal a1 = 0.25;
     static constexpr doublereal b1 = 1. / 6.;
     static constexpr doublereal c1 = 0.5;
     static constexpr doublereal d1 = -2. / 15.;
     static constexpr doublereal e1 = 3. / 40.;
     static constexpr sp_grad::index_type N1 = iNumEvalPointsStiffness;
     static constexpr doublereal r1[] = {a1, b1, b1, b1, c1};
     static constexpr doublereal s1[] = {a1, b1, b1, c1, b1};
     static constexpr doublereal t1[] = {a1, b1, c1, b1, b1};
     static constexpr doublereal w1[] = {d1, e1, e1, e1, e1};

     static constexpr doublereal a2 = 0.25;
     static constexpr doublereal b2_1 = (7. + constexpr_math::sqrt(15.)) / 34.;
     static constexpr doublereal b2_2 = (7. - constexpr_math::sqrt(15.)) / 34.;
     static constexpr doublereal c2_1 = (13. - 3. * constexpr_math::sqrt(15.)) / 34.;
     static constexpr doublereal c2_2 = (13. + 3. * constexpr_math::sqrt(15.)) / 34.;
     static constexpr doublereal d2 = (5. - constexpr_math::sqrt(15.)) / 20.;
     static constexpr doublereal e2 = (5. + constexpr_math::sqrt(15.)) / 20.;
     static constexpr doublereal f2 = 8. / 405.;
     static constexpr doublereal g2 = (2665. - 14. * constexpr_math::sqrt(15.)) / 226800.;
     static constexpr doublereal h2 = (2665. + 14. * constexpr_math::sqrt(15.)) / 226800.;
     static constexpr doublereal i2 = 5. / 567.;
     static constexpr sp_grad::index_type N2 = iNumEvalPointsMass;

     static constexpr doublereal r2[] = {a2, b2_1, b2_1, c2_1, b2_1, b2_2, b2_2, c2_2, b2_2, d2, e2, d2, e2, d2, e2};
     static constexpr doublereal s2[] = {a2, b2_1, c2_1, b2_1, b2_1, b2_2, c2_2, b2_2, b2_2, e2, d2, d2, e2, e2, d2};
     static constexpr doublereal t2[] = {a2, b2_1, b2_1, b2_1, c2_1, b2_2, b2_2, b2_2, c2_2, d2, d2, e2, d2, e2, e2};
     static constexpr doublereal w2[] = {f2, g2, g2, g2, g2, h2, h2, h2, h2, i2, i2, i2, i2, i2, i2};

     static constexpr doublereal a3 = (5. - constexpr_math::sqrt(5.)) / 20.;
     static constexpr doublereal b3 = (5. + 3. * constexpr_math::sqrt(5.)) / 20.;
     static constexpr doublereal c3 = 1. / 24.;
     static constexpr sp_grad::index_type N3 = 4;

     static constexpr doublereal r3[] = {a3, a3, a3, b3};
     static constexpr doublereal s3[] = {a3, a3, b3, a3};
     static constexpr doublereal t3[] = {a3, b3, a3, a3};
     static constexpr doublereal w3[] = {c3, c3, c3, c3};
};

template <typename GaussType_1D, typename IntegLayout_2D>
void
Gauss_2D<GaussType_1D, IntegLayout_2D>::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);

     static_assert(sizeof(ri) / sizeof(ri[0]) == iGaussOrder, "invalid array size");
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints, "invalid array size");
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints, "invalid array size");

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
}

template <typename GaussType_1D, typename IntegLayout_2D>
doublereal
Gauss_2D<GaussType_1D, IntegLayout_2D>::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);

     static_assert(sizeof(alphai) / sizeof(alphai[0]) == iGaussOrder, "invalid array size");
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints, "invalid array size");
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints, "invalid array size");

     return alphai[ridx[i]] * alphai[sidx[i]];
}

void
CollocTria6h::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);

     static_assert(sizeof(zeta) / sizeof(zeta[0]) == iNumEvalPoints, "invalid array size");
     static_assert(sizeof(eta) / sizeof(eta[0]) == iNumEvalPoints, "invalid array size");

     r(1) = zeta[i];
     r(2) = eta[i];
}

doublereal
CollocTria6h::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);

     static_assert(sizeof(w) / sizeof(w[0]) == iNumEvalPoints, "invalid array size");

     return w[i];
}

template <typename GaussType_1D, typename IntegLayout_3D>
void
Gauss_3D<GaussType_1D, IntegLayout_3D>::GetPosition(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(ri) / sizeof(ri[0]) == iGaussOrder, "invalid array size");
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints, "invalid array size");
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints, "invalid array size");
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPoints, "invalid arrray size");

     r(1) = ri[ridx[i]];
     r(2) = ri[sidx[i]];
     r(3) = ri[tidx[i]];
}

template <typename GaussType_1D, typename IntegLayout_3D>
doublereal
Gauss_3D<GaussType_1D, IntegLayout_3D>::dGetWeight(sp_grad::index_type i)
{
     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPoints);
     ASSERT(ridx[i] >= 0);
     ASSERT(ridx[i] < iGaussOrder);
     ASSERT(sidx[i] >= 0);
     ASSERT(sidx[i] < iGaussOrder);
     ASSERT(tidx[i] >= 0);
     ASSERT(tidx[i] < iGaussOrder);

     static_assert(sizeof(alphai) / sizeof(alphai[0]) == iGaussOrder, "invalid array size");
     static_assert(sizeof(ridx) / sizeof(ridx[0]) == iNumEvalPoints, "invalid array size");
     static_assert(sizeof(sidx) / sizeof(sidx[0]) == iNumEvalPoints, "invalid array size");
     static_assert(sizeof(tidx) / sizeof(tidx[0]) == iNumEvalPoints, "invalid array size");

     return alphai[ridx[i]] * alphai[sidx[i]] * alphai[tidx[i]];
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

     static_assert(sizeof(ri) / sizeof(ri[0]) == M, "invalid array size");
     static_assert(sizeof(si) / sizeof(si[0]) == M, "invalid array size");
     static_assert(sizeof(ti) / sizeof(ti[0]) == N, "invalid array size");

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

     static_assert(sizeof(wi) / sizeof(wi[0]) == M, "invalid array size");
     static_assert(sizeof(alphai) / sizeof(alphai[0]) == N, "invalid array size");

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

     static_assert(sizeof(ri_lumped) / sizeof(ri_lumped[0]) == M_lumped, "invalid array size");
     static_assert(sizeof(si_lumped) / sizeof(si_lumped[0]) == M_lumped, "invalid array size");
     static_assert(sizeof(ti_lumped) / sizeof(ti_lumped[0]) == N, "invalid array size");

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

     static_assert(sizeof(wi_lumped) / sizeof(wi_lumped[0]) == M_lumped, "invalid array size");
     static_assert(sizeof(alphai_lumped) / sizeof(alphai_lumped[0]) == N, "invalid array size");

     return 0.5 * wi_lumped[i] * alphai_lumped[j];
}



void
CollocTet10h::GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     static_assert(sizeof(r1) / sizeof(r1[0]) == N1, "invalid array size");
     static_assert(sizeof(s1) / sizeof(s1[0]) == N1, "invalid array size");
     static_assert(sizeof(t1) / sizeof(t1[0]) == N1, "invalid array size");

     static_assert(iNumEvalPointsStiffness == N1, "size does not match");

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);

     r(1) = r1[i];
     r(2) = s1[i];
     r(3) = t1[i];
}

doublereal
CollocTet10h::dGetWeightStiffness(sp_grad::index_type i)
{
     static_assert(sizeof(w1) / sizeof(w1[0]) == N1, "invalid array size");

     static_assert(iNumEvalPointsStiffness == N1, "invalid parameter");

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);

     return w1[i];
}

void
CollocTet10h::GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     static_assert(sizeof(r2) / sizeof(r2[0]) == N2, "invalid array size");
     static_assert(sizeof(s2) / sizeof(s2[0]) == N2, "invalid array size");
     static_assert(sizeof(t2) / sizeof(t2[0]) == N2, "invalid array size");
     static_assert(iNumEvalPointsMass == N2, "invalid parameter");

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMass);

     r(1) = r2[i];
     r(2) = s2[i];
     r(3) = t2[i];
}

doublereal
CollocTet10h::dGetWeightMass(sp_grad::index_type i)
{
     static_assert(sizeof(w2) / sizeof(w2[0]) == N2, "invalid array size");
     static_assert(iNumEvalPointsMass == N2, "invalid parameter");

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
