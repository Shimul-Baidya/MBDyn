/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <ac/lapack.h>
#include <type_traits>
#include "solid.h"
#include "dataman.h"

class NeoHookean: public ConstitutiveLaw<Vec6, Mat6x6> {
public:
     NeoHookean(const doublereal mu, const doublereal lambda)
          :mu(mu), lambda(lambda) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
          return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy() const override {
          NeoHookean* pCL;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 NeoHookean,
                                 NeoHookean(mu, lambda));
          return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
                         sp_grad::SpColVector<doublereal, iDim>& FTmp) override {
          UpdateElasticTpl(Eps, FTmp);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
                         sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& FTmp) override {
          UpdateElasticTpl(Eps, FTmp);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
                         sp_grad::SpColVector<sp_grad::SpGradient, iDim>& FTmp) override {
          UpdateElasticTpl(Eps, FTmp);
     }

     virtual void
     Update(const Vec6& Eps, const Vec6& EpsPrime) override {
          ConstitutiveLaw<Vec6, Mat6x6>::UpdateElasticSparse(this, Eps);
     }

     template <typename VectorType>
     void UpdateElasticTpl(const VectorType& epsilon, VectorType& sigma) {

          typedef typename VectorType::ValueType T;
          using std::pow;
          using namespace sp_grad;

          SpGradExpDofMapHelper<typename VectorType::ValueType> oDofMap;

          oDofMap.GetDofStat(epsilon);
          oDofMap.Reset();
          oDofMap.InsertDof(epsilon);
          oDofMap.InsertDone();

          // Based on Lars Kuebler 2005, chapter 2.2.1.3, page 25-26

          const SpMatrix<T, 3, 3> C{T{2. * epsilon(1) + 1.},        epsilon(4),           epsilon(6),
                                    epsilon(4), T{2. * epsilon(2) + 1.},          epsilon(5),
                                    epsilon(6),          epsilon(5), T{2. * epsilon(3) + 1.}};

          const SpMatrix<T, 3, 3> CC(C * C, oDofMap);

          T IC, IIC, IIIC;

          oDofMap.MapAssign(IC, C(1, 1) + C(2, 2) + C(3, 3));
          oDofMap.MapAssign(IIC, 0.5 * (IC * IC - (CC(1, 1) + CC(2, 2) + CC(3, 3))));

          Det(C, IIIC, oDofMap);

          T gamma;

          oDofMap.MapAssign(gamma, (lambda * (IIIC - sqrt(IIIC)) - mu) / IIIC);

          constexpr index_type i1[] = {1, 2, 3, 1, 2, 3};
          constexpr index_type i2[] = {1, 2, 3, 2, 3, 1};

          for (index_type i = 1; i <= 6; ++i) {
               const index_type j = i1[i - 1];
               const index_type k = i2[i - 1];

               oDofMap.MapAssign(sigma(i), mu * (j == k) + (CC(j, k) - C(j, k) * IC + IIC * (j == k)) * gamma);
          }
     }
private:
     const doublereal mu, lambda;
};

struct NeoHookeanRead: ConstitutiveLawRead<Vec6, Mat6x6> {
        virtual ConstitutiveLaw<Vec6, Mat6x6>*
        Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
             if (!HP.IsKeyWord("E")) {
                  silent_cerr("keyword \"E\" expected at line " << HP.GetLineData() << "\n");
                  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
             }
             const doublereal E = HP.GetReal();

             if (E <= 0.) {
                  silent_cerr("E must be greater than zero at line " << HP.GetLineData() << "\n");
                  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
             }

             if (!HP.IsKeyWord("nu")) {
                  silent_cerr("keyword \"nu\" expected at line " << HP.GetLineData() << "\n");
                  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
             }

             const doublereal nu = HP.GetReal();

             if (nu < 0. || nu >= 0.5) {
                  // FIXME: incompressible case not implemented yet
                  silent_cerr("nu must be between 0 and 0.5 at line " << HP.GetLineData() << "\n");
                  throw ErrGeneric(MBDYN_EXCEPT_ARGS);
             }

             NeoHookean* pCL = nullptr;

             const doublereal mu = (E / (2. * (1. + nu))); // Lame parameters
             const doublereal lambda = (E * nu / ((1. + nu) * (1. - 2. * nu)));

             SAFENEWWITHCONSTRUCTOR(pCL,
                                    NeoHookean,
                                    NeoHookean(mu, lambda));

             CLType = ConstLawType::ELASTIC;

             return pCL;
        }
};

class Hexahedron8 {
public:
     static const char* ElementName() {
          return "hexahedron8";
     }

     static constexpr sp_grad::index_type iNumNodes = 8;
     static constexpr sp_grad::index_type iNumNodesExtrap = iNumNodes;
     static constexpr bool bHaveDiagMass = true;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumNodesExtrap, iNumComp>& taune);
};

class Hexahedron20 {
public:
     static const char* ElementName() {
          return "hexahedron20";
     }

     static constexpr sp_grad::index_type iNumNodes = 20;
     static constexpr sp_grad::index_type iNumNodesExtrap = iNumNodes;
     static constexpr bool bHaveDiagMass = false; // not applicable to incomplete high order elements

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);
};

class Hexahedron20r {
public:
     static const char* ElementName() {
          return "hexahedron20r";
     }

     static constexpr sp_grad::index_type iNumNodes = 20;
     static constexpr sp_grad::index_type iNumNodesExtrap = 8;
     static constexpr bool bHaveDiagMass = false; // not applicable to incomplete high order elements

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);
};

class Pentahedron15 {
public:
     static const char* ElementName() {
          return "pentahedron15";
     }

     static constexpr sp_grad::index_type iNumNodes = 15;
     static constexpr sp_grad::index_type iNumNodesExtrap = iNumNodes;
     static constexpr bool bHaveDiagMass = false; // to be checked

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);
};

class Tetrahedron10h {
public:
     static const char* ElementName() {
          return "tetrahedron10";
     }

     static constexpr sp_grad::index_type iNumNodes = 10;
     static constexpr sp_grad::index_type iNumNodesExtrap = 4;
     static constexpr bool bHaveDiagMass = false; // not supported by collocation rules

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);
};

class Gauss2 {
private:
     static constexpr sp_grad::index_type iGaussOrder = 2;
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
     static constexpr doublereal ri[] = {0.577350269189626, -0.577350269189626};
     static constexpr doublereal alphai[] = {1.0, 1.0};
     static constexpr doublereal ri_lumped[] = {1., -1.};
     static constexpr doublereal alphai_lumped[] = {1.0, 1.0};
};

class Gauss3 {
private:
     static constexpr sp_grad::index_type iGaussOrder = 3;
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
     static constexpr doublereal ri[] = {0.774596669241483, 0., -0.774596669241483};
     static constexpr doublereal alphai[] = {0.555555555555556, 0.888888888888889, 0.555555555555556};
     static constexpr doublereal ri_lumped[] = {1., 0., -1.};
     static constexpr doublereal alphai_lumped[] = {2./3., 2./3., 2./3.};
};

class GaussH20r: public Gauss2, public Gauss3 {
public:
     using Gauss2::iNumEvalPointsStiffness;
     using Gauss3::iNumEvalPointsMass;
     using Gauss3::iNumEvalPointsMassLumped;

     using Gauss2::GetPositionStiffness;
     using Gauss2::dGetWeightStiffness;

     using Gauss3::GetPositionMass;
     using Gauss3::dGetWeightMass;

     using Gauss3::GetPositionMassLumped;
     using Gauss3::dGetWeightMassLumped;
};

class GaussP15 {
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

class GaussT10h {
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

struct SolidMaterialData {
     std::unique_ptr<ConstitutiveLaw6D> pCSL;
     doublereal E = 0.;
     doublereal nu = 0.;
     doublereal beta = 0.;
};

class SolidConstLaw6D {
public:
     SolidConstLaw6D& operator=(SolidMaterialData&& oMaterial) {
          pConstLaw = std::move(oMaterial.pCSL);

          ASSERT(pConstLaw != nullptr);
          ASSERT(oMaterial.E == 0.);
          ASSERT(oMaterial.nu == 0.);
          ASSERT(oMaterial.beta == 0.);

          return *this;
     }

     void AfterConvergence() {
          ASSERT(pConstLaw.get() != nullptr);

          pConstLaw->AfterConvergence(pConstLaw->GetEpsilon(), pConstLaw->GetEpsilonPrime());
     }

protected:
     template <typename T>
     void
     UpdateElastic(const sp_grad::SpMatrix<T, 3, 3>& G, sp_grad::SpColVector<T, 6>& sigma) {
          using namespace sp_grad;

          ASSERT(pConstLaw.get() != nullptr);

          const SpColVector<T, 6> eps{G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1)};

          pConstLaw->Update(eps, sigma);
     }

     template <typename T>
     void
     UpdateViscoElastic(const sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpMatrix<T, 3, 3>& GP, sp_grad::SpColVector<T, 6>& sigma) {
          using namespace sp_grad;

          ASSERT(pConstLaw.get() != nullptr);

          const SpColVector<T, 6> eps{G(1, 1), G(2, 2), G(3, 3), 2. * G(1, 2), 2. * G(2, 3), 2. * G(3, 1)};
          const SpColVector<T, 6> epsP{GP(1, 1), GP(2, 2), GP(3, 3), 2. * GP(1, 2), 2. * GP(2, 3), 2. * GP(3, 1)};

          pConstLaw->Update(eps, epsP, sigma);
     }

private:
     std::unique_ptr<ConstitutiveLaw6D> pConstLaw;
};

class SolidElasticConstLaw6D: private SolidConstLaw6D {
public:
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::ELASTIC;

     using SolidConstLaw6D::AfterConvergence;

     SolidElasticConstLaw6D& operator=(SolidMaterialData&& oMaterial) {
          SolidConstLaw6D::operator=(std::move(oMaterial));
          return *this;
     }

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, sp_grad::SpColVector<T, 6>& sigma,
            const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          // FIXME: pass oDofMap?
          UpdateElastic(G, sigma);
     }

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpMatrix<T, 3, 3>& GP, sp_grad::SpColVector<T, 6>& sigma,
            const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) = delete;
};

class SolidViscoElasticConstLaw6D: public SolidConstLaw6D {
public:
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::VISCOELASTIC;

     using SolidConstLaw6D::AfterConvergence;

     SolidViscoElasticConstLaw6D& operator=(SolidMaterialData&& oMaterial) {
          SolidConstLaw6D::operator=(std::move(oMaterial));
          return *this;
     }

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, sp_grad::SpColVector<T, 6>& sigma,
            const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) = delete;

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpMatrix<T, 3, 3>& GP, sp_grad::SpColVector<T, 6>& sigma,
            const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          // FIXME: pass oDofMap?
          UpdateViscoElastic(G, GP, sigma);
     }
};

class IsotropicLinearElastic {
public:
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::ELASTIC;

     IsotropicLinearElastic()
          :mu{0.}, lambda{0.} {
          }

     IsotropicLinearElastic(const SolidMaterialData& mat)
          :mu{mat.E / (2. * (1. + mat.nu))}, // Lame parameters
           lambda{mat.E * mat.nu / ((1. + mat.nu) * (1. - 2. * mat.nu))} {
                ASSERT(mat.pCSL == nullptr);
                ASSERT(mat.beta == 0.);
           }

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, sp_grad::SpColVector<T, 6>& sigma, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          using namespace sp_grad;

          oDofMap.MapAssign(sigma(1), (2. * mu + lambda) * G(1, 1) + lambda * (G(2, 2) + G(3, 3)));
          oDofMap.MapAssign(sigma(2), (2. * mu + lambda) * G(2, 2) + lambda * (G(1, 1) + G(3, 3)));
          oDofMap.MapAssign(sigma(3), (2. * mu + lambda) * G(3, 3) + lambda * (G(1, 1) + G(2, 2)));
          sigma(4) = 2. * mu * G(1, 2);
          sigma(5) = 2. * mu * G(2, 3);
          sigma(6) = 2. * mu * G(3, 1);
     }

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpMatrix<T, 3, 3>& GP, sp_grad::SpColVector<T, 6>& S) = delete;

     void AfterConvergence() {}

protected:
     doublereal mu;
     doublereal lambda;
};

class IsotropicLinearViscoElastic: private IsotropicLinearElastic {
public:
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::VISCOELASTIC;

     using IsotropicLinearElastic::AfterConvergence;

     IsotropicLinearViscoElastic()
          :beta{0.} {
     }

     IsotropicLinearViscoElastic(const SolidMaterialData& mat)
          :IsotropicLinearElastic{mat},
           beta{mat.beta} {
                ASSERT(mat.beta >= 0.);
           }

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, sp_grad::SpColVector<T, 6>& S) = delete;

     template <typename T>
     void
     Update(const sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpMatrix<T, 3, 3>& GP, sp_grad::SpColVector<T, 6>& sigma, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) {
          using namespace sp_grad;

          oDofMap.MapAssign(sigma(1), (2. * mu + lambda) * (G(1, 1) + beta * GP(1, 1)) + lambda * (G(2, 2) + beta * GP(2, 2) + G(3, 3) + beta * GP(3, 3)));
          oDofMap.MapAssign(sigma(2), (2. * mu + lambda) * (G(2, 2) + beta * GP(2, 2)) + lambda * (G(1, 1) + beta * GP(1, 1) + G(3, 3) + beta * GP(3, 3)));
          oDofMap.MapAssign(sigma(3), (2. * mu + lambda) * (G(3, 3) + beta * GP(3, 3)) + lambda * (G(1, 1) + beta * GP(1, 1) + G(2, 2) + beta * GP(2, 2)));
          oDofMap.MapAssign(sigma(4), 2. * mu * (G(1, 2) + beta * GP(1, 2)));
          oDofMap.MapAssign(sigma(5), 2. * mu * (G(2, 3) + beta * GP(2, 3)));
          oDofMap.MapAssign(sigma(6), 2. * mu * (G(3, 1) + beta * GP(3, 1)));
     }

private:
     doublereal beta;
};

template <typename SolidConstLawType, sp_grad::index_type iSize>
class SolidConstLawArray: public std::array<SolidConstLawType, iSize> {
public:
     SolidConstLawArray(std::array<SolidMaterialData, iSize>&& rgMaterialData) {
          using namespace sp_grad;

          for (index_type i = 0; i < iSize; ++i) {
               (*this)[i] = std::move(rgMaterialData[i]);
          }
     }
};

constexpr sp_grad::index_type Hexahedron8::iNumNodes;
constexpr sp_grad::index_type Hexahedron8::iNumNodesExtrap;

void
Hexahedron8::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     static_assert(iNumNodes == 8);

     h0d1(1,1) = ((r(2)+1)*(r(3)+1))/8.0E+0;
     h0d1(1,2) = ((r(1)+1)*(r(3)+1))/8.0E+0;
     h0d1(1,3) = ((r(1)+1)*(r(2)+1))/8.0E+0;
     h0d1(2,1) = -((r(2)+1)*(r(3)+1))/8.0E+0;
     h0d1(2,2) = ((1-r(1))*(r(3)+1))/8.0E+0;
     h0d1(2,3) = ((1-r(1))*(r(2)+1))/8.0E+0;
     h0d1(3,1) = -((1-r(2))*(r(3)+1))/8.0E+0;
     h0d1(3,2) = -((1-r(1))*(r(3)+1))/8.0E+0;
     h0d1(3,3) = ((1-r(1))*(1-r(2)))/8.0E+0;
     h0d1(4,1) = ((1-r(2))*(r(3)+1))/8.0E+0;
     h0d1(4,2) = -((r(1)+1)*(r(3)+1))/8.0E+0;
     h0d1(4,3) = ((r(1)+1)*(1-r(2)))/8.0E+0;
     h0d1(5,1) = ((r(2)+1)*(1-r(3)))/8.0E+0;
     h0d1(5,2) = ((r(1)+1)*(1-r(3)))/8.0E+0;
     h0d1(5,3) = -((r(1)+1)*(r(2)+1))/8.0E+0;
     h0d1(6,1) = -((r(2)+1)*(1-r(3)))/8.0E+0;
     h0d1(6,2) = ((1-r(1))*(1-r(3)))/8.0E+0;
     h0d1(6,3) = -((1-r(1))*(r(2)+1))/8.0E+0;
     h0d1(7,1) = -((1-r(2))*(1-r(3)))/8.0E+0;
     h0d1(7,2) = -((1-r(1))*(1-r(3)))/8.0E+0;
     h0d1(7,3) = -((1-r(1))*(1-r(2)))/8.0E+0;
     h0d1(8,1) = ((1-r(2))*(1-r(3)))/8.0E+0;
     h0d1(8,2) = -((r(1)+1)*(1-r(3)))/8.0E+0;
     h0d1(8,3) = -((r(1)+1)*(1-r(2)))/8.0E+0;
}

void
Hexahedron8::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                           sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     static_assert(iNumNodes == 8);

     h(1) = ((r(1)+1)*(r(2)+1)*(r(3)+1))/8.0E+0;
     h(2) = ((1-r(1))*(r(2)+1)*(r(3)+1))/8.0E+0;
     h(3) = ((1-r(1))*(1-r(2))*(r(3)+1))/8.0E+0;
     h(4) = ((r(1)+1)*(1-r(2))*(r(3)+1))/8.0E+0;
     h(5) = ((r(1)+1)*(r(2)+1)*(1-r(3)))/8.0E+0;
     h(6) = ((1-r(1))*(r(2)+1)*(1-r(3)))/8.0E+0;
     h(7) = ((1-r(1))*(1-r(2))*(1-r(3)))/8.0E+0;
     h(8) = ((r(1)+1)*(1-r(2))*(1-r(3)))/8.0E+0;
}

void
Hexahedron8::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                 sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     ShapeFunction(r, h);
}

template <sp_grad::index_type iNumComp>
void
Hexahedron8::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                const sp_grad::SpMatrix<doublereal, iNumNodesExtrap, iNumComp>& taune)
{
     tauni = taune;
}

constexpr sp_grad::index_type Hexahedron20::iNumNodes;
constexpr sp_grad::index_type Hexahedron20::iNumNodesExtrap;

void
Hexahedron20::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                 sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r1_2 = r1 * r1;
     const doublereal r2_2 = r2 * r2;
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 20);

     h0d1(1,1) = ((r2+1)*(r3+1))/8.0E+0-(((r2+1)*(1-r3_2))/4.0E+0+((1-r2_2)*(r3+1))/4.0E+0-(r1*(r2+1)*(r3+1))/2.0E+0)/2.0E+0;
     h0d1(1,2) = ((r1+1)*(r3+1))/8.0E+0-(((r1+1)*(1-r3_2))/4.0E+0-((r1+1)*r2*(r3+1))/2.0E+0+((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0;
     h0d1(1,3) = ((r1+1)*(r2+1))/8.0E+0-((-((r1+1)*(r2+1)*r3)/2.0E+0)+((r1+1)*(1-r2_2))/4.0E+0+((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0;
     h0d1(2,1) = (-((-((r2+1)*(1-r3_2))/4.0E+0)-((1-r2_2)*(r3+1))/4.0E+0-(r1*(r2+1)*(r3+1))/2.0E+0)/2.0E+0)-((r2+1)*(r3+1))/8.0E+0;
     h0d1(2,2) = ((1-r1)*(r3+1))/8.0E+0-(((1-r1)*(1-r3_2))/4.0E+0-((1-r1)*r2*(r3+1))/2.0E+0+((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0;
     h0d1(2,3) = ((1-r1)*(r2+1))/8.0E+0-((-((1-r1)*(r2+1)*r3)/2.0E+0)+((1-r1)*(1-r2_2))/4.0E+0+((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0;
     h0d1(3,1) = (-((-((1-r2)*(1-r3_2))/4.0E+0)-((1-r2_2)*(r3+1))/4.0E+0-(r1*(1-r2)*(r3+1))/2.0E+0)/2.0E+0)-((1-r2)*(r3+1))/8.0E+0;
     h0d1(3,2) = (-((-((1-r1)*(1-r3_2))/4.0E+0)-((1-r1)*r2*(r3+1))/2.0E+0-((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0)-((1-r1)*(r3+1))/8.0E+0;
     h0d1(3,3) = ((1-r1)*(1-r2))/8.0E+0-((-((1-r1)*(1-r2)*r3)/2.0E+0)+((1-r1)*(1-r2_2))/4.0E+0+((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0;
     h0d1(4,1) = ((1-r2)*(r3+1))/8.0E+0-(((1-r2)*(1-r3_2))/4.0E+0+((1-r2_2)*(r3+1))/4.0E+0-(r1*(1-r2)*(r3+1))/2.0E+0)/2.0E+0;
     h0d1(4,2) = (-((-((r1+1)*(1-r3_2))/4.0E+0)-((r1+1)*r2*(r3+1))/2.0E+0-((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0)-((r1+1)*(r3+1))/8.0E+0;
     h0d1(4,3) = ((r1+1)*(1-r2))/8.0E+0-((-((r1+1)*(1-r2)*r3)/2.0E+0)+((r1+1)*(1-r2_2))/4.0E+0+((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0;
     h0d1(5,1) = ((r2+1)*(1-r3))/8.0E+0-(((r2+1)*(1-r3_2))/4.0E+0+((1-r2_2)*(1-r3))/4.0E+0-(r1*(r2+1)*(1-r3))/2.0E+0)/2.0E+0;
     h0d1(5,2) = ((r1+1)*(1-r3))/8.0E+0-(((r1+1)*(1-r3_2))/4.0E+0-((r1+1)*r2*(1-r3))/2.0E+0+((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0;
     h0d1(5,3) = (-((-((r1+1)*(r2+1)*r3)/2.0E+0)-((r1+1)*(1-r2_2))/4.0E+0-((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0)-((r1+1)*(r2+1))/8.0E+0;
     h0d1(6,1) = (-((-((r2+1)*(1-r3_2))/4.0E+0)-((1-r2_2)*(1-r3))/4.0E+0-(r1*(r2+1)*(1-r3))/2.0E+0)/2.0E+0)-((r2+1)*(1-r3))/8.0E+0;
     h0d1(6,2) = ((1-r1)*(1-r3))/8.0E+0-(((1-r1)*(1-r3_2))/4.0E+0-((1-r1)*r2*(1-r3))/2.0E+0+((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0;
     h0d1(6,3) = (-((-((1-r1)*(r2+1)*r3)/2.0E+0)-((1-r1)*(1-r2_2))/4.0E+0-((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0)-((1-r1)*(r2+1))/8.0E+0;
     h0d1(7,1) = (-((-((1-r2)*(1-r3_2))/4.0E+0)-((1-r2_2)*(1-r3))/4.0E+0-(r1*(1-r2)*(1-r3))/2.0E+0)/2.0E+0)-((1-r2)*(1-r3))/8.0E+0;
     h0d1(7,2) = (-((-((1-r1)*(1-r3_2))/4.0E+0)-((1-r1)*r2*(1-r3))/2.0E+0-((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0)-((1-r1)*(1-r3))/8.0E+0;
     h0d1(7,3) = (-((-((1-r1)*(1-r2)*r3)/2.0E+0)-((1-r1)*(1-r2_2))/4.0E+0-((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0)-((1-r1)*(1-r2))/8.0E+0;
     h0d1(8,1) = ((1-r2)*(1-r3))/8.0E+0-(((1-r2)*(1-r3_2))/4.0E+0+((1-r2_2)*(1-r3))/4.0E+0-(r1*(1-r2)*(1-r3))/2.0E+0)/2.0E+0;
     h0d1(8,2) = (-((-((r1+1)*(1-r3_2))/4.0E+0)-((r1+1)*r2*(1-r3))/2.0E+0-((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0)-((r1+1)*(1-r3))/8.0E+0;
     h0d1(8,3) = (-((-((r1+1)*(1-r2)*r3)/2.0E+0)-((r1+1)*(1-r2_2))/4.0E+0-((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0)-((r1+1)*(1-r2))/8.0E+0;
     h0d1(9,1) = -(r1*(r2+1)*(r3+1))/2.0E+0;
     h0d1(9,2) = ((1-r1_2)*(r3+1))/4.0E+0;
     h0d1(9,3) = ((1-r1_2)*(r2+1))/4.0E+0;
     h0d1(10,1) = -((1-r2_2)*(r3+1))/4.0E+0;
     h0d1(10,2) = -((1-r1)*r2*(r3+1))/2.0E+0;
     h0d1(10,3) = ((1-r1)*(1-r2_2))/4.0E+0;
     h0d1(11,1) = -(r1*(1-r2)*(r3+1))/2.0E+0;
     h0d1(11,2) = -((1-r1_2)*(r3+1))/4.0E+0;
     h0d1(11,3) = ((1-r1_2)*(1-r2))/4.0E+0;
     h0d1(12,1) = ((1-r2_2)*(r3+1))/4.0E+0;
     h0d1(12,2) = -((r1+1)*r2*(r3+1))/2.0E+0;
     h0d1(12,3) = ((r1+1)*(1-r2_2))/4.0E+0;
     h0d1(13,1) = -(r1*(r2+1)*(1-r3))/2.0E+0;
     h0d1(13,2) = ((1-r1_2)*(1-r3))/4.0E+0;
     h0d1(13,3) = -((1-r1_2)*(r2+1))/4.0E+0;
     h0d1(14,1) = -((1-r2_2)*(1-r3))/4.0E+0;
     h0d1(14,2) = -((1-r1)*r2*(1-r3))/2.0E+0;
     h0d1(14,3) = -((1-r1)*(1-r2_2))/4.0E+0;
     h0d1(15,1) = -(r1*(1-r2)*(1-r3))/2.0E+0;
     h0d1(15,2) = -((1-r1_2)*(1-r3))/4.0E+0;
     h0d1(15,3) = -((1-r1_2)*(1-r2))/4.0E+0;
     h0d1(16,1) = ((1-r2_2)*(1-r3))/4.0E+0;
     h0d1(16,2) = -((r1+1)*r2*(1-r3))/2.0E+0;
     h0d1(16,3) = -((r1+1)*(1-r2_2))/4.0E+0;
     h0d1(17,1) = ((r2+1)*(1-r3_2))/4.0E+0;
     h0d1(17,2) = ((r1+1)*(1-r3_2))/4.0E+0;
     h0d1(17,3) = -((r1+1)*(r2+1)*r3)/2.0E+0;
     h0d1(18,1) = -((r2+1)*(1-r3_2))/4.0E+0;
     h0d1(18,2) = ((1-r1)*(1-r3_2))/4.0E+0;
     h0d1(18,3) = -((1-r1)*(r2+1)*r3)/2.0E+0;
     h0d1(19,1) = -((1-r2)*(1-r3_2))/4.0E+0;
     h0d1(19,2) = -((1-r1)*(1-r3_2))/4.0E+0;
     h0d1(19,3) = -((1-r1)*(1-r2)*r3)/2.0E+0;
     h0d1(20,1) = ((1-r2)*(1-r3_2))/4.0E+0;
     h0d1(20,2) = -((r1+1)*(1-r3_2))/4.0E+0;
     h0d1(20,3) = -((r1+1)*(1-r2)*r3)/2.0E+0;
}

void
Hexahedron20::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                            sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r1_2 = r1 * r1;
     const doublereal r2_2 = r2 * r2;
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 20);

     h(1) = ((r1+1)*(r2+1)*(r3+1))/8.0E+0-(((r1+1)*(r2+1)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(r2+1)*(r3+1))/4.0E+0)/2.0E+0;
     h(2) = ((1-r1)*(r2+1)*(r3+1))/8.0E+0-(((1-r1)*(r2+1)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(r2+1)*(r3+1))/4.0E+0)/2.0E+0;
     h(3) = ((1-r1)*(1-r2)*(r3+1))/8.0E+0-(((1-r1)*(1-r2)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(1-r2)*(r3+1))/4.0E+0)/2.0E+0;
     h(4) = ((r1+1)*(1-r2)*(r3+1))/8.0E+0-(((r1+1)*(1-r2)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(1-r2)*(r3+1))/4.0E+0)/2.0E+0;
     h(5) = ((r1+1)*(r2+1)*(1-r3))/8.0E+0-(((r1+1)*(r2+1)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(r2+1)*(1-r3))/4.0E+0)/2.0E+0;
     h(6) = ((1-r1)*(r2+1)*(1-r3))/8.0E+0-(((1-r1)*(r2+1)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(r2+1)*(1-r3))/4.0E+0)/2.0E+0;
     h(7) = ((1-r1)*(1-r2)*(1-r3))/8.0E+0-(((1-r1)*(1-r2)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(1-r2)*(1-r3))/4.0E+0)/2.0E+0;
     h(8) = ((r1+1)*(1-r2)*(1-r3))/8.0E+0-(((r1+1)*(1-r2)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(1-r2)*(1-r3))/4.0E+0)/2.0E+0;
     h(9) = ((1-r1_2)*(r2+1)*(r3+1))/4.0E+0;
     h(10) = ((1-r1)*(1-r2_2)*(r3+1))/4.0E+0;
     h(11) = ((1-r1_2)*(1-r2)*(r3+1))/4.0E+0;
     h(12) = ((r1+1)*(1-r2_2)*(r3+1))/4.0E+0;
     h(13) = ((1-r1_2)*(r2+1)*(1-r3))/4.0E+0;
     h(14) = ((1-r1)*(1-r2_2)*(1-r3))/4.0E+0;
     h(15) = ((1-r1_2)*(1-r2)*(1-r3))/4.0E+0;
     h(16) = ((r1+1)*(1-r2_2)*(1-r3))/4.0E+0;
     h(17) = ((r1+1)*(r2+1)*(1-r3_2))/4.0E+0;
     h(18) = ((1-r1)*(r2+1)*(1-r3_2))/4.0E+0;
     h(19) = ((1-r1)*(1-r2)*(1-r3_2))/4.0E+0;
     h(20) = ((r1+1)*(1-r2)*(1-r3_2))/4.0E+0;
}

void
Hexahedron20::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                  sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     ShapeFunction(r, h);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Hexahedron20::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                 const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune)
{
     static_assert(iNumRhs >= iNumNodes);
     using namespace sp_grad;

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodes; ++i) {
               tauni(i, j) = taune(i, j);
          }
     }
}

void
Hexahedron20r::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                  sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     static_assert(iNumNodes == 20);

     h0d1(1,1) = 1.25E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0)+1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h0d1(1,2) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0);
     h0d1(1,3) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+r2+r1+2.0E+0);
     h0d1(2,1) = (-1.25E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0))-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h0d1(2,2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0);
     h0d1(2,3) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+r2-r1+2.0E+0);
     h0d1(3,1) = (-1.25E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0))-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(3,2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(3,3) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3-r2-r1+2.0E+0);
     h0d1(4,1) = 1.25E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0)+1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(4,2) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(4,3) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3-r2+r1+2.0E+0);
     h0d1(5,1) = 1.25E-1*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0)+1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(5,2) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-1.25E-1*(r1-1.0E+0)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0);
     h0d1(5,3) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(6,1) = (-1.25E-1*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0))-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(6,2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-1.25E-1*((-r1)-1.0E+0)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0);
     h0d1(6,3) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(7,1) = (-1.25E-1*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0))-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(7,2) = 1.25E-1*((-r1)-1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(7,3) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(8,1) = 1.25E-1*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0)+1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(8,2) = 1.25E-1*(r1-1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(8,3) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(9,1) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h0d1(9,2) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r3);
     h0d1(9,3) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
     h0d1(10,1) = 2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(10,2) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(10,3) = -2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(11,1) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(11,2) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r3);
     h0d1(11,3) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
     h0d1(12,1) = -2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(12,2) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(12,3) = -2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(13,1) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(13,2) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r3+1.0E+0);
     h0d1(13,3) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
     h0d1(14,1) = 2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(14,2) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(14,3) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(15,1) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(15,2) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r3+1.0E+0);
     h0d1(15,3) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
     h0d1(16,1) = -2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(16,2) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(16,3) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(17,1) = -2.5E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(17,2) = -2.5E-1*(1.0E+0-r1)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(17,3) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(18,1) = 2.5E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(18,2) = -2.5E-1*(r1+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(18,3) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(19,1) = 2.5E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(19,2) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(19,3) = 2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(20,1) = -2.5E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(20,2) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(20,3) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
}

void
Hexahedron20r::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                             sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     static_assert(iNumNodes == 20);

     h(1) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0);
     h(2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0);
     h(3) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0);
     h(4) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0);
     h(5) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0);
     h(6) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0);
     h(7) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0);
     h(8) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0);
     h(9) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h(10) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h(11) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h(12) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h(13) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h(14) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h(15) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h(16) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h(17) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h(18) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h(19) = 2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h(20) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
}

void
Hexahedron20r::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                   sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     static_assert(iNumNodesExtrap == 8);

     h(1) = 1.25E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3);
     h(2) = 1.25E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h(3) = 1.25E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h(4) = 1.25E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3);
     h(5) = 1.25E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0);
     h(6) = 1.25E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h(7) = 1.25E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h(8) = 1.25E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Hexahedron20r::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                  const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune)
{
     using namespace sp_grad;

     static_assert(iNumRhs >= iNumNodesExtrap);

     static constexpr struct {
          index_type ico1, ico2, imid;
     } idxint[] = {
          {1, 2,  9},
          {2, 3, 10},
          {3, 4, 11},
          {4, 1, 12},
          {5, 6, 13},
          {6, 7, 14},
          {7, 8, 15},
          {8, 5, 16},
          {1, 5, 17},
          {2, 6, 18},
          {3, 7, 19},
          {4, 8, 20}
     };

     constexpr index_type iNumNodesInterp = sizeof(idxint) / sizeof(idxint[0]);

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodesExtrap; ++i) {
               tauni(i, j) = taune(i, j);
          }
          for (index_type i = 0; i < iNumNodesInterp; ++i) {
               tauni(idxint[i].imid, j) = 0.5 * (taune(idxint[i].ico1, j)
                                              +  taune(idxint[i].ico2, j));
          }
     }
}

constexpr sp_grad::index_type Pentahedron15::iNumNodes;
constexpr sp_grad::index_type Pentahedron15::iNumNodesExtrap;

void
Pentahedron15::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                  sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 15);

     h0d1(1,1) = (-((1-r3)*((-r3)-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(1-r3);
     h0d1(1,2) = (-((1-r3)*((-r3)-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(1-r3);
     h0d1(1,3) = (-(((-r2)-r1+1)*((-r3)-2*r2-2*r1))/2.0E+0)-(((-r2)-r1+1)*(1-r3))/2.0E+0;
     h0d1(2,1) = ((1-r3)*((-r3)+2*r1-2))/2.0E+0+r1*(1-r3);
     h0d1(2,2) = 0;
     h0d1(2,3) = (-(r1*((-r3)+2*r1-2))/2.0E+0)-(r1*(1-r3))/2.0E+0;
     h0d1(3,1) = 0;
     h0d1(3,2) = ((1-r3)*((-r3)+2*r2-2))/2.0E+0+r2*(1-r3);
     h0d1(3,3) = (-(r2*((-r3)+2*r2-2))/2.0E+0)-(r2*(1-r3))/2.0E+0;
     h0d1(4,1) = (-((r3+1)*(r3-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(r3+1);
     h0d1(4,2) = (-((r3+1)*(r3-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(r3+1);
     h0d1(4,3) = (((-r2)-r1+1)*(r3-2*r2-2*r1))/2.0E+0+(((-r2)-r1+1)*(r3+1))/2.0E+0;
     h0d1(5,1) = ((r3+1)*(r3+2*r1-2))/2.0E+0+r1*(r3+1);
     h0d1(5,2) = 0;
     h0d1(5,3) = (r1*(r3+2*r1-2))/2.0E+0+(r1*(r3+1))/2.0E+0;
     h0d1(6,1) = 0;
     h0d1(6,2) = ((r3+1)*(r3+2*r2-2))/2.0E+0+r2*(r3+1);
     h0d1(6,3) = (r2*(r3+2*r2-2))/2.0E+0+(r2*(r3+1))/2.0E+0;
     h0d1(7,1) = 2*((-r2)-r1+1)*(1-r3)-2*r1*(1-r3);
     h0d1(7,2) = -2*r1*(1-r3);
     h0d1(7,3) = -2*r1*((-r2)-r1+1);
     h0d1(8,1) = 2*r2*(1-r3);
     h0d1(8,2) = 2*r1*(1-r3);
     h0d1(8,3) = -2*r1*r2;
     h0d1(9,1) = -2*r2*(1-r3);
     h0d1(9,2) = 2*((-r2)-r1+1)*(1-r3)-2*r2*(1-r3);
     h0d1(9,3) = -2*((-r2)-r1+1)*r2;
     h0d1(10,1) = 2*((-r2)-r1+1)*(r3+1)-2*r1*(r3+1);
     h0d1(10,2) = -2*r1*(r3+1);
     h0d1(10,3) = 2*r1*((-r2)-r1+1);
     h0d1(11,1) = 2*r2*(r3+1);
     h0d1(11,2) = 2*r1*(r3+1);
     h0d1(11,3) = 2*r1*r2;
     h0d1(12,1) = -2*r2*(r3+1);
     h0d1(12,2) = 2*((-r2)-r1+1)*(r3+1)-2*r2*(r3+1);
     h0d1(12,3) = 2*((-r2)-r1+1)*r2;
     h0d1(13,1) = r3_2-1;
     h0d1(13,2) = r3_2-1;
     h0d1(13,3) = -2*((-r2)-r1+1)*r3;
     h0d1(14,1) = 1-r3_2;
     h0d1(14,2) = 0;
     h0d1(14,3) = -2*r1*r3;
     h0d1(15,1) = 0;
     h0d1(15,2) = 1-r3_2;
     h0d1(15,3) = -2*r2*r3;
}

void
Pentahedron15::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                   sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     ShapeFunction(r, h);
}

void
Pentahedron15::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                             sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 15);

     h(1) = (((-r2)-r1+1)*(1-r3)*((-r3)-2*r2-2*r1))/2.0E+0;
     h(2) = (r1*(1-r3)*((-r3)+2*r1-2))/2.0E+0;
     h(3) = (r2*(1-r3)*((-r3)+2*r2-2))/2.0E+0;
     h(4) = (((-r2)-r1+1)*(r3+1)*(r3-2*r2-2*r1))/2.0E+0;
     h(5) = (r1*(r3+1)*(r3+2*r1-2))/2.0E+0;
     h(6) = (r2*(r3+1)*(r3+2*r2-2))/2.0E+0;
     h(7) = 2*r1*((-r2)-r1+1)*(1-r3);
     h(8) = 2*r1*r2*(1-r3);
     h(9) = 2*((-r2)-r1+1)*r2*(1-r3);
     h(10) = 2*r1*((-r2)-r1+1)*(r3+1);
     h(11) = 2*r1*r2*(r3+1);
     h(12) = 2*((-r2)-r1+1)*r2*(r3+1);
     h(13) = ((-r2)-r1+1)*(1-r3_2);
     h(14) = r1*(1-r3_2);
     h(15) = r2*(1-r3_2);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Pentahedron15::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                  const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune) {
     using namespace sp_grad;

     static_assert(iNumRhs >= iNumNodes);

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodes; ++i) {
               tauni(i, j) = taune(i, j);
          }
     }
}

void
Tetrahedron10h::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                   sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     h0d1(1,1) = 0;
     h0d1(1,2) = 4*r2-1;
     h0d1(1,3) = 0;
     h0d1(2,1) = 0;
     h0d1(2,2) = 0;
     h0d1(2,3) = 4*r3-1;
     h0d1(3,1) = 2*r3-2*((-r3)-r2-r1+1)+2*r2+2*r1-1;
     h0d1(3,2) = 2*r3-2*((-r3)-r2-r1+1)+2*r2+2*r1-1;
     h0d1(3,3) = 2*r3-2*((-r3)-r2-r1+1)+2*r2+2*r1-1;
     h0d1(4,1) = 4*r1-1;
     h0d1(4,2) = 0;
     h0d1(4,3) = 0;
     h0d1(5,1) = 0;
     h0d1(5,2) = 4*r3;
     h0d1(5,3) = 4*r2;
     h0d1(6,1) = -4*r3;
     h0d1(6,2) = -4*r3;
     h0d1(6,3) = 4*((-r3)-r2-r1+1)-4*r3;
     h0d1(7,1) = -4*r2;
     h0d1(7,2) = 4*((-r3)-r2-r1+1)-4*r2;
     h0d1(7,3) = -4*r2;
     h0d1(8,1) = 4*r2;
     h0d1(8,2) = 4*r1;
     h0d1(8,3) = 0;
     h0d1(9,1) = 4*r3;
     h0d1(9,2) = 0;
     h0d1(9,3) = 4*r1;
     h0d1(10,1) = 4*((-r3)-r2-r1+1)-4*r1;
     h0d1(10,2) = -4*r1;
     h0d1(10,3) = -4*r1;
}

void
Tetrahedron10h::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                              sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     h(1) = r2*(2*r2-1);
     h(2) = r3*(2*r3-1);
     h(3) = ((-2*r3)-2*r2-2*r1+1)*((-r3)-r2-r1+1);
     h(4) = r1*(2*r1-1);
     h(5) = 4*r2*r3;
     h(6) = 4*((-r3)-r2-r1+1)*r3;
     h(7) = 4*r2*((-r3)-r2-r1+1);
     h(8) = 4*r1*r2;
     h(9) = 4*r1*r3;
     h(10) = 4*r1*((-r3)-r2-r1+1);
}

void
Tetrahedron10h::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                    sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     h(1) = r(2);
     h(2) = r(3);
     h(3) = 1. - r(1) - r(2) - r(3);
     h(4) = r(1);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Tetrahedron10h::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                   const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune) {
     using namespace sp_grad;

     static_assert(iNumRhs >= iNumNodesExtrap);

     static constexpr struct {
          index_type ico1, ico2, imid;
     } idxint[] = {
          {0, 1, 4},
          {1, 2, 5},
          {2, 0, 6},
          {0, 3, 7},
          {1, 3, 8},
          {2, 3, 9}
     };

     constexpr index_type iNumNodesInterp = sizeof(idxint) / sizeof(idxint[0]);

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodesExtrap; ++i) {
               tauni(i, j) = taune(i, j);
          }
          for (index_type i = 0; i < iNumNodesInterp; ++i) {
               tauni(idxint[i].imid + 1, j) = 0.5 * (taune(idxint[i].ico1 + 1, j)
                                                     +  taune(idxint[i].ico2 + 1, j));
          }
     }
}

constexpr sp_grad::index_type Gauss2::iGaussOrder;
constexpr sp_grad::index_type Gauss2::iNumEvalPointsStiffness;
constexpr sp_grad::index_type Gauss2::iNumEvalPointsMass;
constexpr sp_grad::index_type Gauss2::iNumEvalPointsMassLumped;
constexpr sp_grad::index_type Gauss2::ridx[];
constexpr sp_grad::index_type Gauss2::sidx[];
constexpr sp_grad::index_type Gauss2::tidx[];
constexpr doublereal Gauss2::ri[];
constexpr doublereal Gauss2::alphai[];
constexpr doublereal Gauss2::ri_lumped[];
constexpr doublereal Gauss2::alphai_lumped[];

void
Gauss2::GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
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
Gauss2::dGetWeightStiffness(sp_grad::index_type i)
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
Gauss2::GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
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
Gauss2::dGetWeightMassLumped(sp_grad::index_type i)
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

constexpr sp_grad::index_type Gauss3::iGaussOrder;
constexpr sp_grad::index_type Gauss3::iNumEvalPointsStiffness;
constexpr sp_grad::index_type Gauss3::iNumEvalPointsMass;
constexpr sp_grad::index_type Gauss3::iNumEvalPointsMassLumped;
constexpr sp_grad::index_type Gauss3::ridx[];
constexpr sp_grad::index_type Gauss3::sidx[];
constexpr sp_grad::index_type Gauss3::tidx[];
constexpr doublereal Gauss3::ri[];
constexpr doublereal Gauss3::alphai[];
constexpr doublereal Gauss3::ri_lumped[];
constexpr doublereal Gauss3::alphai_lumped[];

void
Gauss3::GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
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
Gauss3::dGetWeightStiffness(sp_grad::index_type i)
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
Gauss3::GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
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
Gauss3::dGetWeightMassLumped(sp_grad::index_type i)
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

constexpr doublereal GaussP15::ri[];
constexpr doublereal GaussP15::si[];
constexpr doublereal GaussP15::wi[];
constexpr doublereal GaussP15::ti[];
constexpr doublereal GaussP15::alphai[];

constexpr doublereal GaussP15::ri_lumped[];
constexpr doublereal GaussP15::si_lumped[];
constexpr doublereal GaussP15::wi_lumped[];
constexpr doublereal GaussP15::ti_lumped[];
constexpr doublereal GaussP15::alphai_lumped[];

void
GaussP15::GetPositionStiffness(sp_grad::index_type idx, sp_grad::SpColVector<doublereal, 3>& r)
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
GaussP15::dGetWeightStiffness(sp_grad::index_type idx)
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
GaussP15::GetPositionMassLumped(sp_grad::index_type idx, sp_grad::SpColVector<doublereal, 3>& r)
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
GaussP15::dGetWeightMassLumped(sp_grad::index_type idx)
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

constexpr sp_grad::index_type GaussT10h::iNumEvalPointsStiffness;
constexpr sp_grad::index_type GaussT10h::iNumEvalPointsMass;
constexpr sp_grad::index_type GaussT10h::iNumEvalPointsMassLumped;
constexpr doublereal GaussT10h::r1[];
constexpr doublereal GaussT10h::s1[];
constexpr doublereal GaussT10h::t1[];
constexpr doublereal GaussT10h::w1[];

constexpr doublereal GaussT10h::r2[];
constexpr doublereal GaussT10h::s2[];
constexpr doublereal GaussT10h::t2[];
constexpr doublereal GaussT10h::w2[];

constexpr doublereal GaussT10h::r3[];
constexpr doublereal GaussT10h::s3[];
constexpr doublereal GaussT10h::t3[];
constexpr doublereal GaussT10h::w3[];

void
GaussT10h::GetPositionStiffness(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
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
GaussT10h::dGetWeightStiffness(sp_grad::index_type i)
{
     static_assert(sizeof(w1) / sizeof(w1[0]) == N1);

     static_assert(iNumEvalPointsStiffness == N1);

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsStiffness);

     return w1[i];
}

void
GaussT10h::GetPositionMass(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
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
GaussT10h::dGetWeightMass(sp_grad::index_type i)
{
     static_assert(sizeof(w2) / sizeof(w2[0]) == N2);
     static_assert(iNumEvalPointsMass == N2);

     ASSERT(i >= 0);
     ASSERT(i < iNumEvalPointsMass);

     return w2[i];
}

void
GaussT10h::GetPositionMassLumped(sp_grad::index_type i, sp_grad::SpColVector<doublereal, 3>& r)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

doublereal
GaussT10h::dGetWeightMassLumped(sp_grad::index_type i)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType = StructDispNodeAd>
class SolidElemStatic: public SolidElem {
public:
     static constexpr ConstLawType::Type eConstLawType = SolidCSLType::eConstLawType;
     static constexpr sp_grad::index_type iNumNodes = ElementType::iNumNodes;
     static constexpr sp_grad::index_type iNumNodesExtrap = ElementType::iNumNodesExtrap;
     static constexpr sp_grad::index_type iNumEvalPointsStiffness = CollocationType::iNumEvalPointsStiffness;
     static constexpr sp_grad::index_type iNumEvalPointsMassLumped = CollocationType::iNumEvalPointsMassLumped;
     static constexpr sp_grad::index_type iNumDof = iNumNodes * 3;

     SolidElemStatic(unsigned uLabel,
                     const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                     const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                     std::array<SolidMaterialData, iNumEvalPointsStiffness>&& rgMaterialData,
                     flag fOut);
     virtual ~SolidElemStatic();

     virtual void Output(OutputHandler& OH) const override;

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                   doublereal dCoef,
                   enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        enum sp_grad::SpFunctionCall func);

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual VariableSubMatrixHandler&
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

     template <typename T>
     inline void
     InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                   const sp_grad::SpGradientVectorHandler<T>& XCurr,
                   enum sp_grad::SpFunctionCall func);

     virtual VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
                   const VectorHandler& XCurr) override;

     virtual SubVectorHandler&
     InitialAssRes(SubVectorHandler& WorkVec,
                   const VectorHandler& XCurr) override;

     virtual void
     InitialWorkSpaceDim(integer* piNumRows,
                         integer* piNumCols) const override;

     virtual void
     AfterConvergence(const VectorHandler& X,
                      const VectorHandler& XP) override;

     virtual int
     GetNumConnectedNodes() const override;
     
     virtual void
     GetConnectedNodes(std::vector<const Node*>& connectedNodes) const override;
protected:
     template <typename T>
     inline void
     GetNodalDeformations(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                          doublereal dCoef,
                          sp_grad::SpFunctionCall func) const;

     template <typename T>
     inline void
     GetNodalVelocities(sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                        doublereal dCoef,
                        sp_grad::SpFunctionCall func) const;

     inline void
     ComputeStressStrain(sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& epsilon,
                         sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& tau) const;

     template <typename T>
     inline void
     AssStiffnessVecElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                            sp_grad::SpColVector<T, iNumDof>& R,
                            doublereal dCoef,
                            sp_grad::SpFunctionCall func,
                            const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     template <typename T>
     inline void
     AssStiffnessVecViscoElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                 const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                 sp_grad::SpColVector<T, iNumDof>& R,
                                 doublereal dCoef,
                                 sp_grad::SpFunctionCall func,
                                 const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     template <typename T>
     inline void
     AssGravityLoadVec(sp_grad::SpColVector<T, iNumDof>& R,
                       doublereal dCoef,
                       sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssVector(sp_grad::SpGradientAssVec<T>& WorkVec,
               sp_grad::SpColVector<T, iNumDof>& R,
               integer (StructDispNode::*pfnGetFirstIndex)(void) const);

     template <sp_grad::index_type iNumComp>
     inline void
     GaussToNodal(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& taun,
                  const sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, iNumComp>& taug) const;

     struct CollocData {
          static constexpr sp_grad::index_type iNumNodes = SolidElemStatic::iNumNodes;
          static constexpr sp_grad::index_type iNumDof = SolidElemStatic::iNumDof;

          inline void
          Init(sp_grad::index_type iColloc,
               const sp_grad::SpMatrix<doublereal, 3, iNumNodes>& x0,
               SolidMaterialData&& oMaterial,
               const SolidElemStatic* pElem);

          template <typename T>
          inline void
          StrainMatrix(const sp_grad::SpMatrix<T, 3, 3>& F,
                       sp_grad::SpMatrix<T, 6, iNumDof>& BL,
                       const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) const;

          template <typename T>
          inline void
          ComputeStressElastic(const sp_grad::SpMatrix<T, 3, 3>& G,
                               const sp_grad::SpMatrix<T, 3, 3>& F,
                               sp_grad::SpColVector<T, 6>& sigma,
                               const sp_grad::SpGradExpDofMapHelper<T>& oDofMap,
                               const SolidElemStatic* pElem);

          template <typename T>
          inline void
          ComputeStressViscoElastic(const sp_grad::SpMatrix<T, 3, 3>& G,
                                    const sp_grad::SpMatrix<T, 3, 3>& GP,
                                    const sp_grad::SpMatrix<T, 3, 3>& F,
                                    sp_grad::SpColVector<T, 6>& sigma,
                                    const sp_grad::SpGradExpDofMapHelper<T>& oDofMap,
                                    const SolidElemStatic* pElem);

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& G_tmp,
                             const sp_grad::SpColVector<doublereal, 6>& sigma_tmp,
                             const sp_grad::SpMatrix<doublereal, 3, 3>& F_tmp,
                             const SolidElemStatic* pElem);

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& G,
                             const sp_grad::SpColVector<sp_grad::SpGradient, 6>& sigma,
                             const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& F,
                             const SolidElemStatic* pElem) {
          }

          inline void
          UpdateStressStrain(const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& G,
                             const sp_grad::SpColVector<sp_grad::GpGradProd, 6>& sigma,
                             const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& F,
                             const SolidElemStatic* pElem) {
          }

          SolidCSLType oMaterialData;
          sp_grad::SpColVectorA<doublereal, iNumNodes> h;
          sp_grad::SpMatrixA<doublereal, iNumNodes, 3> h0d;
          doublereal detJ;
          sp_grad::SpMatrixA<doublereal, 6, iNumDof> BL0;
          sp_grad::SpMatrixA<doublereal, 3, 3> G;
          sp_grad::SpMatrixA<doublereal, 3, 3> F;
          sp_grad::SpColVectorA<doublereal, 6> sigma;
     };

     sp_grad::SpMatrixA<doublereal, 3, iNumNodes> x0;
     std::array<const StructNodeType*, iNumNodes> rgNodes;
     const sp_grad::SpColVectorA<doublereal, iNumNodes> rhon;
     std::array<CollocData, iNumEvalPointsStiffness> rgCollocData;
};

template <typename ElementType, typename CollocationType, typename SolidCSLType>
class SolidElemDynamic: public SolidElemStatic<ElementType, CollocationType, SolidCSLType, DynamicStructDispNodeAd> {
public:
     typedef SolidElemStatic<ElementType, CollocationType, SolidCSLType, DynamicStructDispNodeAd> SolidElemStaticType;
     using SolidElemStaticType::eConstLawType;
     using SolidElemStaticType::iNumNodes;
     using SolidElemStaticType::iNumEvalPointsStiffness;
     using SolidElemStaticType::iNumEvalPointsMassLumped;
     using SolidElemStaticType::iNumDof;

     SolidElemDynamic(unsigned uLabel,
                      const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                      const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                      std::array<SolidMaterialData, iNumEvalPointsStiffness>&& rgMaterialData,
                      const RigidBodyKinematics* pRBK,
                      flag fOut);
     virtual ~SolidElemDynamic();

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     template <typename T>
     inline void
     AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                   doublereal dCoef,
                   enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        enum sp_grad::SpFunctionCall func);

     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            enum sp_grad::SpFunctionCall func);

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual VariableSubMatrixHandler&
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

private:
     template <typename T>
     inline void
     AssInertiaVec(const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                   sp_grad::SpColVector<T, iNumDof>& R,
                   const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     template <typename T>
     inline void
     AssInertiaVecRBK(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                      sp_grad::SpColVector<T, iNumDof>& R,
                      const sp_grad::SpGradExpDofMapHelper<T>& oDofMap);

     inline void
     AssMassMatrix();

     const RigidBodyKinematics* const pRBK;
     sp_grad::SpMatrixA<doublereal, iNumDof, iNumDof> M;
};

template <typename ElementType, typename CollocationType, typename SolidCSLType>
class SolidElemDynamicDiagMass: public SolidElemDynamic<ElementType, CollocationType, SolidCSLType> {
public:
     typedef SolidElemDynamic<ElementType, CollocationType, SolidCSLType> SolidElemDynamicType;
     using SolidElemDynamicType::eConstLawType;
     using SolidElemDynamicType::iNumNodes;
     using SolidElemDynamicType::iNumEvalPointsStiffness;
     using SolidElemDynamicType::iNumEvalPointsMassLumped;
     using SolidElemDynamicType::iNumDof;

     SolidElemDynamicDiagMass(unsigned uLabel,
                              const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                              const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                              std::array<SolidMaterialData, iNumEvalPointsStiffness>&& rgMaterialData,
                              const RigidBodyKinematics* pRBK,
                              flag fOut);

     virtual ~SolidElemDynamicDiagMass();

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

private:
     inline void
     AssDiagMassMatrix();

     sp_grad::SpColVectorA<doublereal, iNumNodes> diagM;
};

SolidElem::SolidElem(unsigned uLabel,
                     flag fOut)
     :Elem{uLabel, fOut},
      ElemGravityOwner{uLabel, fOut},
      InitialAssemblyElem(uLabel, fOut)
{
}

SolidElem::~SolidElem()
{
}

void
SolidElem::SetValue(DataManager *pDM,
                    VectorHandler& X, VectorHandler& XP,
                    SimulationEntity::Hints *ph)
{
}

std::ostream&
SolidElem::Restart(std::ostream& out) const
{
     out << "## solid element: Restart not implemented yet\n";

     return out;
}

Elem::Type SolidElem::GetElemType() const
{
     return SOLID;
}

bool SolidElem::bIsDeformable() const
{
     return true;
}

unsigned int SolidElem::iGetInitialNumDof() const
{
     return 0u;
}

template <ConstLawType::Type SolidCSLType>
struct SolidElemCSLHelper;

template <>
struct SolidElemCSLHelper<ConstLawType::ELASTIC> {
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::ELASTIC;

     template <typename SolidElementType, typename T>
     static inline void AssRes(SolidElementType* pElem,
                               sp_grad::SpGradientAssVec<T>& WorkVec,
                               doublereal dCoef,
                               enum sp_grad::SpFunctionCall func) {
          pElem->AssResElastic(WorkVec, dCoef, func);
     }
};

template <>
struct SolidElemCSLHelper<ConstLawType::VISCOELASTIC> {
     static constexpr ConstLawType::Type eConstLawType = ConstLawType::VISCOELASTIC;

     template <typename SolidElementType, typename T>
     static inline void AssRes(SolidElementType* pElem,
                               sp_grad::SpGradientAssVec<T>& WorkVec,
                               doublereal dCoef,
                               enum sp_grad::SpFunctionCall func) {
          pElem->AssResViscoElastic(WorkVec, dCoef, func);
     }
};

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumNodes;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumEvalPointsStiffness;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumEvalPointsMassLumped;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
constexpr sp_grad::index_type SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::iNumDof;

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::SolidElemStatic(unsigned uLabel,
                                                                                             const std::array<const StructDispNodeAd*, iNumNodes>& rgNodesDisp,
                                                                                             const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                                                                                             std::array<SolidMaterialData, iNumEvalPointsStiffness>&& rgMaterialData,
                                                                                             flag fOut)
:SolidElem{uLabel, fOut},
 Elem{uLabel, fOut},
 rhon{rhon}
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumNodes; ++i) {
          rgNodes[i - 1] = &dynamic_cast<const StructNodeType&>(*rgNodesDisp[i - 1]);

          const Vec3& x0i = rgNodes[i - 1]->GetXCurr();

          for (index_type j = 1; j <= 3; ++j) {
               x0(j, i) = x0i(j);
          }
     }

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          rgCollocData[iColloc].Init(iColloc, x0, std::move(rgMaterialData[iColloc]), this);
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::~SolidElemStatic()
{
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::Output(OutputHandler& OH) const
{
     using namespace sp_grad;

     if (bToBeOutput()) {
          sp_grad::SpMatrixA<doublereal, iNumEvalPointsStiffness, 6> epsilone;
          sp_grad::SpMatrixA<doublereal, iNumEvalPointsStiffness, 6> taue;
          sp_grad::SpMatrixA<doublereal, iNumNodes, 6> epsilonn;
          sp_grad::SpMatrixA<doublereal, iNumNodes, 6> taun;

          ComputeStressStrain(epsilone, taue);
          GaussToNodal(epsilonn, epsilone);
          GaussToNodal(taun, taue);

          if (OH.UseText(OutputHandler::SOLIDS)) {
               std::ostream& of = OH.Solids();

               of << std::setw(8) << GetLabel();

               for (index_type j = 1; j <= iNumNodes; ++j) {
                    for (index_type i = 1; i <= 6; ++i) {
                         of << ' ' << epsilonn(j, i);
                    }

                    for (index_type i = 1; i <= 6; ++i) {
                         of << ' ' << taun(j, i);
                    }
               }

               of << '\n';
          }
     }
}



template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = iNumDof;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                    doublereal dCoef,
                                                                                    const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                                    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                                                                    enum sp_grad::SpFunctionCall func)
{
     SolidElemCSLHelper<eConstLawType>::AssRes(this, WorkVec, dCoef, func);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialAssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                           const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                                           enum sp_grad::SpFunctionCall func)
{
     return SolidElemCSLHelper<eConstLawType>::AssRes(this, WorkVec, 1., func);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                           doublereal dCoef,
                                                                                           enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1);

     GetNodalDeformations(u, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDone();

     constexpr index_type iNumColsR = iNumDof * iNumEvalPointsStiffness;

     SpColVectorA<T, iNumDof, iNumColsR> R;

     AssStiffnessVecElastic(u, R, dCoef, func, oDofMap);

     if (pGravity) {
          AssGravityLoadVec(R, dCoef, func);
     }

     AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                                doublereal dCoef,
                                                                                                enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1), uP(3, iNumNodes, 1);

     GetNodalDeformations(u, dCoef, func);
     GetNodalVelocities(uP, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(uP);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(uP);
     oDofMap.InsertDone();

     constexpr index_type iNumColsR = iNumDof * iNumEvalPointsStiffness;

     SpColVectorA<T, iNumDof, iNumColsR> R;

     AssStiffnessVecViscoElastic(u, uP, R, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() <= iNumColsR);

     if (pGravity) {
          AssGravityLoadVec(R, dCoef, func);
     }

     ASSERT(R.iGetMaxSize() <= iNumColsR);

     AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);
}

template <typename T>
inline void GreenLagrangeStrain(const sp_grad::SpMatrix<T, 3, 3>& F, sp_grad::SpMatrix<T, 3, 3>& G, const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     // G:0.5 * (transpose(F).F - ident(3));
     // exploit symmetry of the strain tensor

     oDofMap.MapAssign(G(1,1), 5.0E-1*(F(3,1)*F(3,1)+F(2,1)*F(2,1)+F(1,1)*F(1,1)-1));
     oDofMap.MapAssign(G(1,2), 5.0E-1*(F(3,1)*F(3,2)+F(2,1)*F(2,2)+F(1,1)*F(1,2)));
     oDofMap.MapAssign(G(1,3), 5.0E-1*(F(3,1)*F(3,3)+F(2,1)*F(2,3)+F(1,1)*F(1,3)));
     G(2,1) = G(1, 2);
     oDofMap.MapAssign(G(2,2), 5.0E-1*(F(3,2)*F(3,2)+F(2,2)*F(2,2)+F(1,2)*F(1,2)-1));
     oDofMap.MapAssign(G(2,3), 5.0E-1*(F(3,2)*F(3,3)+F(2,2)*F(2,3)+F(1,2)*F(1,3)));
     G(3,1) = G(1, 3);
     G(3,2) = G(2, 3);
     oDofMap.MapAssign(G(3,3), 5.0E-1*(F(3,3)*F(3,3)+F(2,3)*F(2,3)+F(1,3)*F(1,3)-1));
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssStiffnessVecElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                                    sp_grad::SpColVector<T, iNumDof>& R,
                                                                                                    const doublereal dCoef,
                                                                                                    const sp_grad::SpFunctionCall func,
                                                                                                    const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     SpMatrixA<T, 6, iNumDof, 2 * iNumDof> BL;
     SpMatrixA<T, 3, 3, iNumDof> F, G;
     SpColVectorA<T, 6, iNumDof> sigma;
     SpColVectorA<T, iNumDof, iNumDof> dR;

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          const doublereal alpha = CollocationType::dGetWeightStiffness(iColloc);
          const auto& h0d = rgCollocData[iColloc].h0d;

          F.MapAssign(u * h0d, oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               F(i, i) += 1.;
          }

          ASSERT(F.iGetMaxSize() <= iNumDof);

          GreenLagrangeStrain(F, G, oDofMap);

          ASSERT(G.iGetMaxSize() <= iNumDof);

          rgCollocData[iColloc].ComputeStressElastic(G, F, sigma, oDofMap, this);
          sigma *= alpha * rgCollocData[iColloc].detJ;

          ASSERT(sigma.iGetMaxSize() <= iNumDof);

          rgCollocData[iColloc].StrainMatrix(F, BL, oDofMap);

          dR.MapAssign(Transpose(BL) * sigma, oDofMap);

          R -= dR;
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssStiffnessVecViscoElastic(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                                         const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                                                                                         sp_grad::SpColVector<T, iNumDof>& R,
                                                                                                         const doublereal dCoef,
                                                                                                         const sp_grad::SpFunctionCall func,
                                                                                                         const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     SpMatrixA<T, 6, iNumDof, 2 * iNumDof> BL;
     SpColVectorA<T, 6, iNumDof> sigma;
     SpMatrixA<T, 3, 3, iNumDof> F, FP, C, invC, G, FP_Tr_F;
     SpMatrixA<T, 3, 3, 2 * iNumDof> GP_detF, GP_scaled;
     SpColVectorA<T, iNumDof, iNumDof> dR;

     T detC;

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          const auto& h0d = rgCollocData[iColloc].h0d;

          F.MapAssign(u * h0d, oDofMap);

          for (index_type i = 1; i <= 3; ++i) {
               F(i, i) += 1.;
          }

          ASSERT(F.iGetMaxSize() <= iNumDof);

          FP.MapAssign(uP * h0d, oDofMap);

          ASSERT(FP.iGetMaxSize() <= iNumDof);

          C.MapAssign(Transpose(F) * F, oDofMap);

          ASSERT(C.iGetMaxSize() <= iNumDof);

          InvSymm(C, invC, detC, oDofMap);

          ASSERT(invC.iGetMaxSize() <= iNumDof);

          G.MapAssign(0.5 * (C - Eye3), oDofMap);

          ASSERT(G.iGetMaxSize() <= iNumDof);

          FP_Tr_F.MapAssign(Transpose(FP) * F, oDofMap);

          ASSERT(FP_Tr_F.iGetMaxSize() <= iNumDof);

          GP_detF.MapAssign((0.5 * (FP_Tr_F + Transpose(FP_Tr_F)) * Det(F)), oDofMap);

          ASSERT(GP_detF.iGetMaxSize() <= iNumDof);

          GP_scaled.MapAssign(invC * GP_detF * invC, oDofMap); // Lars Kuebler 2005, equation 2.92, page 38

          ASSERT(GP_scaled.iGetMaxSize() <= iNumDof);

          rgCollocData[iColloc].ComputeStressViscoElastic(G, GP_scaled, F, sigma, oDofMap, this);

          ASSERT(sigma.iGetMaxSize() <= iNumDof);

          const doublereal alpha = CollocationType::dGetWeightStiffness(iColloc);

          rgCollocData[iColloc].StrainMatrix(F, BL, oDofMap);

          ASSERT(BL.iGetMaxSize() <= 2 * iNumDof);

          sigma *= alpha * rgCollocData[iColloc].detJ;

          ASSERT(sigma.iGetMaxSize() <= iNumDof);

          dR.MapAssign(Transpose(BL) * sigma, oDofMap);

          R -= dR;
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssGravityLoadVec(sp_grad::SpColVector<T, iNumDof>& R,
                                                                                               doublereal dCoef,
                                                                                               sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 X, fg;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);

          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);

          X = Zero3;

          for (index_type i = 1; i <= iNumNodes; ++i) {
               const Vec3& Xi = rgNodes[i - 1]->GetXCurr();

               for (index_type j = 1; j <= 3; ++j) {
                    X(j) += Xi(j) * h(i);
               }
          }

          bGetGravity(X, fg);

          const doublereal rho = Dot(h, rhon);

          fg *= rho * alpha * Det(J);

          for (index_type i = 1; i <= iNumNodes; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    R((i - 1) * 3 + j) += h(i) * fg(j);
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssVector(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                       sp_grad::SpColVector<T, iNumDof>& R,
                                                                                       integer (StructDispNode::*pfnGetFirstIndex)(void) const)
{
     using namespace sp_grad;

     for (index_type i = 1; i <= iNumNodes; ++i) {
          const index_type iEqIndex = (rgNodes[i - 1]->*pfnGetFirstIndex)();

          for (index_type j = 1; j <= 3; ++j) {
               WorkVec.AddItem(iEqIndex + j, R((i - 1) * 3 + j));
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SubVectorHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssRes(SubVectorHandler& WorkVec,
                                                                                    doublereal dCoef,
                                                                                    const VectorHandler& XCurr,
                                                                                    const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemStatic::AssRes");

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
VariableSubMatrixHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssJac(VariableSubMatrixHandler& WorkMat,
                                                                                    doublereal dCoef,
                                                                                    const VectorHandler& XCurr,
                                                                                    const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemStatic::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AssJac(VectorHandler& JacY,
                                                                                    const VectorHandler& Y,
                                                                                    doublereal dCoef,
                                                                                    const VectorHandler& XCurr,
                                                                                    const VectorHandler& XPrimeCurr,
                                                                                    VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
VariableSubMatrixHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                                                                           const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::InitialAssJac(this,
                                                 WorkMat.SetSparseGradient(),
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_JAC);

     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
SubVectorHandler&
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialAssRes(SubVectorHandler& WorkVec,
                                                                                           const VectorHandler& XCurr)
{
     using namespace sp_grad;

     SpGradientAssVec<doublereal>::InitialAssRes(this,
                                                 WorkVec,
                                                 XCurr,
                                                 sp_grad::INITIAL_ASS_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::InitialWorkSpaceDim(integer* piNumRows,
                                                                                                 integer* piNumCols) const
{
     *piNumRows = iNumDof;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::AfterConvergence(const VectorHandler& X,
                                                                                              const VectorHandler& XP)
{
     for (auto& oColloc: rgCollocData) {
          oColloc.oMaterialData.AfterConvergence();
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
int
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetNumConnectedNodes() const
{
     return iNumNodes;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetConnectedNodes(std::vector<const Node*>& connectedNodes) const
{
     using namespace sp_grad;
     
     connectedNodes.resize(iNumNodes);

     for (index_type i = 0; i < iNumNodes; ++i) {
          connectedNodes[i] = rgNodes[i];
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetNodalDeformations(sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                                                  const doublereal dCoef,
                                                                                                  const sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     SpColVectorA<T, 3, 1> Xj;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          rgNodes[j - 1]->GetXCurr(Xj, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               u(i, j) = std::move(Xj(i));
               u(i, j) -= x0(i, j);
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GetNodalVelocities(sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                                                                                const doublereal dCoef,
                                                                                                const sp_grad::SpFunctionCall func) const
{
     using namespace sp_grad;

     SpColVectorA<T, 3, 1> XPj;

     for (index_type j = 1; j <= iNumNodes; ++j) {
          rgNodes[j - 1]->GetVCurr(XPj, dCoef, func);

          for (index_type i = 1; i <= 3; ++i) {
               uP(i, j) = std::move(XPj(i));
          }
     }
}

namespace {
     template <int N>
     struct log2int {
          static constexpr int value = (N / 2 > 0) ? (log2int<N / 2>::value + 1) : 0;
     };

     template <>
     struct log2int<1> {
          static constexpr int value = 0;
     };

     template <>
     struct log2int<0> {
          static constexpr int value = std::numeric_limits<int>::min();
     };

     static_assert(log2int<0>::value < 0);
     static_assert(log2int<1>::value == 0);
     static_assert(log2int<2>::value == 1);
     static_assert(log2int<4>::value == 2);
     static_assert(log2int<8>::value == 3);
     static_assert(log2int<16>::value == 4);
     static_assert(log2int<1024>::value == 10);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <sp_grad::index_type iNumComp>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::GaussToNodal(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& taun,
                                                                                          const sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, iNumComp>& taug) const
{
     using namespace sp_grad;

     static_assert(iNumNodesExtrap <= iNumNodes);
     static_assert(iNumNodesExtrap <= iNumEvalPointsStiffness);

     SpMatrixA<doublereal, iNumEvalPointsStiffness, iNumNodesExtrap> H;
     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodesExtrap> h;

     for (index_type i = 1; i <= iNumEvalPointsStiffness; ++i) {
          CollocationType::GetPositionStiffness(i - 1, r);
          ElementType::ShapeFunctionExtrap(r, h);

          for (index_type j = 1; j <= iNumNodesExtrap; ++j) {
               H(i, j) = h(j); // FIXME: select only a subset of available nodes
          }
     }

     constexpr integer M = iNumEvalPointsStiffness;
     constexpr integer N = iNumNodesExtrap;
     constexpr integer MINMN = M < N ? M : N;
     constexpr integer NRHS = iNumComp;
     constexpr integer LDB = M;
     constexpr integer SMLSIZ = 25; // FIXME: this value should be determined by configure
     constexpr integer NLVLTMP = log2int<static_cast<integer>(MINMN / (SMLSIZ + 1.0) + 0.5)>::value + 1;
     constexpr integer NLVL = NLVLTMP > 0 ? NLVLTMP : 0;
     constexpr integer LIWORK = 3 * MINMN * NLVL + 11 * MINMN;
     constexpr integer LWORK = M >= N
          ? 12 * N + 2 * N * SMLSIZ + 8 * N * NLVL + N * NRHS + std::pow(SMLSIZ + 1, 2)
          : 12 * M + 2 * M * SMLSIZ + 8 * M * NLVL + M * NRHS + std::pow(SMLSIZ + 1, 2);
     doublereal RCOND = -1., WORK[LWORK], S[MINMN];
     integer RANK, IWORK[LIWORK], INFO;

     SpMatrixA<doublereal, LDB, NRHS> B;

     for (index_type j = 1; j <= NRHS; ++j) {
          for (index_type i = 1; i <= iNumEvalPointsStiffness; ++i) {
               B(i, j) = taug(i, j);
          }
     }

     __FC_DECL__(dgelsd)(&M, &N, &NRHS, H.begin(), &M, B.begin(), &LDB, S, &RCOND, &RANK, WORK, &LWORK, IWORK, &INFO);

     if (INFO != 0) {
          silent_cerr("solid(" << GetLabel() << "): dgelsd failed with status " << INFO << "\n");
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     ElementType::GaussToNodalInterp(taun, B);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::ComputeStressStrain(sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& epsilon,
                                                                                                 sp_grad::SpMatrix<doublereal, iNumEvalPointsStiffness, 6>& tau) const
{
     using namespace sp_grad;

     static constexpr index_type idxr[6] = {1, 2, 3, 1, 2, 3};
     static constexpr index_type idxc[6] = {1, 2, 3, 2, 3, 1};

     SpMatrixA<doublereal, 3, 3> S;

     for (index_type iColloc = 0; iColloc < iNumEvalPointsStiffness; ++iColloc) {
          const auto& sigma = rgCollocData[iColloc].sigma;

          for (index_type i = 1; i <= 6; ++i) {
               S(idxr[i - 1], idxc[i - 1]) = sigma(i);
          }

          for (index_type i = 4; i <= 6; ++i) {
               S(idxc[i - 1], idxr[i - 1]) = sigma(i);
          }

          const auto& F = rgCollocData[iColloc].F;

          const SpMatrixA<doublereal, 3, 3> taui = F * S * Transpose(F) / Det(F);

          const auto& G = rgCollocData[iColloc].G;

          for (index_type i = 1; i <= 3; ++i) {
               epsilon(iColloc + 1, i) = sqrt(1. + 2. * G(i, i)) - 1.;
          }

          for (index_type i = 4; i <= 6; ++i) {
               epsilon(iColloc + 1, i) = 2. * G(idxr[i - 1], idxc[i - 1]) / ((1. + epsilon(iColloc + 1, idxr[i - 1])) * (1. + epsilon(iColloc + 1, idxc[i - 1])));
          }

          for (index_type i = 1; i <= 6; ++i) {
               tau(iColloc + 1, i) = taui(idxr[i - 1], idxc[i - 1]);
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::Init(const sp_grad::index_type iColloc,
                                                                                              const sp_grad::SpMatrix<doublereal, 3, iNumNodes>& x0,
                                                                                              SolidMaterialData&& oMaterial,
                                                                                              const SolidElemStatic* const pElem)
{
     using namespace sp_grad;

     oMaterialData = std::move(oMaterial);

     SpColVectorA<doublereal, 3> r;
     SpMatrixA<doublereal, ElementType::iNumNodes, 3> hd;

     CollocationType::GetPositionStiffness(iColloc, r);
     ElementType::ShapeFunction(r, h);
     ElementType::ShapeFunctionDeriv(r, hd);

     const SpMatrix<doublereal, 3, 3> J = Transpose(x0 * hd);
     SpMatrixA<doublereal, 3, 3> invJ;

     Inv(J, invJ, detJ);

     if (detJ <= 0.) {
          silent_cerr("solid(" << pElem << "): Jacobian is singular: det(J)=" << detJ << "\n");
          throw ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     h0d = hd * Transpose(invJ);

     for (index_type k = 1; k <= iNumNodes; ++k) {
          for (index_type i = 1; i <= 3; ++i) {
               BL0(i, (k - 1) * 3 + i) = h0d(k, i);
          }

          BL0(4, (k - 1) * 3 + 1) = h0d(k, 2);
          BL0(4, (k - 1) * 3 + 2) = h0d(k, 1);
          BL0(5, (k - 1) * 3 + 2) = h0d(k, 3);
          BL0(5, (k - 1) * 3 + 3) = h0d(k, 2);
          BL0(6, (k - 1) * 3 + 1) = h0d(k, 3);
          BL0(6, (k - 1) * 3 + 3) = h0d(k, 1);
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::ComputeStressElastic(const sp_grad::SpMatrix<T, 3, 3>& G,
                                                                                                              const sp_grad::SpMatrix<T, 3, 3>& F,
                                                                                                              sp_grad::SpColVector<T, 6>& sigma,
                                                                                                              const sp_grad::SpGradExpDofMapHelper<T>& oDofMap,
                                                                                                              const SolidElemStatic* const pElem)
{
     oMaterialData.Update(G, sigma, oDofMap);
     UpdateStressStrain(G, sigma, F, pElem);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::ComputeStressViscoElastic(const sp_grad::SpMatrix<T, 3, 3>& G,
                                                                                                                   const sp_grad::SpMatrix<T, 3, 3>& GP,
                                                                                                                   const sp_grad::SpMatrix<T, 3, 3>& F,
                                                                                                                   sp_grad::SpColVector<T, 6>& sigma,
                                                                                                                   const sp_grad::SpGradExpDofMapHelper<T>& oDofMap,
                                                                                                                   const SolidElemStatic* const pElem)
{
     oMaterialData.Update(G, GP, sigma, oDofMap);
     UpdateStressStrain(G, sigma, F, pElem);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
inline void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::UpdateStressStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& G_tmp,
                                                                                                            const sp_grad::SpColVector<doublereal, 6>& sigma_tmp,
                                                                                                            const sp_grad::SpMatrix<doublereal, 3, 3>& F_tmp,
                                                                                                            const SolidElemStatic* const pElem)
{
     using namespace sp_grad;

     G = G_tmp;
     sigma = sigma_tmp;
     F = F_tmp;

     const doublereal detF = Det(F);

     if (detF <= 0.) {
          silent_cerr("solid(" << pElem->GetLabel() << "): element became excessively distorted: Det(F) = " << detF << "\n");
          throw NonlinearSolver::ErrSimulationDiverged(MBDYN_EXCEPT_ARGS);
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType, typename StructNodeType>
template <typename T>
void
SolidElemStatic<ElementType, CollocationType, SolidCSLType, StructNodeType>::CollocData::StrainMatrix(const sp_grad::SpMatrix<T, 3, 3>& F,
                                                                                                      sp_grad::SpMatrix<T, 6, iNumDof>& BL,
                                                                                                      const sp_grad::SpGradExpDofMapHelper<T>& oDofMap) const
{
     using namespace sp_grad;

     for (index_type k = 1; k <= iNumNodes; ++k) {
          for (index_type i = 1; i <= 3; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    oDofMap.MapAssign(BL(i, (k - 1) * 3 + j), (F(j , i) - Eye3(j, i)) * h0d(k, i) + BL0(i, (k - 1) * 3 + j));
               }
          }

          static constexpr index_type idx1[] = {2, 3, 3};
          static constexpr index_type idx2[] = {1, 2, 1};

          for (index_type i = 1; i <= 3; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    oDofMap.MapAssign(BL(i + 3, (k - 1) * 3 + j), (F(j, idx2[i - 1]) - Eye3(j, idx2[i - 1])) * h0d(k, idx1[i - 1]) + (F(j, idx1[i - 1]) - Eye3(j, idx1[i - 1])) * h0d(k, idx2[i - 1]) + BL0(i + 3, (k - 1) * 3 + j));
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::SolidElemDynamic(unsigned uLabel,
                                                                               const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                                                                               const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                                                                               std::array<SolidMaterialData, iNumEvalPointsStiffness>&& rgMaterialData,
                                                                               const RigidBodyKinematics* const pRBK,
                                                                               flag fOut)
:Elem{uLabel, fOut},
 SolidElemStaticType{uLabel, rgNodes, rhon, std::move(rgMaterialData), fOut},
 pRBK(pRBK)
{
     AssMassMatrix();
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::~SolidElemDynamic()
{
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
void SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssMassMatrix()
{
     using namespace sp_grad;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);
          const doublereal detJ = Det(J);

          if (detJ <= 0.) {
               silent_cerr("solid(" << this->GetLabel() << "): Jacobian is singular: det(J) = " << detJ << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal rho = Dot(h, this->rhon); // interpolate from nodes to collocation points
          const doublereal dm = rho * detJ * alpha;

          for (index_type k = 1; k <= iNumNodes; ++k) {
               for (index_type j = 1; j <= iNumNodes; ++j) {
                    const doublereal dmhjhk = dm * h(j) * h(k);

                    for (index_type i = 1; i <= 3; ++i) {
                         M((j - 1) * 3 + i, (k - 1) * 3 + i) += dmhjhk;
                    }
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
void SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 2 * iNumDof;
     *piNumCols = 0;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
template <typename T>
inline void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssResElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                            doublereal dCoef,
                                                                            enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1), uP(3, iNumNodes, 1);

     this->GetNodalDeformations(u, dCoef, func);
     this->GetNodalVelocities(uP, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(uP);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(uP);
     oDofMap.InsertDone();

     const index_type iNumColsR = (1 + 3 * (pRBK != nullptr)) * iNumDof * iNumEvalPointsStiffness;

     SpColVector<T, iNumDof> R(iNumDof, iNumColsR);

     this->AssStiffnessVecElastic(u, R, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() <= iNumColsR);

     if (this->pGravity) {
          this->AssGravityLoadVec(R, dCoef, func);
          ASSERT(R.iGetMaxSize() <= iNumColsR);
     }

     if (pRBK) {
          AssInertiaVecRBK(u, R, oDofMap);
          ASSERT(R.iGetMaxSize() <= iNumColsR);
     }

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);

     AssInertiaVec(uP, R, oDofMap);
     ASSERT(R.iGetMaxSize() <= iNumColsR);

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstPositionIndex);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
template <typename T>
inline void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssResViscoElastic(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                                 doublereal dCoef,
                                                                                 enum sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     SpMatrix<T, 3, iNumNodes> u(3, iNumNodes, 1), uP(3, iNumNodes, 1);

     this->GetNodalDeformations(u, dCoef, func);
     this->GetNodalVelocities(uP, dCoef, func);

     SpGradExpDofMapHelper<T> oDofMap;

     oDofMap.GetDofStat(u);
     oDofMap.GetDofStat(uP);
     oDofMap.Reset();
     oDofMap.InsertDof(u);
     oDofMap.InsertDof(uP);
     oDofMap.InsertDone();

     constexpr index_type iNumEvalPointsR = CollocationType::iNumEvalPointsStiffness > CollocationType::iNumEvalPointsMass
          ? CollocationType::iNumEvalPointsStiffness
          : CollocationType::iNumEvalPointsMass;

     const index_type iNumColsR = (1 + 3 * (pRBK != nullptr)) * iNumDof * iNumEvalPointsR;

     SpColVector<T, iNumDof> R(iNumDof, iNumColsR);

     this->AssStiffnessVecViscoElastic(u, uP, R, dCoef, func, oDofMap);

     ASSERT(R.iGetMaxSize() <= iNumColsR);

     if (this->pGravity) {
          this->AssGravityLoadVec(R, dCoef, func);
          ASSERT(R.iGetMaxSize() <= iNumColsR);
     }

     if (pRBK) {
          AssInertiaVecRBK(u, R, oDofMap);
          ASSERT(R.iGetMaxSize() <= iNumColsR);
     }

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstMomentumIndex);

     AssInertiaVec(uP, R, oDofMap);
     ASSERT(R.iGetMaxSize() <= iNumColsR);

     this->AssVector(WorkVec, R, &StructDispNode::iGetFirstPositionIndex);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
template <typename T>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssInertiaVec(const sp_grad::SpMatrix<T, 3, iNumNodes>& uP,
                                                                            sp_grad::SpColVector<T, iNumDof>& R,
                                                                            const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     SpColVector<T, iNumDof> UP(iNumDof, 0);

     static_assert(iNumDof == iNumNodes * 3);

     for (index_type j = 1; j <= iNumNodes; ++j) {
          for (index_type i = 1; i <= 3; ++i) {
               UP((j - 1) * 3 + i) = uP(i, j);
          }
     }

     R.MapAssign(M * UP, oDofMap);
     R *= -1.;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
template <typename T>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssInertiaVecRBK(const sp_grad::SpMatrix<T, 3, iNumNodes>& u,
                                                                               sp_grad::SpColVector<T, iNumDof>& R,
                                                                               const sp_grad::SpGradExpDofMapHelper<T>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(pRBK != nullptr);

     SpColVectorA<T, 3, iNumNodes> Xc;
     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;
     SpColVectorA<T, 3, iNumNodes> frbk;

     for (index_type iColloc = 0; iColloc < CollocationType::iNumEvalPointsMass; ++iColloc) {
          CollocationType::GetPositionMass(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);

          Xc = Zero3;

          for (index_type j = 1; j <= iNumNodes; ++j) {
               Xc += (u.GetCol(j) + this->x0.GetCol(j)) * h(j);
          }

          const doublereal rho = Dot(h, this->rhon);
          const doublereal alpha = CollocationType::dGetWeightMass(iColloc);
          const doublereal dm = rho * alpha * Det(J);

          frbk.MapAssign((pRBK->GetXPP()
                          + Cross(pRBK->GetWP(), Xc)
                          + Cross(pRBK->GetW(), Cross(pRBK->GetW(), Xc))) * dm, oDofMap);

          for (index_type i = 1; i <= iNumNodes; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    R((i - 1) * 3 + j) -= h(i) * frbk(j); // FIXME: Use oDofMap!
               }
          }
     }
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
template <typename T>
inline void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                                                     doublereal dCoef,
                                                                     const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                                                     const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                                                     enum sp_grad::SpFunctionCall func)
{
     SolidElemCSLHelper<eConstLawType>::AssRes(this, WorkVec, dCoef, func);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
SubVectorHandler&
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssRes(SubVectorHandler& WorkVec,
                                                                     doublereal dCoef,
                                                                     const VectorHandler& XCurr,
                                                                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemDynamic::AssRes");

     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
VariableSubMatrixHandler&
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssJac(VariableSubMatrixHandler& WorkMat,
                                                                     doublereal dCoef,
                                                                     const VectorHandler& XCurr,
                                                                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("SolidElemDynamic::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
void
SolidElemDynamic<ElementType, CollocationType, SolidCSLType>::AssJac(VectorHandler& JacY,
                                                                     const VectorHandler& Y,
                                                                     doublereal dCoef,
                                                                     const VectorHandler& XCurr,
                                                                     const VectorHandler& XPrimeCurr,
                                                                     VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
SolidElemDynamicDiagMass<ElementType, CollocationType, SolidCSLType>::SolidElemDynamicDiagMass(unsigned uLabel,
                                                                                               const std::array<const StructDispNodeAd*, iNumNodes>& rgNodes,
                                                                                               const sp_grad::SpColVector<doublereal, iNumNodes>& rhon,
                                                                                               std::array<SolidMaterialData, iNumEvalPointsStiffness>&& rgMaterialData,
                                                                                               const RigidBodyKinematics* pRBK,
                                                                                               flag fOut)
:Elem{uLabel, fOut},
 SolidElemDynamicType{uLabel, rgNodes, rhon, std::move(rgMaterialData), pRBK, fOut}
{
     AssDiagMassMatrix();
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
SolidElemDynamicDiagMass<ElementType, CollocationType, SolidCSLType>::~SolidElemDynamicDiagMass()
{
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
SubVectorHandler&
SolidElemDynamicDiagMass<ElementType, CollocationType, SolidCSLType>::AssRes(SubVectorHandler& WorkVec,
                                                                             doublereal dCoef,
                                                                             const VectorHandler& XCurr,
                                                                             const VectorHandler& XPrimeCurr)
{
     using namespace sp_grad;

     SolidElemDynamicType::AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);

     for (index_type i = 1; i <= iNumNodes; ++i) {
          this->rgNodes[i - 1]->AddInertia(diagM(i));
     }

     return WorkVec;
}

template <typename ElementType, typename CollocationType, typename SolidCSLType>
void
SolidElemDynamicDiagMass<ElementType, CollocationType, SolidCSLType>::AssDiagMassMatrix()
{
     using namespace sp_grad;

     SpColVectorA<doublereal, 3> r;
     SpColVectorA<doublereal, iNumNodes> h;
     SpMatrixA<doublereal, iNumNodes, 3> hd;

     for (index_type iColloc = 0; iColloc < iNumEvalPointsMassLumped; ++iColloc) {
          const doublereal alpha = CollocationType::dGetWeightMassLumped(iColloc);

          CollocationType::GetPositionMassLumped(iColloc, r);
          ElementType::ShapeFunction(r, h);
          ElementType::ShapeFunctionDeriv(r, hd);

          const SpMatrix<doublereal, 3, 3> J = Transpose(this->x0 * hd);
          const doublereal rho = Dot(h, this->rhon); // interpolate from nodes to collocation points
          const doublereal dm = rho * Det(J) * alpha;

          for (index_type k = 1; k <= iNumNodes; ++k) {
               diagM(k) += dm * std::pow(h(k), 2);
          }
     }
}

template <typename ElementType, typename CollocationType>
SolidElem*
ReadSolid(DataManager* const pDM, MBDynParser& HP, const unsigned int uLabel)
{
     DEBUGCOUTFNAME("ReadSolid");

     using namespace sp_grad;

     typedef SolidElemStatic<ElementType, CollocationType, SolidElasticConstLaw6D>      SolidStaticElemTypeElastic;
     typedef SolidElemStatic<ElementType, CollocationType, SolidViscoElasticConstLaw6D> SolidStaticElemTypeViscoElastic;
     typedef SolidElemStatic<ElementType, CollocationType, IsotropicLinearElastic>      SolidStaticElemTypeElasticIsotropic;
     typedef SolidElemStatic<ElementType, CollocationType, IsotropicLinearViscoElastic> SolidStaticElemTypeViscoElasticIsotropic;

     typedef SolidElemDynamic<ElementType, CollocationType, SolidElasticConstLaw6D>      SolidDynamicElemTypeElasticNoDiagMass;
     typedef SolidElemDynamic<ElementType, CollocationType, SolidViscoElasticConstLaw6D> SolidDynamicElemTypeViscoElasticNoDiagMass;
     typedef SolidElemDynamic<ElementType, CollocationType, IsotropicLinearElastic>      SolidDynamicElemTypeElasticIsotropicNoDiagMass;
     typedef SolidElemDynamic<ElementType, CollocationType, IsotropicLinearViscoElastic> SolidDynamicElemTypeViscoElasticIsotropicNoDiagMass;

     typedef SolidElemDynamicDiagMass<ElementType, CollocationType, SolidElasticConstLaw6D>      SolidDynamicElemTypeElasticDiagMass;
     typedef SolidElemDynamicDiagMass<ElementType, CollocationType, SolidViscoElasticConstLaw6D> SolidDynamicElemTypeViscoElasticDiagMass;
     typedef SolidElemDynamicDiagMass<ElementType, CollocationType, IsotropicLinearElastic>      SolidDynamicElemTypeElasticIsotropicDiagMass;
     typedef SolidElemDynamicDiagMass<ElementType, CollocationType, IsotropicLinearViscoElastic> SolidDynamicElemTypeViscoElasticIsotropicDiagMass;

     typedef typename
          std::conditional<ElementType::bHaveDiagMass,
                           SolidDynamicElemTypeElasticNoDiagMass,
                           SolidDynamicElemTypeElasticDiagMass>::type
          SolidDynamicElemTypeElastic;

     typedef typename
          std::conditional<ElementType::bHaveDiagMass,
                           SolidDynamicElemTypeViscoElasticNoDiagMass,
                           SolidDynamicElemTypeViscoElasticDiagMass>::type
          SolidDynamicElemTypeViscoElastic;

     typedef typename
          std::conditional<ElementType::bHaveDiagMass,
                           SolidDynamicElemTypeElasticIsotropicNoDiagMass,
                           SolidDynamicElemTypeElasticIsotropicDiagMass>::type
          SolidDynamicElemTypeElasticIsotropic;

     typedef typename
          std::conditional<ElementType::bHaveDiagMass,
                           SolidDynamicElemTypeViscoElasticIsotropicNoDiagMass,
                           SolidDynamicElemTypeViscoElasticIsotropicDiagMass>::type
          SolidDynamicElemTypeViscoElasticIsotropic;

     constexpr index_type iNumNodes = ElementType::iNumNodes;
     constexpr index_type iNumEvalPointsStiffness = CollocationType::iNumEvalPointsStiffness;

     std::array<const StructDispNodeAd*, iNumNodes> rgNodes;
     SpColVectorA<doublereal, iNumNodes> rhon;
     const bool bStaticModel = pDM->bIsStaticModel();

     for (index_type i = 0; i < iNumNodes; ++i) {
          rgNodes[i] = bStaticModel
               ? pDM->ReadNode<const StructDispNodeAd, Node::STRUCTURAL>(HP)
               : pDM->ReadNode<const DynamicStructDispNodeAd, Node::STRUCTURAL>(HP);
     }

     for (index_type i = 1; i <= iNumNodes; ++i) {
          rhon(i) = HP.GetReal();

          if (rhon(i) <= 0.) {
               silent_cerr("solid(" << uLabel << ") density must be greater than zero at line "
                           << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }
     }

     std::array<SolidMaterialData, iNumEvalPointsStiffness> rgMaterialData;

     ConstLawType::Type eConstLawType = ConstLawType::UNKNOWN;
     bool bIsotropicMaterial = false;

     if (HP.IsKeyWord("linear" "elastic" "isotropic")) {
          eConstLawType = ConstLawType::ELASTIC;
          bIsotropicMaterial = true;
     } else if (HP.IsKeyWord("linear" "viscoelastic" "isotropic")) {
          eConstLawType = ConstLawType::VISCOELASTIC;
          bIsotropicMaterial = true;
     }

     for (index_type i = 0; i < iNumEvalPointsStiffness; ++i) {
          if (bIsotropicMaterial) {
               if (i > 0 && HP.IsKeyWord("same")) {
                    rgMaterialData[i].E = rgMaterialData[i - 1].E;
                    rgMaterialData[i].nu = rgMaterialData[i - 1].nu;
                    rgMaterialData[i].beta = rgMaterialData[i - 1].beta;
                    continue;
               }

               rgMaterialData[i].E = HP.GetReal();

               if (rgMaterialData[i].E <= 0.) {
                    silent_cerr("solid(" << uLabel
                                << "): elastic modulus must be greater than zero at line "
                                << HP.GetLineData() << "\n");
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }

               rgMaterialData[i].nu = HP.GetReal();

               if (rgMaterialData[i].nu < 0. || rgMaterialData[i].nu >= 0.5) {
                    // FIXME: incompressible case not implemented yet
                    silent_cerr("solid(" << uLabel
                                << "): Poisson ratio must be between zero and 0.5 at line "
                                << HP.GetLineData() << "\n");
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }

               if (eConstLawType == ConstLawType::VISCOELASTIC) {
                    rgMaterialData[i].beta = HP.GetReal();

                    if (rgMaterialData[i].beta < 0.) {
                         silent_cerr("solid(" << uLabel
                                     << "): damping coefficient must be greater than or equal to zero at line "
                                     << HP.GetLineData() << "\n");
                         throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
               }
          } else {
               if (i > 0 && HP.IsKeyWord("same")) {
                    rgMaterialData[i].pCSL.reset(rgMaterialData[i - 1].pCSL->pCopy());
                    continue;
               }

               ConstLawType::Type CLType_I = ConstLawType::UNKNOWN;

               rgMaterialData[i].pCSL.reset(HP.GetConstLaw6D(CLType_I));

               if (rgMaterialData[i].pCSL->iGetNumDof() != 0) {
                    silent_cerr("line " << HP.GetLineData()
                                << ": solid(" << uLabel << ") "
                                "does not support dynamic constitutive laws yet"
                                << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }

               if (eConstLawType == ConstLawType::UNKNOWN) {
                    switch (CLType_I) {
                    case ConstLawType::ELASTIC:
                    case ConstLawType::VISCOELASTIC:
                         break;
                    default:
                         silent_cerr("solid(" << uLabel << ") does not support this constitutive law type "
                                     "at line " << HP.GetLineData() << "\n");
                         throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                    eConstLawType = CLType_I;
               } else if (eConstLawType != CLType_I) {
                    silent_cerr("solid(" << uLabel << ") all constitutive laws must have the same type "
                                "at line " << HP.GetLineData() << "\n");
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
               }
          }
     }

     if (bIsotropicMaterial && eConstLawType == ConstLawType::VISCOELASTIC) {
          eConstLawType == ConstLawType::ELASTIC;

          for (index_type i = 0; i < iNumEvalPointsStiffness; ++i) {
               if (rgMaterialData[i].beta != 0.) {
                    eConstLawType = ConstLawType::VISCOELASTIC;
                    break;
               }
          }
     }

     const flag fOut = pDM->fReadOutput(HP, Elem::SOLID);

     std::ostream& out = pDM->GetLogFile();

     out << ElementType::ElementName() << ": " << uLabel;

     for (sp_grad::index_type i = 0; i < iNumNodes; ++i) {
          out << ' ' << rgNodes[i]->GetLabel();
     }

     out << '\n';

     SolidElem* pEl = nullptr;

     if (bStaticModel) {
          switch (eConstLawType) {
          case ConstLawType::ELASTIC:
               if (bIsotropicMaterial) {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidStaticElemTypeElasticIsotropic,
                                           SolidStaticElemTypeElasticIsotropic(uLabel,
                                                                               rgNodes,
                                                                               rhon,
                                                                               std::move(rgMaterialData),
                                                                               fOut));
               } else {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidStaticElemTypeElastic,
                                           SolidStaticElemTypeElastic(uLabel,
                                                                      rgNodes,
                                                                      rhon,
                                                                      std::move(rgMaterialData),
                                                                      fOut));
               }
               break;
          case ConstLawType::VISCOELASTIC:
               if (bIsotropicMaterial) {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidStaticElemTypeViscoElasticIsotropic,
                                           SolidStaticElemTypeViscoElasticIsotropic(uLabel,
                                                                                    rgNodes,
                                                                                    rhon,
                                                                                    std::move(rgMaterialData),
                                                                                    fOut));
               } else {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidStaticElemTypeViscoElastic,
                                           SolidStaticElemTypeViscoElastic(uLabel,
                                                                           rgNodes,
                                                                           rhon,
                                                                           std::move(rgMaterialData),
                                                                           fOut));
               }
               break;
          default:
               ASSERT(0);
          }
     } else {
          switch (eConstLawType) {
          case ConstLawType::ELASTIC:
               if (bIsotropicMaterial) {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidDynamicElemTypeElasticIsotropic,
                                           SolidDynamicElemTypeElasticIsotropic(uLabel,
                                                                                rgNodes,
                                                                                rhon,
                                                                                std::move(rgMaterialData),
                                                                                pDM->pGetRBK(),
                                                                                fOut));
               } else {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidDynamicElemTypeElastic,
                                           SolidDynamicElemTypeElastic(uLabel,
                                                                       rgNodes,
                                                                       rhon,
                                                                       std::move(rgMaterialData),
                                                                       pDM->pGetRBK(),
                                                                       fOut));
               }
               break;
          case ConstLawType::VISCOELASTIC:
               if (bIsotropicMaterial) {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidDynamicElemTypeViscoElasticIsotropic,
                                           SolidDynamicElemTypeViscoElasticIsotropic(uLabel,
                                                                                     rgNodes,
                                                                                     rhon,
                                                                                     std::move(rgMaterialData),
                                                                                     pDM->pGetRBK(),
                                                                                     fOut));
               } else {
                    SAFENEWWITHCONSTRUCTOR(pEl,
                                           SolidDynamicElemTypeViscoElastic,
                                           SolidDynamicElemTypeViscoElastic(uLabel,
                                                                            rgNodes,
                                                                            rhon,
                                                                            std::move(rgMaterialData),
                                                                            pDM->pGetRBK(),
                                                                            fOut));
               }
               break;
          default:
               ASSERT(0);
          }
     }

     if (HP.IsArg()) {
          silent_cerr("semicolon expected "
                      "at line " << HP.GetLineData() << std::endl);
          throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
     }

     return pEl;
}

template SolidElem* ReadSolid<Hexahedron8, Gauss2>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Hexahedron20, Gauss3>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Hexahedron20r, GaussH20r>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Pentahedron15, GaussP15>(DataManager*, MBDynParser&, unsigned int);
template SolidElem* ReadSolid<Tetrahedron10h, GaussT10h>(DataManager*, MBDynParser&, unsigned int);

void InitSolidCL()
{
     SetCL6D("neo" "hookean", new NeoHookeanRead);
}
