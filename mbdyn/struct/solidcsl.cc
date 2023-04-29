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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "constltp.h"
#include "dataman.h"
#include "solidcsl.h"

class NeoHookean: public ConstitutiveLaw<Vec6, Mat6x6> {
protected:
     NeoHookean(const doublereal mu, const doublereal lambda)
          :mu(mu), lambda(lambda) {
     }

     const doublereal mu, lambda;
};

class NeoHookeanElastic: public NeoHookean {
public:
     NeoHookeanElastic(const doublereal mu, const doublereal lambda)
          :NeoHookean(mu, lambda) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
          return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy() const override {
          NeoHookeanElastic* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 NeoHookeanElastic,
                                 NeoHookeanElastic(mu, lambda));
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

          static constexpr index_type i1[] = {1, 2, 3, 1, 2, 3};
          static constexpr index_type i2[] = {1, 2, 3, 2, 3, 1};

          for (index_type i = 1; i <= 6; ++i) {
               const index_type j = i1[i - 1];
               const index_type k = i2[i - 1];
               const bool deltajk = (j == k);

               oDofMap.MapAssign(sigma(i), mu * deltajk + (CC(j, k) - C(j, k) * IC + IIC * deltajk) * gamma);
          }
     }
};

class NeoHookeanViscoelastic: public NeoHookeanElastic {
public:
     NeoHookeanViscoelastic(const doublereal mu, const doublereal lambda, const doublereal beta)
          :NeoHookeanElastic(mu, lambda), beta(beta) {
          ASSERT(beta > 0.);
     }

     virtual ConstLawType::Type GetConstLawType() const override {
          return ConstLawType::VISCOELASTIC;
     }

     virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy() const override {
          NeoHookeanViscoelastic* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 NeoHookeanViscoelastic,
                                 NeoHookeanViscoelastic(mu, lambda, beta));
          return pCL;
     }

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
            const sp_grad::SpColVector<doublereal, iDim>& EpsPrime,
            sp_grad::SpColVector<doublereal, iDim>& FTmp) override {
          UpdateViscoelasticTpl(Eps, EpsPrime, FTmp);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
            const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& EpsPrime,
            sp_grad::SpColVector<sp_grad::SpGradient, iDim>& FTmp) override {
          UpdateViscoelasticTpl(Eps, EpsPrime, FTmp);
     }

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
            const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& EpsPrime,
            sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& FTmp) override {
          UpdateViscoelasticTpl(Eps, EpsPrime, FTmp);
     }

     virtual void
     Update(const Vec6& Eps, const Vec6& EpsPrime) override {
          ConstitutiveLaw<Vec6, Mat6x6>::UpdateViscoelasticSparse(this, Eps, EpsPrime);
     }

     template <typename VectorType>
     void UpdateViscoelasticTpl(const VectorType& epsilon, const VectorType& epsilonP, VectorType& sigma) {

          typedef typename VectorType::ValueType T;
          using namespace sp_grad;

          SpGradExpDofMapHelper<typename VectorType::ValueType> oDofMap;

          oDofMap.GetDofStat(epsilon);
          oDofMap.GetDofStat(epsilonP);
          oDofMap.Reset();
          oDofMap.InsertDof(epsilon);
          oDofMap.InsertDof(epsilonP);
          oDofMap.InsertDone();

          // Based on Lars Kuebler 2005, chapter 2.2.1.3, page 25-26
          // In contradiction to the reference, the effect of damping
          // is assumed proportional to initial stiffness

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

          T traceGP;

          oDofMap.MapAssign(traceGP, epsilonP(1) + epsilonP(2) + epsilonP(3));

          static constexpr index_type i1[] = {1, 2, 3, 1, 2, 3};
          static constexpr index_type i2[] = {1, 2, 3, 2, 3, 1};

          for (index_type i = 1; i <= 6; ++i) {
               const index_type j = i1[i - 1];
               const index_type k = i2[i - 1];
               const bool deltajk = (j == k);

               oDofMap.MapAssign(sigma(i), mu * deltajk + (CC(j, k) - C(j, k) * IC + IIC * deltajk) * gamma
                                 + (deltajk ? 2. : 1.) * mu * beta * epsilonP(i)
                                 + deltajk * beta * lambda * traceGP);
          }
     }
private:
     const doublereal beta;
};

class MooneyRivlinElastic: public ConstitutiveLaw<Vec6, Mat6x6> {
public:
     MooneyRivlinElastic(const doublereal C1, const doublereal C2, const doublereal kappa)
          :C1(C1), C2(C2), kappa(kappa) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
          return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy() const override {
          MooneyRivlinElastic* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 MooneyRivlinElastic,
                                 MooneyRivlinElastic(C1, C2, kappa));
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

          // Based on K.J. Bathe
          // Finite Element Procedures
          // 2nd Edition
          // Chapter 6.6.2, pages 592-594

          const SpMatrix<T, 3, 3> C{T{2. * epsilon(1) + 1.},        epsilon(4),           epsilon(6),
                                    epsilon(4), T{2. * epsilon(2) + 1.},          epsilon(5),
                                    epsilon(6),          epsilon(5), T{2. * epsilon(3) + 1.}};

          const SpMatrix<T, 3, 3> CC(C * C, oDofMap);

          T IC, IIC, IIIC;

          oDofMap.MapAssign(IC, C(1, 1) + C(2, 2) + C(3, 3));
          oDofMap.MapAssign(IIC, 0.5 * (IC * IC - (CC(1, 1) + CC(2, 2) + CC(3, 3))));

          Det(C, IIIC, oDofMap);

          T a1, a2, a3, a4;

          oDofMap.MapAssign(a1, 2. * C1 * pow(IIIC, -1./3.));
          oDofMap.MapAssign(a2, 2. * C2 * pow(IIIC, -2./3.));
          oDofMap.MapAssign(a3, a1 + a2 * IC);
          oDofMap.MapAssign(a4, (kappa * (IIIC - sqrt(IIIC)) - 2./3. * C1 * IC * pow(IIIC, -1./3.) - 4./3. * C2 * IIC * pow(IIIC, -2./3.)) / IIIC);

          static constexpr index_type i1[] = {1, 2, 3, 1, 2, 3};
          static constexpr index_type i2[] = {1, 2, 3, 2, 3, 1};

          for (index_type i = 1; i <= 6; ++i) {
               const index_type j = i1[i - 1];
               const index_type k = i2[i - 1];
               const bool deltajk = (j == k);

               oDofMap.MapAssign(sigma(i), a3 * deltajk - a2 * C(j, k) + a4 * (CC(j, k) - C(j, k) * IC + IIC * deltajk));
          }
     }
private:
     const doublereal C1, C2, kappa;
};

struct NeoHookeanRead: ConstitutiveLawRead<Vec6, Mat6x6> {
     static void
     ReadLameParameters(const DataManager* pDM, MBDynParser& HP, doublereal& mu, doublereal& lambda) {
          if (!HP.IsKeyWord("E")) {
               silent_cerr("keyword \"E\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal E = HP.GetReal(); // Young's modulus

          if (E <= 0.) {
               silent_cerr("E must be greater than zero at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          if (!HP.IsKeyWord("nu")) {
               silent_cerr("keyword \"nu\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal nu = HP.GetReal(); // Poisson's ratio

          if (nu < 0. || nu >= 0.5) {
               // FIXME: incompressible case not implemented yet
               silent_cerr("nu must be between 0 and 0.5 at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          mu = (E / (2. * (1. + nu))); // Lame parameters
          lambda = (E * nu / ((1. + nu) * (1. - 2. * nu)));
     }

     static NeoHookeanElastic*
     pCreateElastic(const doublereal mu, const doublereal lambda) {
          NeoHookeanElastic* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 NeoHookeanElastic,
                                 NeoHookeanElastic(mu, lambda));

          return pCL;
     }

     static NeoHookeanViscoelastic*
     pCreateViscoelastic(const doublereal mu, const doublereal lambda, const doublereal beta) {
          NeoHookeanViscoelastic* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 NeoHookeanViscoelastic,
                                 NeoHookeanViscoelastic(mu, lambda, beta));

          return pCL;
     }
};

struct NeoHookeanReadElastic: NeoHookeanRead {
     virtual ConstitutiveLaw<Vec6, Mat6x6>*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
          doublereal mu, lambda;

          ReadLameParameters(pDM, HP, mu, lambda);

          NeoHookeanElastic* pCL = pCreateElastic(mu, lambda);

          CLType = pCL->GetConstLawType();

          return pCL;
     }
};

struct NeoHookeanReadViscoelastic: NeoHookeanRead {
     virtual ConstitutiveLaw<Vec6, Mat6x6>*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
          doublereal mu, lambda;

          ReadLameParameters(pDM, HP, mu, lambda);

          if (!HP.IsKeyWord("beta")) {
               silent_cerr("keyword \"beta\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal beta = HP.GetReal(); // damping coefficient

          if (beta < 0.) {
               silent_cerr("beta must be greater than or equal to zero at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          NeoHookean* pCL = (beta == 0.)
               ? pCreateElastic(mu, lambda)
               : pCreateViscoelastic(mu, lambda, beta);


          CLType = pCL->GetConstLawType();

          return pCL;
     }
};

struct MooneyRivlinReadElastic: ConstitutiveLawRead<Vec6, Mat6x6> {
     virtual ConstitutiveLaw<Vec6, Mat6x6>*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
          if (!HP.IsKeyWord("C1")) {
               silent_cerr("keyword \"C1\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal C1 = HP.GetReal();

          if (C1 <= 0.) {
               silent_cerr("C1 must be greater than zero at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          if (!HP.IsKeyWord("C2")) {
               silent_cerr("keyword \"C2\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal C2 = HP.GetReal();

          if (C2 < 0.) {
               silent_cerr("C2 must be greater than or equal to zero at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          if (!HP.IsKeyWord("kappa")) {
               silent_cerr("keyword \"kappa\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal kappa = HP.GetReal();

          if (kappa <= 0. || kappa >= std::numeric_limits<doublereal>::max()) {
               silent_cerr("kappa must be between zero and realmax at line " << HP.GetLineData() << std::endl);
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          MooneyRivlinElastic* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 MooneyRivlinElastic,
                                 MooneyRivlinElastic(C1, C2, kappa));

          CLType = pCL->GetConstLawType();

          return pCL;
     }
};

class BilinearIsotropicHardening: public ConstitutiveLaw<Vec6, Mat6x6> {
public:
     BilinearIsotropicHardening(const doublereal E, const doublereal nu, const doublereal ET, const doublereal sigmayv)
          :E(E), nu(nu), ET(ET), sigmayv(sigmayv), sigmay_prev(sigmayv), sigmay_curr(sigmayv), aE((1. + nu) / E), EP(E * ET / (E - ET)) {
     }

     virtual ConstLawType::Type GetConstLawType() const override {
          return ConstLawType::ELASTIC;
     }

     virtual ConstitutiveLaw<Vec6, Mat6x6>* pCopy() const override {
          BilinearIsotropicHardening* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 BilinearIsotropicHardening,
                                 BilinearIsotropicHardening(E, nu, ET, sigmayv)); // Not considering eP_prev
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
          // Small strain elastoplasticity based on
          // K.-J. Bathe Finite Element Procedures second edition 2014 ISBN 978-0-9790049-5-7
          // Chapter 6.6.3, page 597

          typedef typename VectorType::ValueType T;
          using namespace sp_grad;

          SpGradExpDofMapHelper<typename VectorType::ValueType> oDofMap;

          oDofMap.GetDofStat(epsilon);
          oDofMap.Reset();
          oDofMap.InsertDof(epsilon);
          oDofMap.InsertDone();

          T em, sigmam; // mean stress, mean strain

          oDofMap.MapAssign(em, (epsilon(1) + epsilon(2) + epsilon(3)) / 3.); // equation 6.219
          oDofMap.MapAssign(sigmam, (E / (1 - 2 * nu)) * em); // equation 6.215

          static constexpr index_type idx1[] = {1, 2, 3, 1, 2, 3};
          static constexpr index_type idx2[] = {1, 2, 3, 2, 3, 1};
          static constexpr index_type idx_tens[3][3] = {{1, 4, 6},
                                                        {4, 2, 5},
                                                        {6, 5, 3}};

          SpMatrix<T, 3, 3> e2(3, 3, epsilon.iGetMaxSize());

          for (index_type i = 1; i <= 3; ++i) {
               for (index_type j = 1; j <= 3; ++j) {
                    oDofMap.MapAssign(e2(i, j), ((i == j) ? 1. : 0.5) * epsilon(idx_tens[i - 1][j - 1]) - em * (i == j) - eP_prev(i, j)); // equation 6.218, 6.221
               }
          }

          const SpMatrix<T, 3, 3> SE(e2 / aE, oDofMap); // elastic portion of deviatoric stress tensor equation 6.239

          T sum_e2ij_2;

          SpGradientTraits<T>::ResizeReset(sum_e2ij_2, 0., 2 * e2.iGetMaxSize() * e2.iGetNumRows() * e2.iGetNumCols());

          for (const auto& e2ij: e2) {
               sum_e2ij_2 += e2ij * e2ij;
          }

          T d{0.};

          if (sum_e2ij_2 > 0.) { // avoid division by zero
               oDofMap.MapAssign(d, sqrt((3. / 2.) * sum_e2ij_2)); // equation 6.232
          }

          T sigma_barE; // equivalent elastic stress solution

          oDofMap.MapAssign(sigma_barE, d / aE); // equation 6.238

          SpMatrix<T, 3, 3> eP(3, 3, epsilon.iGetMaxSize()); // plastic strain tensor
          SpMatrix<T, 3, 3> S(3, 3, epsilon.iGetMaxSize()); // deviatoric stress tensor

          T sigmay{sigmay_prev}; // equivalent yield stress

          if (sigma_barE > sigmay_prev) {
               // yielding
               T sigma_bar, lambda; // effective stress, yield parameter

               oDofMap.MapAssign(sigma_bar, (2. * EP * d + 3. * sigmay_prev) / (2. * EP * aE + 3.)); // equation 6.236
               oDofMap.MapAssign(lambda, d / sigma_bar - aE); // equation 6.230, 6.231

               const SpMatrix<T, 3, 3> Delta_eP(SE * (aE * lambda / (aE + lambda)), oDofMap); // equation 6.240, 6.241, 6.225

               sigmay = sigma_bar;
               eP.MapAssign(eP_prev + Delta_eP, oDofMap);
               S.MapAssign(SE - Delta_eP / aE, oDofMap); // stress correction equation 6.240
          } else {
               // not yielding
               eP = eP_prev;
               S = SE;
          }

          for (index_type i = 1; i <= 3; ++i) {
               oDofMap.MapAssign(sigma(i), S(i, i) + sigmam); // equation 6.216
          }

          for (index_type i = 4; i <= 6; ++i) {
               oDofMap.MapAssign(sigma(i), S(idx1[i - 1], idx2[i - 1])); // equation 6.216
          }

          UpdatePlasticStrain(eP, sigmay);
     }

     virtual void AfterConvergence(const Vec6& Eps, const Vec6& EpsPrime) override {
          eP_prev = eP_curr;
          sigmay_prev = sigmay_curr;
     }

private:
     void UpdatePlasticStrain(const sp_grad::SpMatrix<doublereal, 3, 3>& eP, doublereal sigmay) {
          eP_curr = eP;
          sigmay_curr = sigmay;
     }

     void UpdatePlasticStrain(const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>&,
                              const sp_grad::SpGradient&) {
     }

     void UpdatePlasticStrain(const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>&,
                              const sp_grad::GpGradProd&) {
     }

     const doublereal E, nu, ET, sigmayv;
     sp_grad::SpMatrixA<doublereal, 3, 3> eP_prev, eP_curr;
     doublereal sigmay_prev, sigmay_curr, aE, EP;
};

struct BilinearIsotropicHardeningRead: ConstitutiveLawRead<Vec6, Mat6x6> {
     virtual ConstitutiveLaw<Vec6, Mat6x6>*
     Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
          if (!HP.IsKeyWord("E")) {
               silent_cerr("keyword \"E\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }
          const doublereal E = HP.GetReal(); // Young's modulus

          if (E <= 0.) {
               silent_cerr("E must be greater than zero at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          if (!HP.IsKeyWord("nu")) {
               silent_cerr("keyword \"nu\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal nu = HP.GetReal(); // Poisson's ratio

          if (nu < 0. || nu >= 0.5) {
               // FIXME: incompressible case not implemented yet
               silent_cerr("nu must be between 0 and 0.5 at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          if (!HP.IsKeyWord("ET")) {
               silent_cerr("keyword \"ET\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal ET = HP.GetReal(); // Plastic tangent modulus

          if (ET < 0. || ET >= E) {
               silent_cerr("ET must be between zero and E at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          if (!HP.IsKeyWord("sigmayv")) {
               silent_cerr("keyword \"sigmayv\" expected at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          const doublereal sigmayv = HP.GetReal(); // Equivalent initial yield stress

          if (sigmayv <= 0.) {
               silent_cerr("sigmayv must be greater than zero at line " << HP.GetLineData() << "\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }

          BilinearIsotropicHardening* pCL = nullptr;

          SAFENEWWITHCONSTRUCTOR(pCL,
                                 BilinearIsotropicHardening,
                                 BilinearIsotropicHardening(E, nu, ET, sigmayv));

          CLType = ConstLawType::ELASTIC;

          return pCL;
     }
};

void InitSolidCSL()
{
     SetCL6D("neo" "hookean" "elastic", new NeoHookeanReadElastic);
     SetCL6D("neo" "hookean" "viscoelastic", new NeoHookeanReadViscoelastic);
     SetCL6D("mooney" "rivlin" "elastic", new MooneyRivlinReadElastic);
     SetCL6D("bilinear" "isotropic" "hardening", new BilinearIsotropicHardeningRead);
}
