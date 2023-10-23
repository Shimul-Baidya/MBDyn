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

/* Legami costitutivi */

#ifndef CONSTLTP_H
#define CONSTLTP_H

#include "withlab.h"
#include "simentity.h"
#include "tpldrive.h"

#include "matvec3.h"
#include "matvec6.h"
#include "sp_gradient.h"
#include "sp_matrix_base.h"

typedef sp_grad::SpColVector<doublereal, 7> Vec7;
typedef sp_grad::SpMatrix<doublereal, 7, 7> Mat7x7;

template <typename T>
struct ConstLawHelper {
     static constexpr sp_grad::index_type iDim = T::iNumRowsStatic;
};

template <>
struct ConstLawHelper<doublereal> {
     static constexpr sp_grad::index_type iDim = 1;
};

/* Tipi di cerniere deformabili */
class ConstLawType {
public:
        enum Type {
                UNKNOWN        = 0x0,

                ELASTIC        = 0x1,
                VISCOUS        = 0x2,
                INCOMPRESSIBLE = 0x4,
                VISCOELASTIC = (ELASTIC | VISCOUS),
                ELASTICINCOMPR = (ELASTIC | INCOMPRESSIBLE)
        };
};

/* ConstitutiveLaw - begin */
template <class T, class Tder>
class ConstitutiveLawBase : public WithLabel, public SimulationEntity {
public:
	using SimulationEntity::AfterConvergence;
	using SimulationEntity::Update;

	class ErrNotAvailable : public MBDynErrBase {
	public:
		ErrNotAvailable(MBDYN_EXCEPT_ARGS_DECL)
		: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU)
		{
			silent_cerr("Constitutive law not available "
				"for this dimensionality"
				<< std::endl);
		};
		ErrNotAvailable(std::ostream& out, MBDYN_EXCEPT_ARGS_DECL)
		: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU)
		{
			out << "Constitutive law not available "
				"for this dimensionality"
				<< what()
				<< std::endl;
		};
	};

        typedef typename ConstitutiveLawBase<T, Tder>::ErrNotAvailable Err;   

        typedef T StrainType;
        typedef T StressType;
        typedef Tder StressDerStrainType;
protected:
        T Epsilon;		/* strain */
        T EpsilonPrime;		/* strain rate */

        T F;			/* force */
        Tder FDE;		/* force strain derivative */
        Tder FDEPrime;		/* force strain rate derivative */

public:
        static constexpr sp_grad::index_type iDim = ConstLawHelper<T>::iDim;

        ConstitutiveLawBase(void)
        : WithLabel(0),
        Epsilon(mb_zero<T>()), EpsilonPrime(mb_zero<T>()),
        F(mb_zero<T>()), FDE(mb_zero<Tder>()), FDEPrime(mb_zero<Tder>()) {
                NO_OP;
        };

        virtual ~ConstitutiveLawBase(void) {
                NO_OP;
        };

        virtual ConstLawType::Type GetConstLawType(void) const = 0;

        virtual std::ostream& Restart(std::ostream& out) const {
                return out;
        };

        // Main interface used by all the elements without automatic differentiation
        virtual void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) = 0;
     
        virtual void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
                NO_OP;
        };

        virtual const T& GetEpsilon(void) const {
                return Epsilon;
        };

        virtual const T& GetEpsilonPrime(void) const {
                return EpsilonPrime;
        };

        virtual const T& GetF(void) const {
                return F;
        };

        virtual const Tder& GetFDE(void) const {
                return FDE;
        };

        virtual const Tder& GetFDEPrime(void) const {
                return FDEPrime;
        };

        /* simentity */
        virtual unsigned int iGetNumDof(void) const {
                return 0;
        };

        virtual std::ostream& DescribeDof(std::ostream& out,
                        const char *prefix = "",
                        bool bInitial = false) const
        {
                return out;
        };
        virtual void DescribeDof(std::vector<std::string>& desc,
                        bool bInitial = false, int i = -1) const
        {
                ASSERT(i <= 0);
                desc.resize(0);
        };

        virtual std::ostream& DescribeEq(std::ostream& out,
                        const char *prefix = "",
                        bool bInitial = false) const
        {
                return out;
        };

        virtual void DescribeEq(std::vector<std::string>& desc,
                        bool bInitial = false, int i = -1) const
        {
                ASSERT(i <= 0);
                desc.resize(0);
        };

        virtual DofOrder::Order GetDofType(unsigned int i) const {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        };
};

template <typename T, typename Tder>
class ConstitutiveLawAd: public ConstitutiveLawBase<T, Tder> {
public:
     using ConstitutiveLawBase<T, Tder>::iDim;
     using ConstitutiveLawBase<T, Tder>::Update;
     using ConstitutiveLawBase<T, Tder>::F;
     using ConstitutiveLawBase<T, Tder>::FDE;
     using ConstitutiveLawBase<T, Tder>::FDEPrime;
     
     // The following functions provide a default implementation of the interfaces,
     // used by elements with automatic differentiation, for all the constitutive laws
     // without explicit support for automatic differentiation.
     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
            sp_grad::SpColVector<doublereal, iDim>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
            sp_grad::SpColVector<sp_grad::SpGradient, iDim>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
            sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
            const sp_grad::SpColVector<doublereal, iDim>& EpsPrime,
            sp_grad::SpColVector<doublereal, iDim>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
            const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& EpsPrime,
            sp_grad::SpColVector<sp_grad::SpGradient, iDim>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap);

     virtual void
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
            const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& EpsPrime,
            sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& FTmp,
            const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap);

protected:
     // Default implementation of the constitutive law matrix for constitutive laws with
     // automatic differentiation, used by elements without automatic differentiation.
     // Those functions may be instantiated by derived classes.
     template <typename ConstLaw>
     static inline void
     UpdateElasticSparse(ConstLaw* pCl, const T& Eps);

     template <typename ConstLaw>
     static inline void
     UpdateViscoelasticSparse(ConstLaw* pCl, const T& Eps, const T& EpsPrime);
};

template <>
class ConstitutiveLawAd<doublereal, doublereal>: public ConstitutiveLawBase<doublereal, doublereal> {    
     // TODO: add dedicated Update functions for scalar types
};

template <typename T, typename Tder>
class ConstitutiveLaw: public ConstitutiveLawAd<T, Tder> {
     typedef ConstitutiveLawBase<T, Tder> ConstLawBaseType;
     typedef ConstitutiveLawAd<T, Tder> ConstLawAdType;
public:
     virtual ConstitutiveLaw* pCopy() const = 0;

     using ConstLawBaseType::iDim;
     using ConstLawBaseType::Update;
     using ConstLawAdType::Update;
     using ConstLawBaseType::GetConstLawType;
     using ConstLawBaseType::Restart;
     using ConstLawBaseType::AfterConvergence;
     using ConstLawBaseType::GetEpsilon;
     using ConstLawBaseType::GetEpsilonPrime;
     using ConstLawBaseType::GetF;
     using ConstLawBaseType::GetFDE;
     using ConstLawBaseType::GetFDEPrime;
     using ConstLawBaseType::iGetNumDof;
     using ConstLawBaseType::DescribeDof;
     using ConstLawBaseType::DescribeEq;
     using ConstLawBaseType::GetDofType;
     using ConstLawBaseType::F;
     using ConstLawBaseType::FDE;
     using ConstLawBaseType::FDEPrime;
};

typedef ConstitutiveLaw<doublereal, doublereal> ConstitutiveLaw1D;
typedef ConstitutiveLaw<Vec3, Mat3x3> ConstitutiveLaw3D;
typedef ConstitutiveLaw<Vec6, Mat6x6> ConstitutiveLaw6D;
typedef ConstitutiveLaw<Vec7, Mat7x7> ConstitutiveLaw7D;

template <typename T, typename Tder>
void
ConstitutiveLawAd<T, Tder>::Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
                                   sp_grad::SpColVector<doublereal, iDim>& FTmp,
                                   const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     // FIXME: Need to copy the data because vectors are not compatible
     this->Epsilon = Eps;

     this->Update(this->Epsilon, this->EpsilonPrime);

     // FIXME: Return by output argument because move constructors will not be most efficient with static vector size
     FTmp = this->GetF();
}

template <typename T, typename Tder>
void
ConstitutiveLawAd<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
                                   sp_grad::SpColVector<sp_grad::SpGradient, iDim>& FTmp,
                                   const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).ResizeReset(this->F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());

          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).AddDeriv(FTmp(i), this->FDE(i, j), oDofMap.GetDofMap());
          }
     }
}

template <typename T, typename Tder>
void
ConstitutiveLawAd<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
                                   sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& FTmp,
                                   const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) == 0);

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).Reset(this->F(i));

          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).InsertDeriv(FTmp(i), this->FDE(i, j));
          }
     }
}

template <typename T, typename Tder>
void
ConstitutiveLawAd<T, Tder>::Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
                                   const sp_grad::SpColVector<doublereal, iDim>& EpsPrime,
                                   sp_grad::SpColVector<doublereal, iDim>& FTmp,
                                   const sp_grad::SpGradExpDofMapHelper<doublereal>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);

     this->Epsilon = Eps; // FIXME: Need to copy the data because vectors are not compatible
     this->EpsilonPrime = EpsPrime;

     this->Update(this->Epsilon, this->EpsilonPrime);

     FTmp = this->GetF();
}

template <typename T, typename Tder>
void
ConstitutiveLawAd<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
                                   const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& EpsPrime,
                                   sp_grad::SpColVector<sp_grad::SpGradient, iDim>& FTmp,
                                   const sp_grad::SpGradExpDofMapHelper<sp_grad::SpGradient>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).ResizeReset(this->F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap.GetDofMap());

          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).AddDeriv(FTmp(i), this->FDE(i, j), oDofMap.GetDofMap());
               EpsPrime(j).AddDeriv(FTmp(i), this->FDEPrime(i, j), oDofMap.GetDofMap());
          }
     }
}

template <typename T, typename Tder>
void
ConstitutiveLawAd<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
                                   const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& EpsPrime,
                                   sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& FTmp,
                                   const sp_grad::SpGradExpDofMapHelper<sp_grad::GpGradProd>& oDofMap)
{
     using namespace sp_grad;

     ASSERT(this->iGetNumDof() == 0);
     ASSERT((this->GetConstLawType() & ConstLawType::VISCOUS) != 0);

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).Reset(this->F(i));

          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).InsertDeriv(FTmp(i), this->FDE(i, j));
               EpsPrime(j).InsertDeriv(FTmp(i), this->FDEPrime(i, j));
          }
     }
}

template <class T, class Tder>
template <typename ConstLaw>
void
ConstitutiveLawAd<T, Tder>::UpdateViscoelasticSparse(ConstLaw* pCl, const T& Eps, const T& EpsPrime)
{
        using namespace sp_grad;
        constexpr index_type N = ConstLawHelper<T>::iDim;
        SpColVectorA<SpGradient, N> gEps, gEpsPrime, gF;

        pCl->Epsilon = Eps;
        pCl->EpsilonPrime = EpsPrime;

        for (index_type i = 1; i <= N; ++i)
        {
             gEps(i).Reset(Eps(i), i, 1.);
             gEpsPrime(i).Reset(EpsPrime(i), i + N, 1.);
        }

        SpGradExpDofMapHelper<SpGradient> oDofMap;

        oDofMap.GetDofStat(gEps);
        oDofMap.GetDofStat(gEpsPrime);
        oDofMap.Reset();
        oDofMap.InsertDof(gEps);
        oDofMap.InsertDof(gEpsPrime);
        oDofMap.InsertDone();

        pCl->UpdateViscoelasticTpl(gEps, gEpsPrime, gF, oDofMap);

        for (index_type j = 1; j <= N; ++j) {
             for (index_type i = 1; i <= N; ++i) {
                  pCl->FDE(i, j) = pCl->FDEPrime(i, j) = 0.;
             }
        }

        for (index_type i = 1; i <= N; ++i) {
                pCl->F(i) = gF(i).dGetValue();

                for (const auto& oDer: gF(i)) {
                     ASSERT(oDer.iDof >= 1);
                     ASSERT(oDer.iDof <= 2 * N);

                     if (oDer.iDof <= N) {
                          pCl->FDE(i, oDer.iDof) += oDer.dDer;
                     } else {
                          pCl->FDEPrime(i, oDer.iDof - N) += oDer.dDer;
                     }
                }
        }
}

template <class T, class Tder>
template <typename ConstLaw>
void
ConstitutiveLawAd<T, Tder>::UpdateElasticSparse(ConstLaw* pCl, const T& Eps)
{
        using namespace sp_grad;
        constexpr index_type N = ConstLawHelper<T>::iDim;
        SpColVectorA<SpGradient, N> gEps, gF;

        pCl->Epsilon = Eps;

        for (index_type i = 1; i <= N; ++i)
        {
             gEps(i).Reset(Eps(i), i, 1.);
        }

        SpGradExpDofMapHelper<SpGradient> oDofMap;

        oDofMap.GetDofStat(gEps);
        oDofMap.Reset();
        oDofMap.InsertDof(gEps);
        oDofMap.InsertDone();

        pCl->UpdateElasticTpl(gEps, gF, oDofMap);

        for (index_type j = 1; j <= N; ++j) {
             for (index_type i = 1; i <= N; ++i) {
                  pCl->FDE(i, j) = 0.;
             }
        }

        for (index_type i = 1; i <= N; ++i) {
                pCl->F(i) = gF(i).dGetValue();

                for (const auto& oDer: gF(i)) {
                     ASSERT(oDer.iDof >= 1);
                     ASSERT(oDer.iDof <= N);
                     pCl->FDE(i, oDer.iDof) += oDer.dDer;
                }
        }
}

/* ConstitutiveLaw - end */


/* ConstitutiveLawOwner - begin */

template <class T, class Tder>
class ConstitutiveLawOwner : public SimulationEntity {
protected:
        mutable ConstitutiveLaw<T, Tder>* pConstLaw;

public:
        ConstitutiveLawOwner(const ConstitutiveLaw<T, Tder>* pCL)
        : pConstLaw(const_cast<ConstitutiveLaw<T, Tder> *>(pCL)) {
                ASSERT(pCL != NULL);
        };

        virtual ~ConstitutiveLawOwner(void) {
                ASSERT(pConstLaw != NULL);
                if (pConstLaw != NULL) {
                        SAFEDELETE(pConstLaw);
                }
        };

        inline ConstitutiveLaw<T, Tder>* pGetConstLaw(void) const {
                ASSERT(pConstLaw != NULL);
                return pConstLaw;
        };

	using SimulationEntity::Update;
        inline void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
                ASSERT(pConstLaw != NULL);
                pConstLaw->Update(Eps, EpsPrime);
        };

	using SimulationEntity::AfterConvergence;
        inline void AfterConvergence(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
                ASSERT(pConstLaw != NULL);
                pConstLaw->AfterConvergence(Eps, EpsPrime);
        };

        inline const T& GetF(void) const {
                ASSERT(pConstLaw != NULL);
                return pConstLaw->GetF();
        };

        inline const Tder& GetFDE(void) const {
                ASSERT(pConstLaw != NULL);
                return pConstLaw->GetFDE();
        };

        inline const Tder& GetFDEPrime(void) const {
                ASSERT(pConstLaw != NULL);
                return pConstLaw->GetFDEPrime();
        };

        /* simentity */
        virtual unsigned int iGetNumDof(void) const override {
                ASSERT(pConstLaw != NULL);
                return pConstLaw->iGetNumDof();
        };

        virtual std::ostream& DescribeDof(std::ostream& out,
                        const char *prefix = "",
                        bool bInitial = false) const override
        {
                return out;
        };

        virtual void DescribeDof(std::vector<std::string>& desc,
                        bool bInitial = false, int i = -1) const override
        {
                ASSERT(i <= 0);
                desc.resize(0);
        };

        virtual std::ostream& DescribeEq(std::ostream& out,
                        const char *prefix = "",
                        bool bInitial = false) const override
        {
                return out;
        };

        virtual void DescribeEq(std::vector<std::string>& desc,
                        bool bInitial = false, int i = -1) const override
        {
                ASSERT(i <= 0);
                desc.resize(0);
        };

        virtual DofOrder::Order GetDofType(unsigned int i) const override {
                ASSERT(pConstLaw != NULL);
                return pConstLaw->GetDofType(i);
        };

        /*
         * Metodi per l'estrazione di dati "privati".
         * Si suppone che l'estrattore li sappia interpretare.
         * Come default non ci sono dati privati estraibili
         */
        virtual unsigned int iGetNumPrivData(void) const override {
                return pConstLaw->iGetNumPrivData();
        };

        /*
         * Maps a string (possibly with substrings) to a private data;
         * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0
         * in case of unrecognized data; error must be handled by caller
         */
        virtual unsigned int iGetPrivDataIdx(const char *s) const override {
                return pConstLaw->iGetPrivDataIdx(s);
        };

        /*
         * Returns the current value of a private data
         * with 0 < i <= iGetNumPrivData()
         */
        virtual doublereal dGetPrivData(unsigned int i) const override {
                return pConstLaw->dGetPrivData(i);
        };

        virtual std::ostream& OutputAppend(std::ostream& out) const override {
                return pConstLaw->OutputAppend(out);
        };

        virtual void NetCDFOutputAppend(OutputHandler& OH) const override {
                return pConstLaw->NetCDFOutputAppend(OH);
        };

        virtual void OutputAppendPrepare(OutputHandler& OH, const std::string& name) override {
                pConstLaw->OutputAppendPrepare(OH, name);
        };

};

typedef ConstitutiveLawOwner<doublereal, doublereal> ConstitutiveLaw1DOwner;
typedef ConstitutiveLawOwner<Vec3, Mat3x3> ConstitutiveLaw3DOwner;
typedef ConstitutiveLawOwner<Vec6, Mat6x6> ConstitutiveLaw6DOwner;

/* ConstitutiveLawOwner - end */

/* functions that read a constitutive law */
extern ConstitutiveLaw<doublereal, doublereal> *
ReadCL1D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
extern ConstitutiveLaw<Vec3, Mat3x3> *
ReadCL3D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
extern ConstitutiveLaw<Vec6, Mat6x6> *
ReadCL6D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);
extern ConstitutiveLaw<Vec7, Mat7x7> *
ReadCL7D(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType);

/* prototype of the template functional object: reads a constitutive law */
template <class T, class Tder>
struct ConstitutiveLawRead {
        virtual ~ConstitutiveLawRead<T, Tder>( void ) { NO_OP; };
        virtual ConstitutiveLaw<T, Tder> *
        Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) = 0;
};

/* constitutive law registration functions: call to register one */
extern bool
SetCL1D(const char *name, ConstitutiveLawRead<doublereal, doublereal> *rf);
extern bool
SetCL3D(const char *name, ConstitutiveLawRead<Vec3, Mat3x3> *rf);
extern bool
SetCL6D(const char *name, ConstitutiveLawRead<Vec6, Mat6x6> *rf);
extern bool
SetCL7D(const char *name, ConstitutiveLawRead<Vec7, Mat7x7> *rf);
/* create/destroy */
extern void InitCL(void);
extern void DestroyCL(void);

#endif /* CONSTLTP_H */
