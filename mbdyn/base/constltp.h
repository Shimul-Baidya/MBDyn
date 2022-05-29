/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

template <typename T>
struct ConstLawHelper;

template <>
struct ConstLawHelper<doublereal> {
     static constexpr sp_grad::index_type iDim = 1;
};

template <>
struct ConstLawHelper<Vec3> {
     static constexpr sp_grad::index_type iDim = Vec3::iNumRowsStatic;
};

template <>
struct ConstLawHelper<Vec6> {
     static constexpr sp_grad::index_type iDim = Vec6::iNumRowsStatic;
};

/* Tipi di cerniere deformabili */
class ConstLawType {
public:
	enum Type {
		UNKNOWN = 0x0,
		
		ELASTIC = 0x1,
		VISCOUS = 0x2,
		VISCOELASTIC = (ELASTIC | VISCOUS),
	
		LASTCONSTLAWTYPE
	};
};

/* ConstitutiveLaw - begin */
template <class T, class Tder>
class ConstitutiveLaw : public WithLabel, public SimulationEntity {
public:
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

	typedef typename ConstitutiveLaw<T, Tder>::ErrNotAvailable Err;   
   
protected:
	T Epsilon;		/* strain */
	T EpsilonPrime;		/* strain rate */

	T F;			/* force */
	Tder FDE;		/* force strain derivative */
	Tder FDEPrime;		/* force strain rate derivative */

public:
	ConstitutiveLaw(void)
	: WithLabel(0),
	Epsilon(mb_zero<T>()), EpsilonPrime(mb_zero<T>()),
	F(mb_zero<T>()), FDE(mb_zero<Tder>()), FDEPrime(mb_zero<Tder>()) {
		NO_OP;
	};

	virtual ~ConstitutiveLaw(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const = 0;

	virtual ConstitutiveLaw<T, Tder>* pCopy(void) const = 0;
	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};
	
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

     static constexpr sp_grad::index_type iDim = ConstLawHelper<T>::iDim;

     inline sp_grad::SpColVector<doublereal, iDim>
     Update(const sp_grad::SpColVector<doublereal, iDim>& Eps);

     inline sp_grad::SpColVector<sp_grad::SpGradient, iDim>
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps);

     inline sp_grad::SpColVector<sp_grad::GpGradProd, iDim>
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps);
     
     inline sp_grad::SpColVector<doublereal, iDim>
     Update(const sp_grad::SpColVector<doublereal, iDim>& Eps,
            const sp_grad::SpColVector<doublereal, iDim>& EpsPrime);

     inline sp_grad::SpColVector<sp_grad::SpGradient, iDim>
     Update(const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& Eps,
            const sp_grad::SpColVector<sp_grad::SpGradient, iDim>& EpsPrime);

     inline sp_grad::SpColVector<sp_grad::GpGradProd, iDim>
     Update(const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& Eps,
            const sp_grad::SpColVector<sp_grad::GpGradProd, iDim>& EpsPrime);
protected:
     template <typename ConstLaw>
     static inline void UpdateViscoelasticSparse(ConstLaw* pCl, const T& Eps, const T& EpsPrime);

     template <typename ConstLaw>
     static inline void UpdateElasticSparse(ConstLaw* pCl, const T& Eps);
};

typedef ConstitutiveLaw<doublereal, doublereal> ConstitutiveLaw1D;
typedef ConstitutiveLaw<Vec3, Mat3x3> ConstitutiveLaw3D;
typedef ConstitutiveLaw<Vec6, Mat6x6> ConstitutiveLaw6D;

template <typename T, typename Tder>
sp_grad::SpColVector<doublereal, ConstitutiveLaw<T, Tder>::iDim>
ConstitutiveLaw<T, Tder>::Update(const sp_grad::SpColVector<doublereal, ConstitutiveLaw<T, Tder>::iDim>& Eps)
{
     using namespace sp_grad;

     static_assert(iDim == T::iNumRowsStatic);
     static_assert(iDim >= 1);
     static_assert(T::iNumColsStatic == 1);

     ASSERT(iGetNumDof() == 0);
     ASSERT((GetConstLawType() & ConstLawType::VISCOUS) == 0);
     
     Update(T(Eps.begin()));

     return GetF();
}

template <typename T, typename Tder>
sp_grad::SpColVector<sp_grad::SpGradient, ConstitutiveLaw<T, Tder>::iDim>
ConstitutiveLaw<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::SpGradient, ConstitutiveLaw<T, Tder>::iDim>& Eps)
{
     using namespace sp_grad;

     static_assert(iDim == T::iNumRowsStatic);
     static_assert(iDim >= 1);
     static_assert(T::iNumColsStatic == 1);

     ASSERT(iGetNumDof() == 0);
     ASSERT((GetConstLawType() & ConstLawType::VISCOUS) == 0);

     SpGradDofStat oDofStat;

     for (const SpGradient& g: Eps) {
          g.GetDofStat(oDofStat);
     }
     
     SpGradExpDofMap oDofMap(oDofStat);

     for (const SpGradient& g: Eps) {
          g.InsertDof(oDofMap);
     }

     oDofMap.InsertDone();
     
     SpColVector<SpGradient, iDim> FTmp(iDim, oDofMap.iGetLocalSize());

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).ResizeReset(F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap);
          
          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).AddDeriv(FTmp(i), FDE(i, j), oDofMap);
          }
     }

     return FTmp;
}


template <typename T, typename Tder>
sp_grad::SpColVector<sp_grad::GpGradProd, ConstitutiveLaw<T, Tder>::iDim>
ConstitutiveLaw<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::GpGradProd, ConstitutiveLaw<T, Tder>::iDim>& Eps)
{
     using namespace sp_grad;

     static_assert(iDim == T::iNumRowsStatic);
     static_assert(iDim >= 1);
     static_assert(T::iNumColsStatic == 1);

     ASSERT(iGetNumDof() == 0);
     ASSERT((GetConstLawType() & ConstLawType::VISCOUS) == 0);

     SpColVector<GpGradProd, iDim> FTmp(iDim, 1);

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).Reset(F(i));
          
          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).InsertDeriv(FTmp(i), FDE(i, j));
          }
     }

     return FTmp;
}

template <typename T, typename Tder>
sp_grad::SpColVector<doublereal, ConstitutiveLaw<T, Tder>::iDim>
ConstitutiveLaw<T, Tder>::Update(const sp_grad::SpColVector<doublereal, ConstitutiveLaw<T, Tder>::iDim>& Eps,
                                 const sp_grad::SpColVector<doublereal, ConstitutiveLaw<T, Tder>::iDim>& EpsPrime)
{
     using namespace sp_grad;

     static_assert(iDim == T::iNumRowsStatic);
     static_assert(iDim >= 1);
     static_assert(T::iNumColsStatic == 1);

     ASSERT(iGetNumDof() == 0);
     ASSERT((GetConstLawType() & ConstLawType::VISCOUS) != 0);
     
     Update(T(Eps.begin()), T(EpsPrime.begin()));

     return GetF();
}

template <typename T, typename Tder>
sp_grad::SpColVector<sp_grad::SpGradient, ConstitutiveLaw<T, Tder>::iDim>
ConstitutiveLaw<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::SpGradient, ConstitutiveLaw<T, Tder>::iDim>& Eps,
                                 const sp_grad::SpColVector<sp_grad::SpGradient, ConstitutiveLaw<T, Tder>::iDim>& EpsPrime)
{
     using namespace sp_grad;

     static_assert(iDim == T::iNumRowsStatic);
     static_assert(iDim >= 1);
     static_assert(T::iNumColsStatic == 1);

     ASSERT(iGetNumDof() == 0);
     ASSERT((GetConstLawType() & ConstLawType::VISCOUS) != 0);
     
     SpGradDofStat oDofStat;

     for (const SpGradient& g: Eps) {
          g.GetDofStat(oDofStat);
     }

     for (const SpGradient& g: EpsPrime) {
          g.GetDofStat(oDofStat);
     }
     
     SpGradExpDofMap oDofMap(oDofStat);

     for (const SpGradient& g: Eps) {
          g.InsertDof(oDofMap);
     }

     for (const SpGradient& g: EpsPrime) {
          g.InsertDof(oDofMap);
     }

     oDofMap.InsertDone();
     
     SpColVector<SpGradient, iDim> FTmp(iDim, oDofMap.iGetLocalSize());

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).ResizeReset(F(i), oDofMap.iGetLocalSize());
          FTmp(i).InitDeriv(oDofMap);
          
          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).AddDeriv(FTmp(i), FDE(i, j), oDofMap);
               EpsPrime(j).AddDeriv(FTmp(i), FDEPrime(i, j), oDofMap);
          }
     }

     return FTmp;
}

template <typename T, typename Tder>
sp_grad::SpColVector<sp_grad::GpGradProd, ConstitutiveLaw<T, Tder>::iDim>
ConstitutiveLaw<T, Tder>::Update(const sp_grad::SpColVector<sp_grad::GpGradProd, ConstitutiveLaw<T, Tder>::iDim>& Eps,
                                 const sp_grad::SpColVector<sp_grad::GpGradProd, ConstitutiveLaw<T, Tder>::iDim>& EpsPrime)
{
     using namespace sp_grad;

     static_assert(iDim == T::iNumRowsStatic);
     static_assert(iDim >= 1);
     static_assert(T::iNumColsStatic == 1);

     ASSERT(iGetNumDof() == 0);
     ASSERT((GetConstLawType() & ConstLawType::VISCOUS) != 0);
     
     SpColVector<GpGradProd, iDim> FTmp(iDim, 1);

     for (index_type i = 1; i <= iDim; ++i) {
          FTmp(i).Reset(F(i));
          
          for (index_type j = 1; j <= iDim; ++j) {
               Eps(j).InsertDeriv(FTmp(i), FDE(i, j));
               EpsPrime(j).InsertDeriv(FTmp(i), FDEPrime(i, j));
          }
     }

     return FTmp;
}

template <class T, class Tder>
template <typename ConstLaw>
inline void ConstitutiveLaw<T, Tder>::UpdateViscoelasticSparse(ConstLaw* pCl, const T& Eps, const T& EpsPrime)
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

	pCl->UpdateViscoelasticTpl(gEps, gEpsPrime, gF);

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
inline void ConstitutiveLaw<T, Tder>::UpdateElasticSparse(ConstLaw* pCl, const T& Eps)
{
	using namespace sp_grad;
	constexpr index_type N = ConstLawHelper<T>::iDim;
	SpColVectorA<SpGradient, N> gEps, gF;

	pCl->Epsilon = Eps;

	for (index_type i = 1; i <= N; ++i)
	{
             gEps(i).Reset(Eps(i), i, 1.);
	}

	pCl->UpdateElasticTpl(gEps, gF);

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

	inline void Update(const T& Eps, const T& EpsPrime = mb_zero<T>()) {
		ASSERT(pConstLaw != NULL);
		pConstLaw->Update(Eps, EpsPrime);
	};

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
	virtual unsigned int iGetNumDof(void) const {
		ASSERT(pConstLaw != NULL);
		return pConstLaw->iGetNumDof();
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
		ASSERT(pConstLaw != NULL);
		return pConstLaw->GetDofType(i);
	};

	/*
	 * Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili
	 */
	virtual unsigned int iGetNumPrivData(void) const {
		return pConstLaw->iGetNumPrivData();
	};

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return pConstLaw->iGetPrivDataIdx(s);
	};

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const {
		return pConstLaw->dGetPrivData(i);
	};

	virtual std::ostream& OutputAppend(std::ostream& out) const {
		return pConstLaw->OutputAppend(out);
	};
	
	virtual void NetCDFOutputAppend(OutputHandler& OH) const {
		return pConstLaw->NetCDFOutputAppend(OH);
	};
	
	virtual void OutputAppendPrepare(OutputHandler& OH, const std::string& name) {
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

/* create/destroy */
extern void InitCL(void);
extern void DestroyCL(void);

#endif /* CONSTLTP_H */

