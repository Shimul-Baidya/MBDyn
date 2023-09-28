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
        Copyright (C) 2020(-2023) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef __SP_GRADIENT_FWD_H__INCLUDED__
#define __SP_GRADIENT_FWD_H__INCLUDED__

#include <initializer_list>
#include <vector>

#include "sp_gradient_base.h"
#include "sp_matrix_base_fwd.h"

namespace sp_grad {
     namespace util {
          struct SpGradExpDofMapHelperBase {
               static inline void Reset() {}

               static inline void ResetDofStat() {}

               static inline void InsertDone() {}
          };

          template <typename T, typename Base = SpGradExpDofMapHelperBase>
          class SpGradExpDofMapHelperGeneric: public Base {
          public:
               // Allow us to reuse a SpGradExpDofMap accross multiple expressions.
               // For some elements this may reduce the overhead of compressed
               // evaluation of expressions significantly.
               //
               // Usage:
               // SpGradient var1, var2, var3;
               // var1 = ...
               // var2 = ...
               // SpGradExpDofMapHelper<SpGradient> oDofMap;
               // oDofMap.ResetDofStat();
               // oDofMap.GetDofStat(var1);
               // oDofMap.GetDofStat(var2);
               // oDofMap.Reset();
               // oDofMap.InsertDof(var1);
               // oDofMap.InsertDof(var2);
               // oDofMap.InsertDone();
               // oDofMap.MapAssign(var3, var1 + var2); /* equivalent to var3 = var1 + var2 */

               static inline void GetDofStat(const T&) {}

               static inline void GetDofStat(const T*, const T*) {}

               template <index_type iNumRows, index_type iNumCols>
               static inline void GetDofStat(const SpMatrixBase<T, iNumRows, iNumCols>& A) {}

               static inline void InsertDof(const T&) {}

               static inline void InsertDof(const T*, const T*) {}

               template <index_type iNumRows, index_type iNumCols>
               static inline void InsertDof(const SpMatrixBase<T, iNumRows, iNumCols>& A) {}

               static inline void InitDofMap(T&) {}

               static inline T& MapAssign(T& g, const T& expr) {
                    g = expr;

                    return g;
               }

               static inline const T& MapEval(const T& expr) {
                    return expr;
               }

               template <typename Func>
               static inline T& MapAssignOper(T& g, const T& expr) {
                    g = Func::f(g, expr);

                    return g;
               }

               static inline T& Add(T& g, const T& expr) {
                    g += expr;

                    return g;
               }

               static inline T& Sub(T& g, const T& expr) {
                    g -= expr;

                    return g;
               }

               static inline void GetDofStat() = delete;
               static inline void GetDofMap() = delete;
          };
     }

     template <typename T>
     class SpGradExpDofMapHelper;

     template <>
     class SpGradExpDofMapHelper<doublereal>: public util::SpGradExpDofMapHelperGeneric<doublereal> {
     public:
          static inline index_type iGetLocalSize() { return 0; }
     };

     template <>
     class SpGradExpDofMapHelper<GpGradProd>: public util::SpGradExpDofMapHelperGeneric<GpGradProd, SpGradExpDofMapHelper<doublereal> >  {
     public:
          static inline index_type iGetLocalSize() { return 1; }
     };

     template <>
     class SpGradExpDofMapHelper<SpGradient>: public SpGradExpDofMapHelper<doublereal> {
     public:
          using SpGradExpDofMapHelper<doublereal>::GetDofStat; // Required for expressions with operand types including SpGradient and doublereal
          using SpGradExpDofMapHelper<doublereal>::InsertDof; // Required for expressions with operand types including SpGradient and doublereal

          inline void ResetDofStat();
          inline void Reset();

          inline void InsertDone();

          inline void GetDofStat(const SpGradient& g);

          inline void InsertDof(const SpGradient& g);

          inline void GetDofStat(const SpGradient* pFirst, const SpGradient* const pLast);

          inline void InsertDof(const SpGradient* pFirst, const SpGradient* const pLast);

          inline void InitDofMap(SpGradient& g) const;

          inline index_type iGetLocalSize() const;

          template <index_type iNumRows, index_type iNumCols>
          inline void GetDofStat(const SpMatrixBase<SpGradient, iNumRows, iNumCols>& A);

          template <index_type iNumRows, index_type iNumCols>
          inline void InsertDof(const SpMatrixBase<SpGradient, iNumRows, iNumCols>& A);

          template <typename Expr>
          inline SpGradient& MapAssign(SpGradient& g, const SpGradBase<Expr>& expr) const;

          template <typename Expr>
          inline SpGradient MapEval(const SpGradBase<Expr>& expr) const;

          template <typename Func, typename Expr>
          inline SpGradient& MapAssignOper(SpGradient& g, const SpGradBase<Expr>& expr) const;

          template <typename Expr>
          inline SpGradient& Add(SpGradient& g, const SpGradBase<Expr>& expr) const;

          template <typename Expr>
          inline SpGradient& Sub(SpGradient& g, const SpGradBase<Expr>& expr) const;

          const SpGradDofStat& GetDofStat() const { return oDofStat; }
          const SpGradExpDofMap& GetDofMap() const { return oDofMap; }

     private:
          SpGradDofStat oDofStat;
          SpGradExpDofMap oDofMap;
     };

     template <typename T>
     struct SpGradientTraits;

     template <>
     struct SpGradientTraits<doublereal> {
          static inline void ResizeReset(doublereal& g, doublereal dVal, index_type iSize);
          inline static void InsertDof(doublereal, SpGradExpDofMap&);
          static inline void AddDeriv(doublereal f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) {}

          template <typename Expr, index_type NumRows, index_type NumCols>
          static constexpr inline bool bHaveRefTo(doublereal g, const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

          static constexpr inline doublereal
          dGetValue(doublereal a);
          static constexpr inline index_type
          iGetSize(doublereal a);
          static inline void InsertDeriv(const doublereal& f, doublereal& g, doublereal dCoef) noexcept {}
          static inline void Sort(doublereal);
          static bool bIsUnique(doublereal) { return true; }
          inline static void GetDofStat(doublereal, SpGradDofStat&);
          inline static constexpr doublereal dGetDeriv(doublereal, index_type);
          inline static void ZeroInit(doublereal& g) { g = 0.; }
     };

     template <>
     struct SpGradientTraits<SpGradient> {
          static inline void ResizeReset(SpGradient& g, doublereal dVal, index_type iSize);

          template <typename Expr, index_type NumRows, index_type NumCols>
          static constexpr inline bool bHaveRefTo(const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

          template <index_type NumRows, index_type NumCols>
          static inline bool bHaveRefTo(const SpGradient& g, const SpMatrixBase<SpGradient, NumRows, NumCols>& A);

          static inline void InsertDof(const SpGradient& g, SpGradExpDofMap& oDofMap);

          static inline void AddDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap);

          template <typename Expr>
          static constexpr inline doublereal
          dGetValue(const SpGradBase<Expr>& a);

          template <typename Expr>
          static inline index_type
          iGetSize(const SpGradBase<Expr>& a);

          static inline index_type
          iGetSize(const SpGradient& a);

          static inline void InsertDeriv(const SpGradient& f, SpGradient& g, doublereal dCoef);
          static inline void InsertDeriv(const doublereal& f, SpGradient& g, doublereal dCoef) noexcept {}

          inline static void Sort(SpGradient& g);

          inline static void GetDofStat(const SpGradient& g, SpGradDofStat& s);

          inline static doublereal dGetDeriv(const SpGradient&g, index_type iDof);

          inline static bool bIsUnique(const SpGradient& g);

          inline static void ZeroInit(SpGradient&) {}
     };

     template <>
     struct SpGradientTraits<GpGradProd> {
          static inline doublereal
          dGetValue(const GpGradProd& a);

          static inline index_type
          iGetSize(const GpGradProd& a) { return 1; }

          static inline void InsertDeriv(const GpGradProd& f, GpGradProd& g, doublereal dCoef);
          static inline void InsertDeriv(const doublereal& f, GpGradProd& g, doublereal dCoef) noexcept {}
          static inline void ResizeReset(GpGradProd& g, doublereal dVal, index_type);

          template <typename Expr, index_type NumRows, index_type NumCols>
          static constexpr inline bool bHaveRefTo(const SpMatrixBase<Expr, NumRows, NumCols>& A) { return false; }

          template <index_type NumRows, index_type NumCols>
          static inline bool bHaveRefTo(const GpGradProd& g, const SpMatrixBase<GpGradProd, NumRows, NumCols>& A);

          inline static bool bIsUnique(const GpGradProd& g) { return true; }
          inline static void ZeroInit(GpGradProd&) {}
     };

     class SpGradient: public SpGradBase<SpGradient> {
     public:
          friend util::SpMatrixDataTraits<SpGradient>;

          static constexpr ExprEvalFlags eExprEvalFlags = ExprEvalDuplicate;

          inline SpGradient();

          inline explicit SpGradient(doublereal d);

          inline SpGradient(const SpGradient& g);

          inline SpGradient(SpGradient&& g);

          inline SpGradient(const SpGradExpDofMap& oDofMap);

          inline SpGradient(doublereal dVal, const std::initializer_list<SpDerivRec>& rgDer);

          inline SpGradient(doublereal dVal, const std::vector<SpDerivRec>& rgDer);

          template <typename Expr>
          inline SpGradient(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient(const SpGradBase<Expr>& g, const SpGradExpDofMapHelper<SpGradient>& oDofMap);

          inline ~SpGradient();

          inline SpGradient& operator=(const SpGradient& g);

          inline SpGradient& operator=(SpGradient&& g);

          template <typename Expr>
          inline SpGradient& operator=(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& operator+=(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& operator-=(const SpGradBase<Expr>& g);

          inline SpGradient& operator+=(doublereal b);

          inline SpGradient& operator-=(doublereal b);

          template <typename Expr>
          inline SpGradient& operator*=(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& operator/=(const SpGradBase<Expr>& g);

          inline SpGradient& operator*=(doublereal b);

          inline SpGradient& operator/=(doublereal b);

          inline void Reset(doublereal dVal, const std::initializer_list<SpDerivRec>& rgDer);

          inline void Reset(doublereal dVal, const std::vector<SpDerivRec>& rgDer);

          inline void Reset(doublereal dVal, index_type iDof, doublereal dDer);

          inline void ResetNumeric();

          inline void ResizeReset(doublereal dVal, index_type iSize);

          inline void Scale(doublereal dRowScale, const std::vector<doublereal>& oColScale);

          template <typename Expr>
          inline bool bHaveRefTo(const SpGradBase<Expr>&) const;

          inline bool bHaveRefTo(const SpGradBase<SpGradient>& g) const;

          template <index_type NumRows, index_type NumCols>
          inline bool bHaveRefTo(const SpMatrixBase<SpGradient, NumRows, NumCols>& A) const;

          inline doublereal dGetValue() const;

          inline doublereal dGetDeriv(index_type iDof) const;

          inline void InsertDeriv(SpGradient& g, doublereal dCoef) const;

          inline void InsertDof(SpGradExpDofMap& oExpDofMap) const;

          inline void AddDeriv(SpGradient& g, doublereal dCoef, const SpGradExpDofMap& oDofMap) const;

          inline const SpDerivRec* begin() const;

          inline const SpDerivRec* end() const;

          inline index_type iGetSize() const;

          inline void GetDofStat(SpGradDofStat& s) const;

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset);

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset,
                          const SpGradExpDofMap& oDofMap);

          template <typename AITER, typename BITER>
          inline void
          InnerProduct(AITER pAFirst,
                       AITER pALast,
                       index_type iAOffset,
                       BITER pBFirst,
                       BITER pBLast,
                       index_type iBOffset);

          template <typename Expr>
          inline void Assign(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& MapAssign(const SpGradBase<Expr>& g);

          template <typename Expr>
          inline SpGradient& MapAssign(const SpGradBase<Expr>& g, const SpGradExpDofMap& oDofMap);

          template <typename Func, typename Expr>
          inline SpGradient& AssignOper(const SpGradBase<Expr>& g);

          template <typename Func, typename Expr>
          inline SpGradient& MapAssignOper(const SpGradBase<Expr>& g);

          template <typename Func, typename Expr>
          inline SpGradient& MapAssignOper(const SpGradBase<Expr>& g, const SpGradExpDofMap& oDofMap);

          template <typename Func>
          inline void InitDerivAssign(doublereal f, doublereal df_du, const SpGradExpDofMap& oExpDofMap);

          template <typename Func>
          inline void InitDerivAssign(doublereal f, doublereal df_du, index_type iSizeRes);

          inline void InitDeriv(const SpGradExpDofMap& oExpDofMap);

          inline void InitDofMap(const SpGradExpDofMap& oExpDofMap);

          void Sort();

          inline bool bIsSorted() const;
          inline bool bIsUnique() const;

          void MakeUnique();

#ifdef SP_GRAD_DEBUG
          bool bValid() const;
          bool bCheckUnique() const;
          void PrintValue(std::ostream& os) const;
          void PrintDeriv(std::ostream& os, doublereal dCoef) const;
          static index_type iGetRefCntNullData() { return pGetNullData()->iRefCnt; }
#endif
     private:
          inline void UniqueOwner();

          inline explicit SpGradient(SpDerivData* pData);

          inline void SetValuePreserve(doublereal dVal);

          inline static size_t uGetAllocSize(index_type iSizeRes);

          inline static SpDerivData* pAllocMem(SpDerivData* ptr, index_type iSize);

          void Allocate(index_type iSizeRes, index_type iSizeInit, unsigned uFlags);

          inline void Free();
          inline void Cleanup();

          inline constexpr static bool bRecCompareWithDof(const SpDerivRec& a, index_type b);

          inline constexpr static doublereal AssignMulConst(doublereal a, doublereal b);

          inline constexpr static doublereal AssignDivConst(doublereal a, doublereal b);

          inline void MaybeSort() const;

          inline SpDerivRec* pFindRec(index_type iDof) const;

          inline SpDerivRec* pInsertRec(index_type iDof, doublereal dDer);

          template <typename CONT_TYPE>
          inline void CopyDeriv(doublereal dVal, const CONT_TYPE& rgDer);

          template <doublereal AssOpFn(doublereal, doublereal, doublereal&, doublereal&), typename Expr>
          inline void AssignOper(const SpGradBase<Expr>& g);

          template<doublereal AssOpFn(doublereal, doublereal)>
          inline void AssignOper(doublereal b);

          inline void InnerProductAddDer(const SpGradient& g, const doublereal dVal);

          inline void InnerProductAddDer(doublereal, doublereal) {}

          template <typename ITER>
          inline static index_type InnerProductSize(ITER pFirst,
                                                    ITER pLast,
                                                    index_type iOffset);

          template <typename ITER>
          inline static void InnerProductInsertDof(ITER pFirst,
                                                   ITER pLast,
                                                   index_type iOffset,
                                                   SpGradExpDofMap& oDofMap);

          inline void InnerProductAddDer(const SpGradient& g,
                                         doublereal dVal,
                                         const SpGradExpDofMap& oDofMap);

          inline static void InnerProductAddDer(doublereal,
                                                doublereal,
                                                const SpGradExpDofMap&) {}

          template <typename ITER>
          inline static void InnerProductDofStat(ITER pFirst,
                                                 ITER pLast,
                                                 index_type iOffset,
                                                 SpGradDofStat& s);

          inline static SpDerivData* pGetNullData();

          SpDerivData* pData;
          static SpDerivData oNullData;
     };

     class GpGradProd {
     public:
          static constexpr SpGradCommon::ExprEvalFlags eExprEvalFlags = SpGradCommon::ExprEvalDuplicate;

          constexpr explicit GpGradProd(doublereal dVal = 0., doublereal dDer = 0.)
               :dVal(dVal), dDer(dDer) {
          }

          void Reset(doublereal dNewVal = 0., doublereal dNewDer = 0.) {
               dVal = dNewVal;
               dDer = dNewDer;
          }

          template <typename AITER, typename BITER>
          inline void
          MapInnerProduct(AITER pAFirst,
                          AITER pALast,
                          index_type iAOffset,
                          BITER pBFirst,
                          BITER pBLast,
                          index_type iBOffset);

          constexpr doublereal dGetValue() const {
               return dVal;
          }

          static inline constexpr doublereal dGetValue(doublereal dVal) {
               return dVal;
          }

          static inline doublereal dGetValue(const GpGradProd& g) {
               return g.dGetValue();
          }

          constexpr doublereal dGetDeriv() const {
               return dDer;
          }

          GpGradProd& operator+=(doublereal dCoef) {
               dVal += dCoef;

               return *this;
          }

          GpGradProd& operator-=(doublereal dCoef) {
               dVal -= dCoef;

               return *this;
          }

          GpGradProd& operator*=(doublereal dCoef) {
               dVal *= dCoef;
               dDer *= dCoef;

               return *this;
          }

          GpGradProd& operator/=(doublereal dCoef) {
               dVal /= dCoef;
               dDer /= dCoef;

               return *this;
          }

          inline GpGradProd& operator+=(const GpGradProd& oExpr);

          inline GpGradProd& operator-=(const GpGradProd& oExpr);

          inline GpGradProd& operator*=(const GpGradProd& oExpr);

          inline GpGradProd& operator/=(const GpGradProd& oExpr);

          static constexpr bool bIsScalarConst = false;

          inline void InsertDeriv(GpGradProd& g, doublereal dCoef) const;

          template <index_type NumRows, index_type NumCols>
          static inline bool bHaveRefTo(const SpMatrixBase<GpGradProd, NumRows, NumCols>& A) {
               return false;
          }

          template <typename BinFunc>
          inline GpGradProd& AssignOper(const GpGradProd& oExpr);

#ifdef SP_GRAD_DEBUG
          inline bool bValid() const;
#endif
     private:
          inline void InnerProductAddDer(const GpGradProd& g, const doublereal dCoef);

          inline void InnerProductAddDer(doublereal, doublereal) {}

          doublereal dVal;
          doublereal dDer;
     };
}
#endif
