/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Paolo Mantegazza     <paolo.mantegazza@polimi.it>
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

#ifndef TPLDRIVE_IMPL_H
#define TPLDRIVE_IMPL_H

#include <typeinfo>
#include "tpldrive.h"
#include "mbpar.h"

/* ZeroTplDriveCaller - begin */

template <class T>
class ZeroTplDriveCaller : public TplDriveCaller<T> {
public:
        explicit ZeroTplDriveCaller(const DriveHandler* pDH = nullptr)
             :TplDriveCaller<T>(pDH) {
                NO_OP;
        }

        ~ZeroTplDriveCaller(void) {
                NO_OP;
        }

        /* copia */
        virtual TplDriveCaller<T>* pCopy(void) const override {
                typedef ZeroTplDriveCaller<T> dc;
                TplDriveCaller<T>* pDC = 0;

                SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(this->pDrvHdl));

                return pDC;
        }

        /* Scrive il contributo del DriveCaller al file di restart */
        virtual std::ostream& Restart(std::ostream& out) const override {
                return out << "zero";
        }

        virtual std::ostream& Restart_int(std::ostream& out) const override {
                return out;
        }

        inline T Get(const doublereal&) const override {
                return Get();
        }

        inline T Get(void) const override {
                return mb_zero<T>();
        }

        /* this is about drives that are differentiable */
        inline bool bIsDifferentiable(void) const override {
                return true;
        }

        inline T GetP(const doublereal&) const override {
                return GetP();
        }

        inline T GetP() const override {
                return mb_zero<T>();
        }

        inline int getNDrives(void) const override {
                return 0;
        }
};

/* ZeroTplDriveCaller - end */

/* SingleTplDriveCaller - begin */

template <class T>
class SingleTplDriveCaller : public TplDriveCaller<T>, public DriveOwner {
protected:
        T t;

public:
        SingleTplDriveCaller(const DriveHandler* pDH, const DriveCaller* pDC, const T& x)
        : TplDriveCaller<T>(pDH), DriveOwner(pDC), t(const_cast<T&>(x)) {
                NO_OP;
        }

        ~SingleTplDriveCaller(void) {
                NO_OP;
        }

        /* copia */
        virtual TplDriveCaller<T>* pCopy(void) const override {
                typedef SingleTplDriveCaller<T> dc;
                TplDriveCaller<T>* pDC = 0;

                SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(this->pDrvHdl, pGetDriveCaller()->pCopy(), t));

                return pDC;
        }

        /* Scrive il contributo del DriveCaller al file di restart */
        virtual std::ostream& Restart(std::ostream& out) const override {
                out << "single, ",
                Write(out, t, ", ") << ", ";
                return pGetDriveCaller()->Restart(out);
        }

        virtual std::ostream& Restart_int(std::ostream& out) const override {
                Write(out, t, ", ") << ", ";
                return pGetDriveCaller()->Restart(out);
        }

        inline T Get(const doublereal& dVar) const override {
                return t*dGet(dVar);
        }

        /* this is about drives that are differentiable */
        inline bool bIsDifferentiable(void) const override {
                return DriveOwner::bIsDifferentiable();
        }

        inline T GetP(const doublereal& dVar) const override {
                return t * dGetP(dVar);
        }

        inline int getNDrives(void) const override {
                return 1;
        }
};

/* Nota: in caso scalare, viene semplificata la classe in modo da
 *       usare solo il drive senza pesatura che viene assunta unitaria */

template<>
class SingleTplDriveCaller<doublereal>
: public TplDriveCaller<doublereal>, public DriveOwner {
public:
        SingleTplDriveCaller(const DriveHandler* pDH, const DriveCaller* pDC, const doublereal& = 0.)
        : TplDriveCaller<doublereal>(pDH), DriveOwner(pDC) {
                NO_OP;
        }

        ~SingleTplDriveCaller(void) {
                NO_OP;
        }

        /* copia */
        virtual TplDriveCaller<doublereal>* pCopy(void) const override {
                TplDriveCaller<doublereal>* pDC = 0;

                typedef SingleTplDriveCaller<doublereal> dc;
                SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(pDrvHdl, pGetDriveCaller()->pCopy()));

                return pDC;
        }

        /* Scrive il contributo del DriveCaller al file di restart */
        virtual std::ostream& Restart(std::ostream& out) const override {
                out << "single, ";
                return pGetDriveCaller()->Restart(out);
        }

        virtual std::ostream& Restart_int(std::ostream& out) const override {
                return pGetDriveCaller()->Restart(out);
        }

        inline doublereal Get(const doublereal& dVar) const override {
                return dGet(dVar);
        }

        inline bool bIsDifferentiable(void) const override {
                return DriveOwner::bIsDifferentiable();
        }

        inline doublereal GetP(const doublereal& dVar) const override {
                return dGetP(dVar);
        }

        inline int getNDrives(void) const override {
                return 1;
        };
};

template <class T>
TplDriveCaller<T> *
DC2TDC(const DriveHandler* pDH, DriveCaller *pDC, const T& t)
{
        typedef SingleTplDriveCaller<T> STDC;

        TplDriveCaller<T> *p = 0;
        SAFENEWWITHCONSTRUCTOR(p, STDC, STDC(pDH, pDC, t));
        return p;
}


/* SingleTplDriveCaller - end */

#endif // TPLDRIVE_IMPL_H
