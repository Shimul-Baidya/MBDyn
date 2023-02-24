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

/* bi-stop drive caller */

#ifndef BISTOPDRIVE_H
#define BISTOPDRIVE_H

/* include generali */

/* include per il debug */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

#include "dataman.h"

/* include del programma */

/* BiStopDriveCaller - begin */

class BiStopDriveCaller : public DriveCaller {
private:
	enum Status { INACTIVE, ACTIVE };

	mutable enum Status m_status;
	const DriveCaller *m_pActivatingCondition;
	const DriveCaller *m_pDeactivatingCondition;

public:
	BiStopDriveCaller(
			const DriveHandler* pDH,
			bool bInitialStatus,
			const DriveCaller *pA,
			const DriveCaller *pD
	) : DriveCaller(pDH),
	m_status(bInitialStatus ? ACTIVE : INACTIVE),
	m_pActivatingCondition(pA), m_pDeactivatingCondition(pD)
	{
		ASSERT(m_pActivatingCondition != 0);
		ASSERT(m_pDeactivatingCondition != 0);
	};

	virtual
	~BiStopDriveCaller(void) {
		SAFEDELETE(m_pActivatingCondition);
		SAFEDELETE(m_pDeactivatingCondition);
	};

	virtual std::ostream&
	Restart(std::ostream& out) const;

	/* Copia */
	virtual DriveCaller* pCopy(void) const;

	inline doublereal dGet(void) const;
	inline doublereal dGet(const doublereal& dVar) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
};

doublereal
BiStopDriveCaller::dGet(void) const
{
	return dGet(pDrvHdl->dGetTime());
}

doublereal
BiStopDriveCaller::dGet(const doublereal& dVar) const
{
	switch (m_status) {
	case INACTIVE:
		if (m_pActivatingCondition->dGet(dVar) == 0.) {
			/* remains inactive: nothing to do */
			return 0.;
		}

		/* activates: change data and ask for jacobian rigeneration */
		m_status = ACTIVE;
		return 1.;

	case ACTIVE:
		if (m_pDeactivatingCondition->dGet(dVar) == 0.) {
			/* remains active: nothing to do */
			return 1.;
		}

		/* disactivates: reset data and ask for jacobian rigeneration */
		m_status = INACTIVE;
		return 0.;

	default:
		/* impossible */
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* BiStopDriveCaller - end */

/* BiStopDCR - start */

struct BiStopDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

/* BiStopDCR - end */

#endif // BISTOPDRIVE_H

