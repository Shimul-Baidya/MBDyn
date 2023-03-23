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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "bistopdrive.h"

/* BiStopDriveCaller - begin */

/* Copia */
BiStopDriveCaller::DriveCaller*
BiStopDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC,
		BiStopDriveCaller,
		BiStopDriveCaller(pDrvHdl, (m_status == ACTIVE ? true : false),
			m_pActivatingCondition->pCopy(),
			m_pDeactivatingCondition->pCopy()));

	return pDC;
	
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
BiStopDriveCaller::Restart(std::ostream& out) const
{
	return out << "bistop, initial status, "
		<< (m_status == ACTIVE ? "active" : "inactive" ) << ", ",
		m_pActivatingCondition->Restart(out) << ", ",
		m_pDeactivatingCondition->Restart(out);
}

/* this is about drives that are differentiable */
bool
BiStopDriveCaller::bIsDifferentiable(void) const
{
	return false;
}

doublereal
BiStopDriveCaller::dGetP(const doublereal& dVar) const
{
	return 0.;
}

/* BiStopDriveCaller - end */

/* BiStopDCR - start */

DriveCaller *
BiStopDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "bistop");

	/* driver legato ai driver */
	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no driver dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const DriveHandler* pDrvHdl = pDM->pGetDrvHdl();
	DriveCaller *pDC = 0;

	bool bStatus(0);
	if (HP.IsKeyWord("initial" "status")) {
		if (HP.IsKeyWord("active")) {
			bStatus = true;
		} else if (HP.IsKeyWord("inactive")) {
			bStatus = false;
		} else {
		}
	}

	const DriveCaller *pDCa = HP.GetDriveCaller();
	const DriveCaller *pDCd = HP.GetDriveCaller();

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		BiStopDriveCaller,
		BiStopDriveCaller(pDrvHdl, bStatus, pDCa, pDCd));

	return pDC;
}

/* BiStopDCR - end */

