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

/* Drive che usa  i gradi di liberta' */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "privdrive.h"
#include "elem.h"
#include "node.h"


PrivDriveCaller::PrivDriveCaller(const DriveHandler* pDH, 
		const DriveCaller* pDC,
		const SimulationEntity *p, unsigned int i, const std::string& s)
: DriveCaller(pDH), DriveOwner(pDC), pSE(p), iIndex(i), sIndexName(s), bWarnMsg(true)
{
	NO_OP;
};     	

PrivDriveCaller::~PrivDriveCaller(void)
{
	NO_OP;
}

void PrivDriveCaller::PrintWarningMsgOnce() const
{
	const char* pszSimEntityType = dynamic_cast<const Elem*>(pSE) ? "element" : "node";

	char buf[BUFSIZ] = "";
	const WithLabel* pWL = dynamic_cast<const WithLabel*>(pSE);
	if (pWL) {
		snprintf(buf, sizeof(buf), "(%u)", pWL->GetLabel());
	}

	const char *type = 0;
	if (dynamic_cast<const Elem*>(pSE)) {
		type = psElemNames[dynamic_cast<const Elem*>(pSE)->GetElemType()];
	} else {
		type = psNodeNames[dynamic_cast<const Node*>(pSE)->GetNodeType()];
	}

	// should this be "pedantic"?
	pedantic_cerr("warning, possible improper call of "
		<< pszSimEntityType << " drive with real argument; \"dVar\" is ignored and private data \""
		<< sIndexName << "\" of " << type << buf << " is returned instead\n");
	bWarnMsg = false; // Let us assume that it is sufficient to print this annoying warning only once for each drive caller!
}

/* Copia */
DriveCaller*
PrivDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = NULL;

	SAFENEWWITHCONSTRUCTOR(pDC, 
			PrivDriveCaller,
			PrivDriveCaller(pDrvHdl,
				pGetDriveCaller()->pCopy(),
				pSE, iIndex, sIndexName));
   
	return pDC;
}


/* Restart */
std::ostream&
PrivDriveCaller::Restart(std::ostream& out) const
{
	const Elem *pElem = dynamic_cast<const Elem *>(pSE);
	const Node *pNode = dynamic_cast<const Node *>(pSE);

	if (pElem != 0) {
		out << " element, " 
			<< pElem->GetLabel() << ", "
			<< psReadElemsElems[pElem->GetElemType()] << ", ";

	} else if (pNode != 0) {
		out << " node, " 
			<< pNode->GetLabel() << ", "
			<< psReadNodesNodes[pNode->GetNodeType()] << ", ";

	} else {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!sIndexName.empty()) {
		out << "string, " << sIndexName;

	} else if (pElem->iGetNumPrivData() > 1) {
		out << "index, " << iIndex;
	}

	out << ", ";

	return DriveOwner::pGetDriveCaller()->Restart(out);
}

unsigned int
PrivDriveCaller::iGetIndex(void) const {
	return iIndex;
}

const SimulationEntity *
PrivDriveCaller::iGetSE(void) const {
	return pSE;
}


