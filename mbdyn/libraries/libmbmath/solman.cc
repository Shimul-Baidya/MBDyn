/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* solution manager */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <string.h>	/* for memset() */

#include "ac/iostream"
#include "ac/iomanip"

#include "solman.h"
#include "matvec3.h"
#include "ls.h"

/* Zero for sparse vector and matrix handlers */
const doublereal dZero = 0.;

/* SolutionDataManager - begin */

SolutionDataManager::~SolutionDataManager(void)
{
	NO_OP;
}

/* SolutionDataManager - end */


/* SolutionManager - begin */

SolutionManager::SolutionManager(void)
: pLS(0)
{
	NO_OP;
}

SolutionManager::~SolutionManager(void)
{
   	if (pLS != NULL) {	
      		SAFEDELETE(pLS);
   	}
}

/* Inizializzatore "speciale" */
void
SolutionManager::MatrInitialize(const doublereal& d)
{
	MatrInit(d);
}

void
SolutionManager::LinkToSolution(const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) {
	NO_OP;
}

/* sposta il puntatore al vettore del residuo */
doublereal *
SolutionManager::ChangeResPoint(doublereal* pd)
{
	ASSERT(pLS);
	return pLS->ChangeResPoint(pd);
}
   
/* sposta il puntatore al vettore della soluzione */
doublereal *
SolutionManager::ChangeSolPoint(doublereal* pd)
{
	ASSERT(pLS);
	return pLS->ChangeSolPoint(pd);
}

/* SolutionManager - end */

