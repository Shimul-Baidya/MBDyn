/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>
#endif /* HAVE_CONFIG */

#include <stdlib.h>
#include <unistd.h>
#include "ac/iostream"
#include "ac/iomanip"

#include "solman.h"
#include "spmapmh.h"
#include "ccmh.h"
#include "dirccmh.h"
#include "y12wrap.h"
#include "harwrap.h"
#include "mschwrap.h"
#include "umfpackwrap.h"

char *solvers[] = {
#ifdef USE_UMFPACK
	"umfpack",
#endif /* USE_UMFPACK */
#ifdef USE_Y12
	"y12",
#endif /* USE_Y12 */
	NULL
};

static void
usage(void)
{
	std::cerr << "usage: cctest [-d] [-m <solver>]" << std::endl;

	if (!solvers[0]) {
		std::cerr << "\tno solvers available!!!" << std::endl;

	} else {
		std::cerr << "\t<solver>={" << solvers[0];
		for (unsigned i = 1; solvers[i]; i++) {
			std::cerr << "|" << solvers[i];
		}
		std::cerr << "}" << std::endl;
	}

	exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
	SolutionManager *pSM = NULL;
	char *solver = NULL;
	bool dir(false);
	const int size(3);

	while (1) {
		int opt = getopt(argc, argv, "dm:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'd':
			dir = true;
			break;

		case 'm':
			solver = optarg;
			break;

		default:
			usage();
		}
	}

	if (solver == NULL) {
		usage();
	}

	if (strcasecmp(solver, "y12") == 0) {
#ifdef USE_Y12
		if (dir) {
			typedef Y12SparseCCSolutionManager<DirCColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			typedef Y12SparseCCSolutionManager<CColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));
		}
#else /* !USE_Y12 */
		std::cerr << "need --with-y12 to use y12m library" 
			<< std::endl;
		usage();
#endif /* !USE_Y12 */

	} else if (strcasecmp(solver, "umfpack") == 0
			|| strcasecmp(solver, "umfpack3") == 0) {
#ifdef USE_UMFPACK
		if (dir) {
			typedef UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));
		}
#else /* !USE_UMFPACK */
		std::cerr << "need --with-umfpack to use Umfpack library" 
			<< std::endl;
		usage();
#endif /* !USE_UMFPACK */

	} else {
		std::cerr << "unknown solver '" << solver << "'" << std::endl;
		usage();
	}

	std::cout << "using " << solver << " solver..." << std::endl;

	MatrixHandler *pM = pSM->pMatHdl();
	VectorHandler *pV = pSM->pResHdl();
	VectorHandler *px = pSM->pSolHdl();

	doublereal d = 0.;
	int count = 0;

	int p = 2;
	int w = 7 + p;
	std::cout.setf(std::ios::scientific);
	std::cout.precision(p);

retry:;
	pSM->MatrInit(0.);

	try {
		pM = pSM->pMatHdl();

		pM->PutCoef(1, 1, 2.);
		pM->PutCoef(2, 2, 2.);
		pM->PutCoef(3, 3, 1.);
		if (d) {
			pM->PutCoef(1, 2, d);
			pM->PutCoef(2, 1, d);
			pM->PutCoef(2, 3, d);
			pM->PutCoef(3, 2, d);
		}
		
		pV->PutCoef(1, 0.);
		pV->PutCoef(2, 0.);
		pV->PutCoef(3, 1.);

	} catch (MatrixHandler::ErrRebuildMatrix) {
		std::cerr << "need to rebuild matrix..." << std::endl;
		pSM->MatrInitialize();
		goto retry;

	} catch (...) {
		std::cerr << "build failure" << std::endl;
		exit(EXIT_FAILURE);
	}

	try {
		count++;
		std::cout << "solution " << count << "..." << std::endl;
		
		pSM->Solve();

	} catch (...) {
		std::cerr << "solution failure" << std::endl;
		exit(EXIT_FAILURE);
	}

	std::cout
		<< "{x1}   [" << std::setw(w) << 2.
			<< "," << std::setw(w) << d 
			<< "," << std::setw(w) << 0.
		<< "]^-1 {" << std::setw(w) << 0.
		<< "}   {" << std::setw(w) << (*px)(1) << "}" << std::endl
		<< "{x2} = ["<< std::setw(w) << d
			<< "," << std::setw(w) << 2.
			<< "," << std::setw(w) << d 
		<< "]    {" << std::setw(w) << 0.
		<< "} = {" << std::setw(w) << (*px)(2) << "}" << std::endl	
		<< "{x3}   [" << std::setw(w) << 0.
			<< "," << std::setw(w) << d 
			<< "," << std::setw(w) << 1.
		<< "]    {" << std::setw(w) << 1.
		<< "}   {" << std::setw(w) << (*px)(3) << "}" << std::endl;

	switch (count) {
	case 1:
		d = -1.;
		break;

	case 2:
		d = 0.;
		break;

	case 3:
		return 0;
	}

	goto retry;
}

