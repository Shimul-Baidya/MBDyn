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
#include <ac/iostream>

#include <solman.h>
#include <spmapmh.h>
#include <ccmh.h>
#include <dirccmh.h>
#include <y12wrap.h>
#include <harwrap.h>
#include <mschwrap.h>
#include <umfpackwrap.h>
#include <superluwrap.h>
#include <lapackwrap.h>
#include <naivewrap.h>

static void
usage(void)
{
	std::cerr << "usage: wraptest [-c] [-d] [-m <solver>] [-s]" << std::endl
		<< "\t<solver>={y12|harwell|meschach|umfpack|superlu|lapack|naive}" << std::endl;
	exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
	SolutionManager *pSM = NULL;
	char *solver =
#if defined(USE_UMFPACK)
		"umfpack"
#elif defined(USE_Y12)
		"y12"
#elif defined(USE_SUPERLU)
		"superlu"
#elif defined(USE_HARWELL)
		"harwell"
#elif defined(USE_MESCHACH)
		"meschach"
#elif defined(USE_LAPACK)
		"lapack"
#else
		"naive"
#if 0
		"no solver!!!"
#endif
#endif /* NO SOLVER !!! */
		;
	bool cc(false);
	bool dir(false);
	unsigned nt = 1;
	bool singular(false);
	const int size(3);

	while (1) {
		int opt = getopt(argc, argv, "cdm:st:");

		if (opt == EOF) {
			break;
		}

		switch (opt) {
		case 'c':
			cc = true;
			break;

		case 'd':
			dir = true;
			break;

		case 'm':
			solver = optarg;
			break;

		case 's':
			singular = true;
			break;

		case 't':
			nt = atoi(optarg);
			if (nt < 1) {
				nt = 1;
			}
			break;

		default:
			usage();
		}
	}

	if (strcasecmp(solver, "lapack") == 0) {
#ifdef USE_LAPACK
		SAFENEWWITHCONSTRUCTOR(pSM, LapackSolutionManager,
				LapackSolutionManager(size));
#else /* !USE_LAPACK */
		std::cerr << "need --with-lapack to use Lapack library dense solver" 
			<< std::endl;
		usage();
#endif /* !USE_LAPACK */

	} else if (strcasecmp(solver, "superlu") == 0) {
#ifdef USE_SUPERLU
		if (dir) {
			typedef SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));

		} else if (cc) {
			typedef SuperLUSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(nt, size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM, SuperLUSparseSolutionManager,
					SuperLUSparseSolutionManager(nt, size));
		}
#else /* !USE_SUPERLU */
		std::cerr << "need --with-superlu to use SuperLU library" 
			<< std::endl;
		usage();
#endif /* !USE_SUPERLU */

	} else if (strcasecmp(solver, "y12") == 0) {
#ifdef USE_Y12
		if (dir) {
			typedef Y12SparseCCSolutionManager<DirCColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else if (cc) {
			typedef Y12SparseCCSolutionManager<CColMatrixHandler<1> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM, Y12SparseSolutionManager,
					Y12SparseSolutionManager(size));
		}
#else /* !USE_Y12 */
		std::cerr << "need --with-y12 to use y12m library" 
			<< std::endl;
		usage();
#endif /* !USE_Y12 */

	} else if (strcasecmp(solver, "harwell") == 0) {
#ifdef USE_HARWELL
		SAFENEWWITHCONSTRUCTOR(pSM, HarwellSparseSolutionManager,
				HarwellSparseSolutionManager(size));
#else /* !USE_HARWELL */
		std::cerr << "need --with-harwell to use HSL library" 
			<< std::endl;
		usage();
#endif /* !USE_HARWELL */

	} else if (strcasecmp(solver, "meschach") == 0) {
#ifdef USE_MESCHACH
		SAFENEWWITHCONSTRUCTOR(pSM, MeschachSparseSolutionManager,
				MeschachSparseSolutionManager(size));
#else /* !USE_MESCHACH */
		std::cerr << "need --with-meschach to use Meschach library" 
			<< std::endl;
		usage();
#endif /* !USE_MESCHACH */
	} else if (strcasecmp(solver, "umfpack") == 0
			|| strcasecmp(solver, "umfpack3") == 0) {
#ifdef USE_UMFPACK
		if (dir) {
			typedef UmfpackSparseCCSolutionManager<DirCColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else if (cc) {
			typedef UmfpackSparseCCSolutionManager<CColMatrixHandler<0> > CCMH;
			SAFENEWWITHCONSTRUCTOR(pSM, CCMH, CCMH(size));

		} else {
			SAFENEWWITHCONSTRUCTOR(pSM,
					UmfpackSparseSolutionManager,
					UmfpackSparseSolutionManager(size));
		}
#else /* !USE_UMFPACK */
		std::cerr << "need --with-umfpack to use Umfpack library" 
			<< std::endl;
		usage();
#endif /* !USE_UMFPACK */

	} else if (strcasecmp(solver, "naive") == 0) {
		SAFENEWWITHCONSTRUCTOR(pSM,
				NaiveSparseSolutionManager,
				NaiveSparseSolutionManager(size));

	} else {
		std::cerr << "unknown solver '" << solver << "'" << std::endl;
		usage();
	}

	pSM->MatrInit(0.);
	
	MatrixHandler *pM = pSM->pMatHdl();
	VectorHandler *pV = pSM->pResHdl();
	VectorHandler *px = pSM->pSolHdl();

	pM->PutCoef(1, 1, 1.);
	pM->PutCoef(2, 2, 2.);
	pM->PutCoef(2, 3, -1.);
	pM->PutCoef(3, 2, -1.);
	if (singular) {
		pM->PutCoef(3, 3, 0.);

	} else {
		pM->PutCoef(3, 3, 3.);
	}
	
	pV->PutCoef(1, 1.);
	pV->PutCoef(2, 1.);
	pV->PutCoef(3, 1.);
	
#ifdef USE_EXCEPTIONS
	try {
#endif /* USE_EXCEPTIONS */
	pSM->Solve();
#ifdef USE_EXCEPTIONS
	} catch (...) {
		exit(EXIT_FAILURE);
	}
#endif /* USE_EXCEPTIONS */
	
	for (int i = 1; i <= size; i++) {
		std::cout << "\tsol[" << i << "] = " << px->dGetCoef(i) 
			<< std::endl;
	}
	
	return 0;
}

