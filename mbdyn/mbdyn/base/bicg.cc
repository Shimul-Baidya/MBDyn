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
 
 /* 
  *
  * Copyright (C) 2004
  * Giuseppe Quaranta	<quaranta@aero.polimi.it>
  *
  * classi che implementano la risoluzione del sistema nonlineare 
  */
  
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */
  
#include <solver.h>
#include <bicg.h>  
#ifdef USE_MPI
#include <mbcomm.h>
#include <schsolman.h>
#endif /* USE_MPI */

#include <dofown.h>
#include <unistd.h>
#include <output.h>
#include <ac/math.h>

BiCGStab::BiCGStab(const Preconditioner::PrecondType PType, 
		const integer iPStep,
		doublereal ITol,
		integer MaxIt,
		doublereal etaMx,
		doublereal T) 
: MatrixFreeSolver(PType, iPStep, ITol, MaxIt, etaMx, T)
{
	NO_OP;
}
	
BiCGStab::~BiCGStab(void)
{
	NO_OP;
}

void
BiCGStab::Solve(const NonlinearProblem* pNLP,
		Solver* pS,
		const integer iMaxIter,
		const doublereal& Tol,
		integer& iIterCnt,
		doublereal& dErr ,
		const doublereal& SolTol,
		doublereal& dSolErr)
{
	ASSERT(pNLP != NULL);
	ASSERT(pS != NULL);

	SolutionManager *pSM = pS->pGetSolutionManager();
	
	iIterCnt = 0;
	dSolErr = 0.;

	/* external nonlinear iteration */	
	
	/* riassembla sempre lo jacobiano se l'integratore e' nuovo */
	if (pNLP != pPrevNLP) {
		bBuildMat = true;
	}

	if (!PrecondIter) {
		bBuildMat = true;
	}
	
	pPrevNLP = pNLP;
        pRes = pSM->pResHdl();
	Size = pRes->iGetSize();

	doublereal eta = etaMax;
	doublereal rateo = 0.;
	doublereal Fnorm = 1.;
        doublereal resid;
        doublereal rho_1; 
	doublereal rho_2 = 0.; 
	doublereal alpha;
	doublereal beta;
	doublereal omega;
	VectorHandler* pr;

	/*
	 * these will be resized (actually allocated)
	 * only the first time they're used, unless
	 * the size of the problem changes
	 *
	 * FIXME: need to review this code.
	 */
	rHat.Resize(Size);
	p.Resize(Size);
	pHat.Resize(Size);
	s.Resize(Size);
	sHat.Resize(Size);
	t.Resize(Size);
	v.Resize(Size);
	dx.Resize(Size); 

	integer TotalIter = 0;

#ifdef DEBUG_ITERATIVE
	std::cout << " BiCGStab New Step " <<std::endl;
#endif /* DEBUG_ITERATIVE */

	while (true) {

#ifdef 	USE_EXTERNAL 	
		SendExternal();
#endif /* USE_EXTERNAL */
		
		pRes->Reset(0.);
      		pNLP->Residual(pRes);
		
      		if (outputRes()) {
	 		std::cout << "Residual (" << iIterCnt
				<< "):" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof " << std::setw(8)
					<< iTmpCnt << ": " 
					<< pRes->dGetCoef(iTmpCnt) << std::endl;
			}
      		}

		dErr = MakeResTest(pS, *pRes) * pNLP->TestScale(pResTest);

#ifdef DEBUG_ITERATIVE
		std::cerr << "dErr " << dErr << std::endl;
#endif /* DEBUG_ITERATIVE */

		if (outputIters()) {
#ifdef USE_MPI
			if (MBDynComm.Get_rank() == 0) {
#endif /* USE_MPI */
				std::cout << "\tIteration " << iIterCnt
					<< " " << dErr;
				if (bBuildMat && dErr >= Tol) {
					std::cout << " J";
				}
				std::cout << std::endl;
#ifdef USE_MPI
			}
#endif /* USE_MPI */
		}
		
		if (dErr < Tol) {
	 		return;
      		}
      		if (!isfinite(dErr)) {
			THROW(ErrSimulationDiverged());
		}
		if (iIterCnt > iMaxIter) {
			THROW(NoConvergence());
		}
          	rateo = dErr*dErr/Fnorm;
		Fnorm = dErr*dErr;
		
      		iIterCnt++;

		/* inner iteration to solve the linear system */	
	
		/* BiCGSTAB Iterative solver */
		DEBUGCOUT("Using BiCGStab iterative solver" << std::endl);

        	/* r0 = b- A*x0  but we choose  (x0 = 0)   => r0 = b */
        	/* N.B. *pRes = -F(0) */ 
		
		pr = pRes;
        	doublereal LocTol = eta * dErr;

        	rHat = *pr;

#ifdef DEBUG_ITERATIVE		
		std::cerr << "LocTol " << LocTol << std::endl;
#endif /* DEBUG_ITERATIVE */
		
		rho_1 = dErr*dErr;   /*rhat.InnerProd(r); */
		resid = dErr;
		v.Reset(0.);
		t.Reset(0.);
		p.Reset(0.);
		dx.Reset(0.);
		
		if (bBuildMat) {
			pSM->MatrInit(0.);

rebuild_matrix:;
			try {
      				pNLP->Jacobian(pSM->pMatHdl());

			} catch (MatrixHandler::ErrRebuildMatrix) {
				silent_cout("NewtonRaphsonSolver: "
						"rebuilding matrix..."
						<< std::endl);

				/* need to rebuild the matrix... */
      				pSM->MatrInitialize(0.);
				goto rebuild_matrix;

			} catch (...) {
				throw;
			}

			bBuildMat = false;
			TotalIter = 0;
			TotJac++;

#ifdef DEBUG_ITERATIVE			
			std::cerr << "Jacobian " << std::endl;
#endif /* DEBUG_ITERATIVE */

		}

#ifdef DEBUG_ITERATIVE
		std::cerr << "rho_1 " << rho_1 << std::endl;
#endif /* DEBUG_ITERATIVE */

		int It = 0;
        	while ((resid > LocTol) && (It++ < MaxLinIt)) {
			if (It == 1) {
				p = *pr;
			} else {
				rho_1 = rHat.InnerProd(*pr);

#ifdef DEBUG_ITERATIVE
				std::cerr << "rho_1 " << rho_1 << std::endl;
#endif /* DEBUG_ITERATIVE */

	               		if (fabs(rho_1) < DBL_EPSILON) {
                        		std::cout << "Bi-CGStab Iterative Solver breakdown" 
						<<  " rho_1 = 0 " << std::endl;
					break;
				}
				beta = (rho_1/rho_2) * (alpha/omega);

#ifdef DEBUG_ITERATIVE
				std::cerr << "beta " << beta << std::endl;
#endif /* DEBUG_ITERATIVE */

				p.ScalarAddMul(*pr, p.ScalarAddMul(v, -omega), beta);
			}
			/* right preconditioning */
			pPM->Precond(p, pHat, pSM);
			pNLP->EvalProd(Tau, rHat, pHat, v);
#if 0			
			(pSM->pMatHdl())->MatVecMul(v,pHat);
#endif			
#ifdef DEBUG_ITERATIVE
			std::cout << "v:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof" << std::setw(4) << iTmpCnt << ": " 
					<< v.dGetCoef(iTmpCnt) << std::endl;
			}
#endif /* DEBUG_ITERATIVE */

			alpha = rHat.InnerProd(v);
			alpha = rho_1 / alpha;

#ifdef DEBUG_ITERATIVE
			std::cerr << "alpha " << alpha << std::endl;
#endif /* DEBUG_ITERATIVE */

			s.ScalarAddMul(*pr, v, -alpha);

#ifdef DEBUG_ITERATIVE
			std::cerr << "s.Norm() " << s.Norm() << std::endl;
#endif /* DEBUG_ITERATIVE */

			if ((resid = s.Norm()) < LocTol) {
				dx.ScalarAddMul(pHat, alpha);
				TotalIter++;
				break;
			}
			pPM->Precond(s, sHat, pSM);
			pNLP->EvalProd(Tau, rHat, sHat, t);
#if 0
			(pSM->pMatHdl())->MatVecMul(t,sHat);
#endif
			omega = t.Norm();
			omega = t.InnerProd(s) / omega;

#ifdef DEBUG_ITERATIVE
			std::cerr << "omega " << omega << std::endl;
#endif /* DEBUG_ITERATIVE */

			dx.ScalarAddMul(pHat, alpha);
			dx.ScalarAddMul(sHat, omega);
			pr->ScalarAddMul(s, t, -omega);
			rho_2 = rho_1;
			resid = pr->Norm();

#ifdef DEBUG_ITERATIVE
			std::cerr << "resid " << resid << std::endl;
#endif /* DEBUG_ITERATIVE */

			TotalIter++;
                	if (fabs(omega) < DBL_EPSILON) {
                        	std::cout << "Bi-CGStab Iterative Solver breakdown" 
					<<  " omega = 0 " << std::endl;
				break;
			}
			if (It == MaxLinIt) {
                        	std::cerr << "Iterative inner solver didn't converge."
					<< " Continuing..." << std::endl;
			}
		}
		/* se ha impiegato troppi passi riassembla lo jacobiano */
		
		if (TotalIter >= PrecondIter && PrecondIter) {
			bBuildMat = true;
		}
		/* calcola il nuovo eta */
		
		doublereal etaNew = gamma * rateo;
		doublereal etaBis;
		if ((etaBis = gamma*eta*eta) > .1) {
			etaNew = (etaNew > etaBis) ? etaNew : etaBis;
		}
		eta = (etaNew < etaMax) ? etaNew : etaMax;
		/* prevent oversolving */
		etaBis = .5*Tol/dErr;
		eta = (eta > etaBis) ? eta : etaBis;

#ifdef DEBUG_ITERATIVE
		std::cerr << "eta " << eta << std::endl;
#endif /* DEBUG_ITERATIVE */
		
		if (outputSol()) {      
	 		std::cout << "Solution:" << std::endl;
	 		for (int iTmpCnt = 1; iTmpCnt <= Size; iTmpCnt++) {
	    			std::cout << "Dof " << std::setw(8)
					<< iTmpCnt << ": "
					<< dx.dGetCoef(iTmpCnt) << std::endl;
			}
		}		
		
      		pNLP->Update(&dx);

		dSolErr = MakeSolTest(pS, dx);
		if (outputIters()) {
#ifdef USE_MPI
			if (MBDynComm.Get_rank() == 0) {
#endif /* USE_MPI */
				std::cout << "\t\tSolErr "
					<< dSolErr << std::endl;
#ifdef USE_MPI
			}
#endif /* USE_MPI */
		}
		if (dSolErr < SolTol) {
			THROW(ConvergenceOnSolution());
		}
	}
}

