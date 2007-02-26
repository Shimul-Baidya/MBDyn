/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef HAVE_RUNTIME_LOADING

#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <ac/getopt.h>
#include <ac/iostream>

#include <myassert.h>
#include <solman.h>
#include <y12wrap.h>
#include <harwrap.h>
#include <mschwrap.h>

struct integration_data {
   	doublereal ti;
   	doublereal tf;
   	doublereal dt;
   	doublereal tol;
   	int maxiter;
   	doublereal rho;
};

int method_multistep(const char*, integration_data*, void*, const char*);
int method_hope(const char*, integration_data*, void*, const char*);
int method_cubic(const char*, integration_data*, void*, const char*);
int method_cn(const char*, integration_data*, void*, const char*);

void* get_method_data(int, const char*);

int 
main(int argn, char *const argv[])
{
   	enum {
      		METHOD_UNKNOWN,
		METHOD_MULTISTEP,
		METHOD_HOPE,
		METHOD_CUBIC,
		METHOD_CRANK_NICHOLSON
   	};
	
   	struct s_method {
      		const char* s;
      		int m;
   	} s_m[] = {
		{ "ms",                METHOD_MULTISTEP         },
		{ "hope",              METHOD_HOPE              },
		{ "cubic",             METHOD_CUBIC             },
		{ "crank-nicholson",   METHOD_CRANK_NICHOLSON   },
	
		{ NULL,                METHOD_UNKNOWN           }
   	};
   	int curr_method = METHOD_UNKNOWN;   
   	char* module = "intg-default.so";
   	char* user_defined = NULL;
   	void* method_data = NULL;
   	integration_data d = { 0., 1., 1.e-3, 1.e-6, 10 };

   	/* opzioni */
   	const char optstring[] = "i:m:t:T:n:r:u:h";
   	const struct option longopts[] = { 
		{ "integrator",     required_argument, NULL, int('i') },
		{ "method-data",    required_argument, NULL, int('m') },
		{ "timing",         required_argument, NULL, int('t') },
		{ "tolerance",      required_argument, NULL, int('T') },
		{ "iterations",     required_argument, NULL, int('n') },
		{ "rho",            required_argument, NULL, int('r') },
		{ "user-defined",   required_argument, NULL, int('u') },
		{ "help",           no_argument,       NULL, int('h') },
	
		{ NULL,             no_argument,       NULL, int(0)   }  
   	};

   	while (true) {
		int curropt;
		
		curropt = getopt_long(argn, argv, optstring, longopts, NULL);
		
		if (curropt == EOF) {
			break;
		}
		
      		switch(curropt) {
       		case int('?'):
	  		/* 
			 * std::cerr << "unknown option -" << char(optopt) << std::endl;
			 */
	  		break;
			
       		case int(':'):
	  		std::cerr << "missing parameter for option -" 
				<< optopt << std::endl;
	  		break;
			
       		case int('h'):
	  		std::cerr << "usage: int -[imtTnruh] [module]" << std::endl
	    			<< std::endl
	    			<< " -i,--integrator   : integration method" << std::endl
	    			<< " -m,--method-data  : method-dependent data" << std::endl
	    			<< " -t,--timing       : ti:dt:tf" << std::endl
	    			<< " -T,--tolerance" << std::endl
	    			<< " -n,--maxiter" << std::endl
	    			<< " -r,--rho          : asymptotic radius" << std::endl
	    			<< " -u,--user         : user-defined data" << std::endl
	    			<< std::endl
	    			<< " -h,--help         : print this message and exit" << std::endl;
	  		exit(EXIT_SUCCESS);

       		case int('i'): {
	  		s_method* p = s_m;
	  		while (p->s != NULL) {
	     			if (strcmp(p->s, optarg) == 0) {
					curr_method = p->m;
					break;
	     			}
	     			p++;
	  		}
	  		if (p->s == NULL) {
	     			std::cerr << "unknown integrator " 
					<< optarg << std::endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		break;
       		}
		
       		case int('m'):
	  		method_data = get_method_data(curr_method, optarg);
	  		break;

       		case int('t'): {
	  		char* s = new char[strlen(optarg)+1];
	  		strcpy(s, optarg);
	  		char* p = strrchr(s, ':');
	  		if (p == NULL) {
	     			std::cerr << "syntax: ti:dt:tf" << std::endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		d.tf = atof(p+1);
	  		p[0] = '\0';
	  		p = strrchr(s, ':');
	  		if (p == NULL) {
	     			std::cerr << "syntax: ti:dt:tf" << std::endl;
	     			exit(EXIT_FAILURE);
	  		}
	  		d.dt = atof(p+1);
	  		p[0] = '\0';
	  		d.ti = atof(s);
	  		delete[] s;
	  		break;
       		}
		
       		case int ('T'):
	  		d.tol = atof(optarg);
	  		break;

       		case int ('n'):
	  		d.maxiter = atoi(optarg);
	  		break;

       		case int ('r'):
	  		d.rho = atof(optarg);
	  		break;

       		case int('u'):
	  		user_defined = new char[strlen(optarg)+1];
	  		strcpy(user_defined, optarg);
	  		break;
			
       		default:
	  		/* std::cerr << "unknown option" << std::endl; */
	  		break;
      		}
   	}
   
   	if (argv[optind] != NULL) {
      		module = argv[optind];
      		/* std::cerr << "using module " << module << std::endl; */
   	}
   
   	int rc = 0;
   	switch (curr_method) {
    	case METHOD_MULTISTEP :
      		rc = method_multistep(module, &d, method_data, user_defined);
      		break;
    	case METHOD_HOPE:
      		rc = method_hope(module, &d, method_data, user_defined);
      		break;
    	case METHOD_CUBIC:
      		rc = method_cubic(module, &d, method_data, user_defined);
      		break;
    	case METHOD_CRANK_NICHOLSON:
      		rc = method_cn(module, &d, method_data, user_defined);
      		break;
    	default:
      		std::cerr << "you must select an integration method" << std::endl;
      		exit(EXIT_FAILURE);
   	}
   
   	if (user_defined != NULL) {
      		delete[] user_defined;
      		user_defined = NULL;
   	}
   	return 0;
}

void* 
get_method_data(int curr_method, const char* optarg)
{
   switch (curr_method) {
    default:
      std::cerr << "not implemented yet" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   return NULL;
}


typedef int (*pread_t)(void**, const char*);
typedef int (*pinit_t)(void*, VectorHandler&);
typedef int (*psize_t)(void*);
typedef int (*pgrad_t)(void*, MatrixHandler&, const VectorHandler&, const doublereal&);
typedef int (*pfunc_t)(void*, VectorHandler&, const VectorHandler&, const doublereal&);
typedef std::ostream& (*pout_t)(void*, std::ostream&, const VectorHandler&, const VectorHandler&);
typedef int (*pdestroy_t)(void**);


enum {
   F_READ = 0,
     F_INIT,
     F_SIZE,
     F_GRAD,
     F_FUNC,
     F_OUT,
     F_DESTROY,
     
     F_LAST
};

const char* fnames[] = {
   "read",
     "init",
     "size",
     "grad",
     "func",
     "out",
     "destroy",
     
     NULL
};

const char* fnames_mangled[] = {
     "read__FPPvPCc",
     "init__FPvR13VectorHandler",
     "size__FPv",
     "grad__FPvR13MatrixHandlerRC13VectorHandlerRCd",
     "func__FPvR13VectorHandlerRC13VectorHandlerRCd",
     "out__FPvR7ostreamRC13VectorHandlerT2",
     "destroy__FPPv",
   
     NULL
};

static pread_t ff_read = NULL;
static pinit_t ff_init = NULL;
static psize_t ff_size = NULL;
static pgrad_t ff_grad = NULL;
static pfunc_t ff_func = NULL;
static pout_t ff_out = NULL;
static pdestroy_t ff_destroy = NULL;

int open_module(const char* module) 
{
   void* handle = NULL;
   const char* err = NULL;
   if ((handle = dlopen(module, RTLD_NOW)) == NULL) {
      err = dlerror();
      std::cerr << "dlopen(\"" << module << "\") returned \"" << err << "\"" 
	<< std::endl;
      exit(EXIT_FAILURE);
   }
   
   if ((::ff_read = (pread_t)dlsym(handle, fnames_mangled[F_READ])) == NULL) {
      err = dlerror();
      std::cerr << "dlsym(\"read\") returned \"" << err << "\"" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   if ((::ff_init = (pinit_t)dlsym(handle, fnames_mangled[F_INIT])) == NULL) {
      err = dlerror();
      std::cerr << "dlsym(\"init\") returned \"" << err << "\"" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   if ((::ff_size = (psize_t)dlsym(handle, fnames_mangled[F_SIZE])) == NULL) {
      err = dlerror();
      std::cerr << "dlsym(\"size\") returned \"" << err << "\"" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   if ((::ff_grad = (pgrad_t)dlsym(handle, fnames_mangled[F_GRAD])) == NULL) {
      err = dlerror();
      std::cerr << "dlsym(\"grad\") returned \"" << err << "\"" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   if ((::ff_func = (pfunc_t)dlsym(handle, fnames_mangled[F_FUNC])) == NULL) {
      err = dlerror();
      std::cerr << "dlsym(\"func\") returned \"" << err << "\"" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   if ((::ff_out = (pout_t)dlsym(handle, fnames_mangled[F_OUT])) == NULL) {
      err = dlerror();
      std::cerr << "dlsym(\"out\") returned \"" << err << "\"" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   if ((::ff_destroy = (pdestroy_t)dlsym(handle, fnames_mangled[F_DESTROY])) == NULL) {
      err = dlerror();
      std::cerr << "dlsym(\"destroy\") returned \"" << err << "\"" << std::endl;
      exit(EXIT_FAILURE);
   }
   
   return 0;
}

void flip(VectorHandler** ppX, VectorHandler** ppXP,
	  VectorHandler** ppXm1, VectorHandler** ppXPm1,
	  VectorHandler** ppXm2, VectorHandler** ppXPm2)
{
   VectorHandler* p = *ppXm2;
   *ppXm2 = *ppXm1;
   *ppXm1 = *ppX;
   *ppX = p;
   
   p = *ppXPm2;
   *ppXPm2 = *ppXPm1;
   *ppXPm1 = *ppXP;
   *ppXP = p;
}

int method_multistep(const char* module, integration_data* d, 
		     void* method_data, const char* user_defined)
{ 
   open_module(module);
 
   // prepara i dati
   void* p_data = NULL;
   (*::ff_read)(&p_data, user_defined);
   
   // prepara le strutture dati per il calcolo
   int size = (*::ff_size)(p_data);
   MyVectorHandler v0(size);
   MyVectorHandler v1(size);
   MyVectorHandler v2(size);
   MyVectorHandler v3(size);
   MyVectorHandler v4(size);
   MyVectorHandler v5(size);
   
   VectorHandler* pX = &v0;
   VectorHandler* pXP = &v1;
   VectorHandler* pXm1 = &v2;
   VectorHandler* pXPm1 = &v3;
   VectorHandler* pXm2 = &v4;
   VectorHandler* pXPm2 = &v5;
   pX->Reset();
   pXP->Reset();
   pXm1->Reset();
   pXPm1->Reset();
   pXm2->Reset();
   pXPm2->Reset();
   
   doublereal* pd = new doublereal[size*size];
   doublereal** ppd = new doublereal*[size];
   FullMatrixHandler J(pd, ppd, size*size, size, size);
   MyVectorHandler R(size);

#if defined(USE_Y12)
   Y12SparseSolutionManager sm(size, size*(size+1)+1, 1.);
#elif defined(USE_HARWELL)
   HarwellSparseSolutionManager sm(size, size*(size+1)+1, 1.);
#elif defined(USE_MESCHACH)
   MeschachSparseSolutionManager sm(size, size*(size+1)+1, 1.);
#endif
   
   MatrixHandler& Jac = *sm.pMatHdl();
   VectorHandler& Res = *sm.pResHdl();   
   VectorHandler& Sol = *sm.pSolHdl();   

   // paramteri di integrazione
   doublereal ti = d->ti;
   doublereal dt = d->dt;
   doublereal tf = d->tf;
   
   doublereal tol = d->tol;
   int maxiter = d->maxiter;
   
   // coefficienti del metodo
   doublereal rho = d->rho;
   doublereal beta = (4.*rho*rho-(1.-rho)*(1.-rho))/(4.-(1.-rho)*(1.-rho));
   doublereal delta = (1.-rho)*(1.-rho)/(2.*(4.-(1.-rho)*(1.-rho)));
   doublereal a1 = 1.-beta;
   doublereal a2 = beta;
   doublereal b0 = delta+.5;
   doublereal b1 = .5*beta+.5-2.*delta;
   doublereal b2 = .5*beta+delta;

   doublereal t = ti;
   
   // inizializza la soluzione
   (*::ff_init)(p_data, *pX);
   (*::ff_func)(p_data, *pXP, *pX, t);
   for (int k = 1; k <= size; k++) {
      doublereal x = pX->dGetCoef(k);
      doublereal xp = pXP->dGetCoef(k);
      pXPm1->PutCoef(k, xp);
      pXm1->PutCoef(k, x-dt*xp);
   }
   
   // output iniziale
   std::cout << ti << " " << 0. << " ";
   (*::ff_out)(p_data, std::cout, *pX, *pXP) << std::endl;
   
   flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
   
   while (t < tf) {
      t += dt;
      // predict
      for (int k = size; k > 0; k--) {
	 doublereal xm1 = pXm1->dGetCoef(k);
	 doublereal xPm1 = pXPm1->dGetCoef(k);
	 doublereal xm2 = pXm2->dGetCoef(k);
	 doublereal xPm2 = pXPm2->dGetCoef(k);
	 doublereal x = -4.*xm1+5.*xm2+dt*(4.*xPm1+2.*xPm2);
	 pX->PutCoef(k, x);
	 R.PutCoef(k, a1*xm1+a2*xm2+dt*(b1*xPm1+b2*xPm2));
      }
      
      // test
      int j = 0;
      doublereal test;
      doublereal coef = dt*b0;
      do {
	 (*::ff_func)(p_data, *pXP, *pX, t);
	 for (int k = 1; k <= size; k++) {
	    doublereal x = pX->dGetCoef(k);
	    doublereal xP = pXP->dGetCoef(k);
	    Res.PutCoef(k, R.dGetCoef(k)-x+coef*xP);
	 }

	 test = Res.Norm();
	 if (test < tol) {
	    break;
	 }      
	 if (++j > maxiter) {
	    std::cerr << "current iteration " << j 
	      << " exceeds max iteration number " << maxiter << std::endl;
	    exit(EXIT_FAILURE);
	 }
	 
	 // correct
	 sm.MatrReset();
	 J.Reset();
	 (*::ff_grad)(p_data, J, *pX, t);
	 for (int k = 1; k <= size; k++) {
	    for (int l = 1; l <= size; l++) {
	       Jac.PutCoef(k, l, -coef*J.dGetCoef(k, l));
	    }
	    Jac.IncCoef(k, k, 1.);
	 }      
	 sm.Solve();
	 
	 // update
	 for (int k = size; k > 0; k--) {
	    doublereal dx = Sol.dGetCoef(k);	
	    pX->IncCoef(k, dx);
	 }	 
      } while (true);
      
      // output
      std::cout << t << " " << test << " ";
      (*::ff_out)(p_data, std::cout, *pX, *pXP) << std::endl;
      
      flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);            
   }
   
   (*::ff_destroy)(&p_data);
   delete[] pd;
   delete[] ppd;

   return 0;
}

int method_hope(const char* module, integration_data* d, 
		void* method_data, const char* user_defined)
{
   open_module(module);
   
   std::cerr << __FUNCTION__ << "not implemented yet!" << std::endl;
   exit(EXIT_FAILURE);
   
   return 0;
}

int method_cubic(const char* module, integration_data* d, 
		 void* method_data, const char* user_defined)
{
   open_module(module);

   // prepara i dati
   void* p_data = NULL;
   (*::ff_read)(&p_data, user_defined);
   
   // prepara le strutture dati per il calcolo
   int size = (*::ff_size)(p_data);
   MyVectorHandler v0(size);
   MyVectorHandler v1(size);
   MyVectorHandler v2(size);
   MyVectorHandler v3(size);
   MyVectorHandler v4(size);
   MyVectorHandler v5(size);
   
   VectorHandler* pX = &v0;
   VectorHandler* pXP = &v1;
   VectorHandler* pXm1 = &v2;
   VectorHandler* pXPm1 = &v3;
   VectorHandler* pXm2 = &v4;
   VectorHandler* pXPm2 = &v5;
   pX->Reset();
   pXP->Reset();
   pXm1->Reset();
   pXPm1->Reset();
   pXm2->Reset();
   pXPm2->Reset();
   
   doublereal* pd = new doublereal[2*size*size];
   doublereal** ppd = new doublereal*[2*size];
   FullMatrixHandler Jz(pd, ppd, size*size, size, size);
   FullMatrixHandler J0(pd+size*size, ppd+size, size*size, size, size);
   MyVectorHandler Xz(size);
   MyVectorHandler XPz(size);
   
#if defined(USE_Y12)
   Y12SparseSolutionManager sm(size, size*(size+1)+1, 1.);
#elif defined(USE_HARWELL)
   HarwellSparseSolutionManager sm(size, size*(size+1)+1, 1.);
#elif defined(USE_MESCHACH)
   MeschachSparseSolutionManager sm(size, size*(size+1)+1, 1.);
#endif
   
   MatrixHandler& Jac = *sm.pMatHdl();
   VectorHandler& Res = *sm.pResHdl();   
   VectorHandler& Sol = *sm.pSolHdl();   

   // paramteri di integrazione
   doublereal ti = d->ti;
   doublereal dt = d->dt;
   doublereal tf = d->tf;
   
   doublereal tol = d->tol;
   int maxiter = d->maxiter;
   
   // coefficienti del metodo
   doublereal rho = d->rho;
   doublereal z = -rho/(1.+rho);
   doublereal w1 = (2.+3.*z)/(6.*(1.+z));
   doublereal wz = -1./(6.*z*(1.+z));
   doublereal w0 = (1.+3.*z)/(6.*z);
   doublereal m0 = 1.-z*z*(3.+2.*z);
   doublereal m1 = z*z*(3.+2.*z);
   doublereal n0 = z*(1.+z)*(1.+z);
   doublereal n1 = z*z*(1.+z);

   doublereal t = ti;
   
   // inizializza la soluzione
   (*::ff_init)(p_data, *pX);
   (*::ff_func)(p_data, *pXP, *pX, t);
   for (int k = 1; k <= size; k++) {
      doublereal x = pX->dGetCoef(k);
      doublereal xp = pXP->dGetCoef(k);
      pXPm1->PutCoef(k, xp);
      pXm1->PutCoef(k, x-dt*xp);
   }
   
   // output iniziale
   std::cout << ti << " " << 0. << " ";
   (*::ff_out)(p_data, std::cout, *pX, *pXP) << std::endl;
   
   flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);
   
   while (t < tf) {
      t += dt;
      // predict
      for (int k = 1; k <= size; k++) {
	 doublereal xm1 = pXm1->dGetCoef(k);
	 doublereal xPm1 = pXPm1->dGetCoef(k);
	 doublereal xm2 = pXm2->dGetCoef(k);
	 doublereal xPm2 = pXPm2->dGetCoef(k);
	 doublereal x = -4.*xm1+5.*xm2+dt*(4.*xPm1+2.*xPm2);
	 pX->PutCoef(k, x);
      }
      
      // test
      int j = 0;
      doublereal test;      
      do {
	 pXP->Reset();
	 (*::ff_func)(p_data, *pXP, *pX, t);
	 for (int k = 1; k <= size; k++) {
	    doublereal x = pX->dGetCoef(k);
	    doublereal xP = pXP->dGetCoef(k);
	    doublereal xm1 = pXm1->dGetCoef(k);
	    doublereal xPm1 = pXPm1->dGetCoef(k);
	    doublereal xz = m0*x+m1*xm1+dt*(n0*xP+n1*xPm1);
	    Xz.PutCoef(k, xz);
	 }
	 XPz.Reset();
	 (*::ff_func)(p_data, XPz, Xz, t+z*dt);
	 for (int k = 1; k <= size; k++) {
	    doublereal d = dt*(
			       w1*pXPm1->dGetCoef(k)
			       +wz*XPz.dGetCoef(k)
			       +w0*pXP->dGetCoef(k)
			       )-(
				  pX->dGetCoef(k)
				  -pXm1->dGetCoef(k)
				  );
	    Res.PutCoef(k, d);
	 }

	 test = Res.Norm();
	 if (test < tol) {
	    break;
	 }      
	 if (++j > maxiter) {
	    std::cerr << "current iteration " << j 
	      << " exceeds max iteration number " << maxiter << std::endl;
	    exit(EXIT_FAILURE);
	 }
	 
	 // correct
	 sm.MatrReset();
	 Jz.Reset();
	 J0.Reset();      
	 (*::ff_grad)(p_data, Jz, Xz, t+z*dt);
	 (*::ff_grad)(p_data, J0, *pX, t);	 
	 for (int k = 1; k <= size; k++) {
	    for (int l = 1; l <= size; l++) {
	       doublereal d = 0.;
	       for (int m = 1; m <= size; m++) {
		  d += Jz.dGetCoef(k, m)*J0.dGetCoef(m, l);
	       }
	       d = -dt*(wz*(Jz.dGetCoef(k, l)+dt*n0*d)+w0*J0.dGetCoef(k, l));
	       Jac.PutCoef(k, l, d);
	    }
	    Jac.IncCoef(k, k, 1.);
	 }
	 sm.Solve();
	 
	 // update
	 for (int k = size; k > 0; k--) {
	    doublereal dx = Sol.dGetCoef(k);
	    pX->IncCoef(k, dx);
	 }
      } while (true);
      
      // output
      std::cout << t << " " << test << " ";
      (*::ff_out)(p_data, std::cout, *pX, *pXP) << std::endl;
      
      flip(&pX, &pXP, &pXm1, &pXPm1, &pXm2, &pXPm2);            
   }
   
   (*::ff_destroy)(&p_data);
   delete[] pd;
   delete[] ppd;

   return 0;
}

int method_cn(const char* module, integration_data* d,
	      void* method_data, const char* user_defined)
{
   open_module(module);
   
   std::cerr << __FUNCTION__ << "not implemented yet!" << std::endl;
   exit(EXIT_FAILURE);

   return 0;
}

#else /* ! HAVE_RUNTIME_LOADING */

#include <ac/iostream>
#include <stdlib.h>

int
main(void)
{
	std::cerr << "Need dynamic load capabilities" << std::endl;
   	exit(EXIT_FAILURE);
}

#endif /* ! HAVE_RUNTIME_LOADING */
