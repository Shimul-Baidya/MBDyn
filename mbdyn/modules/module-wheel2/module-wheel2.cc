/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include <ac/iostream>
#include <ac/float.h>

#include <loadable.h>

/*
 * Usage:
 *
 *	loadable: <label>, <module name>, help [ , <module data> ] ;
 */

struct module_wheel {
	/*
	 * Ruota
	 */
	StructNode *pWheel;

	/*
	 * Direzione asse ruota
	 */
	Vec3 WheelAxle;
	
	/*
	 * Terreno
	 */
	StructNode *pGround;

	/*
	 * Posizione e orientazione del terreno
	 */
	Vec3 GroundPosition;
	Vec3 GroundDirection;

	/*
	 * Geometria ruota
	 */
	doublereal dRadius;
	doublereal dInternalRadius;
	doublereal dVolCoef;

	doublereal dRefArea;
	doublereal dRNP;		/* R+nG'*pG */
	doublereal dV0;

	/*
	 * Proprieta' pneumatico
	 */
	doublereal dP0;
	doublereal dGamma;
	doublereal dHystVRef;

	/*
	 * Attrito
	 */
	flag fSlip;
	DriveCaller *pMuX0;
	DriveCaller *pMuY0;
	DriveCaller *pMuY1;
	doublereal dvThreshold;

	/*
	 * Output
	 */
	Vec3 F;
	Vec3 M;
	doublereal dInstRadius;
	doublereal dDeltaL;
	doublereal dVn;
	doublereal dSr;
	doublereal dAlpha;
	doublereal dAlphaThreshold;
	doublereal dMuX;
	doublereal dMuY;
};

/* funzioni di default */
static void*
read(LoadableElem* pEl,
		DataManager* pDM,
		MBDynParser& HP,
		const DriveHandler* pDH)
{
	DEBUGCOUTFNAME("read");
	
	/* allocation of user-defined struct */
	module_wheel* p = NULL;
	SAFENEW(p, module_wheel);

	/*
	 * help
	 */
	if (HP.IsKeyWord("help")) {
		std::cout <<
"									\n"
"Module: 	wheel2							\n"
"Author: 	Stefania Gualdi <gualdi@aero.polimi.it>			\n"
"		Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"Connects 2 structural nodes:						\n"
"     -	Wheel								\n"
"     -	Ground								\n"
"									\n"
"Note: 									\n"
"     -	The Axle and the Wheel structural nodes must be connected 	\n"
"	by a joint that allows relative rotations only about 		\n"
"	one axis (the axle)						\n"
"     -	The center of the wheel is assumed coincident with 		\n"
"	the position of the wheel structural node			\n"
"     -	The Ground structural node supports a plane defined		\n"
"	a point and a direction orthogonal to the plane (future 	\n"
"	versions might use an arbitrary, deformable surface)		\n"
"     -	The forces are applied at the \"contact point\", that 		\n"
"	is defined according to geometrical properties 			\n"
"	of the system and according to the relative position 		\n"
"	and orientation of the Wheel and Ground structural nodes	\n"
"									\n"
"     -	Input:								\n"
"		<wheel structural node label> ,				\n"
"		<wheel axle direction> ,				\n"
"		<ground structural node label> ,			\n"
"		<reference point position of the ground plane> ,	\n"
"		<direction orthogonal to the ground plane> ,		\n"
"		<wheel radius> ,					\n"
"		<torus radius> ,					\n"
"		<volume coefficient (black magic?)> ,			\n"
"		<tire pressure> ,					\n"
"		<tire polytropic exponent> ,				\n"
"		<reference velocity for tire hysteresis>		\n"
"		[ slip ,						\n"
"		<longitudinal friction coefficient drive>		\n"
"		<lateral friction coefficient drive for s.r.=0>		\n"
"		<lateral friction coefficient drive for s.r.=1>		\n"
"		[ , threshold , <slip ratio velocity threshold> , 	\n"
"			<slip angle velocity threshold> ] ]		\n"
"									\n"
"     -	Output:								\n"
"		1)	element label					\n"
"		2-4)	tire force in global reference frame		\n"
"		5-7)	tire couple in global reference frame		\n"
"		8)	effective radius				\n"
"		9)	tire radial deformation				\n"
"		10)	tire radial deformation velocity		\n"
"		11)	slip ratio					\n"
"		12)	slip angle					\n"
"		13)	longitudinal friction coefficient		\n"
"		14)	lateral friction coefficient			\n"
			<< std::endl;

		if (!HP.fIsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			THROW(NoErr());
		}
	}

	/*
	 * leggo la ruota
	 */
	p->pWheel = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);

	/*
	 * leggo l'orientazione dell'asse ruota nel sistema locale
	 */
	ReferenceFrame RF = ReferenceFrame(p->pWheel);
	p->WheelAxle = HP.GetVecRel(RF);
	
	/*
	 * leggo il terreno
	 */
	p->pGround = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);
	
	/*
	 * leggo posizione ed orientazione del terreno nel sistema del nodo
	 */
	RF = ReferenceFrame(p->pGround);
	p->GroundPosition = HP.GetPosRel(RF);
	p->GroundDirection = HP.GetVecRel(RF);

	/*
	 * normalizzo l'orientazione del terreno
	 */
	doublereal d = p->GroundDirection.Dot();
	if (d <= DBL_EPSILON) {
		std::cerr << "null direction at line " << HP.GetLineData() << std::endl;
		THROW(DataManager::ErrGeneric());
	}
	p->GroundDirection /= sqrt(d);

	/*
	 * Geometria ruota
	 */
	p->dRadius = HP.GetReal();
	p->dInternalRadius = HP.GetReal();
	p->dVolCoef = HP.GetReal();
	
	/*
	 * Area di riferimento
	 */
	p->dRefArea = 3.7*p->dInternalRadius*sqrt(
			p->dInternalRadius*(2.*p->dRadius-p->dInternalRadius)
			);

	/*
	 * termine per il calcolo di Delta L
	 */
	p->dRNP = p->dRadius+p->GroundPosition*p->GroundDirection;

	/*
	 * Volume di riferimento
	 */
	p->dV0 = 2.*M_PI*(p->dRadius-p->dInternalRadius)
		*M_PI*p->dInternalRadius*p->dInternalRadius*p->dVolCoef;

	/*
	 * Proprieta' pneumatico
	 */
	p->dP0 = HP.GetReal();
	p->dGamma = HP.GetReal();
	p->dHystVRef = HP.GetReal();

	/*
	 * Attrito
	 */
	p->fSlip = 0;
	if (HP.IsKeyWord("slip")) {
		p->fSlip = 1;

		/*
		 * Parametri di attrito
		 */
		p->pMuX0 = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
		p->pMuY0 = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
		p->pMuY1 = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
	
		p->dvThreshold = 0.;
		p->dAlphaThreshold = 0.;
		if (HP.IsKeyWord("threshold")) {
			p->dvThreshold = HP.GetReal();
			if (p->dvThreshold < 0.) {
				std::cerr << "illegal velocity threshold "
					<< p->dvThreshold << " at line "
					<< HP.GetLineData() << std::endl;
				p->dvThreshold = fabs(p->dvThreshold);
			}

			p->dAlphaThreshold = HP.GetReal();
			if (p->dvThreshold < 0.) {
				std::cerr << "illegal slip angle threshold "
					<< p->dAlphaThreshold << " at line "
					<< HP.GetLineData() << std::endl;
				p->dAlphaThreshold = fabs(p->dAlphaThreshold);
			}
		}
	}
	
	return (void *)p;
}

static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_dof");
	return 0;
}

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("output");
	
	if (pEl->fToBeOutput()) {
		module_wheel* p = (module_wheel *)pEl->pGetData();      
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << pEl->GetLabel() << " ",
			p->F.Write(out, " ") << " ",
			p->M.Write(out, " ") << " ";
		if (p->dDeltaL > 0.) {
			out << p->dInstRadius << " "
				<< p->dDeltaL << " ";
		} else {
			out << p->dRadius << " "
				<< 0. << " ";
		}
		out << p->dVn << " " 
			<< p->dSr << " "
			<< 180./M_PI*p->dAlpha << " "
			<< p->dMuX << " "
			<< p->dMuY << std::endl;
	}
}

static std::ostream&
restart(const LoadableElem* pEl, std::ostream& out)
{
	DEBUGCOUTFNAME("restart");
	return out << "not implemented yet;" << std::endl;
}

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
	DEBUGCOUTFNAME("work_space_dim");
	*piNumRows = 12;
	*piNumCols = 12;
}

static VariableSubMatrixHandler& 
ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
	DEBUGCOUTFNAME("ass_jac");   
	WorkMat.SetNullMatrix();
	
	return WorkMat;
}

static SubVectorHandler& 
ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ass_res");
	
	module_wheel* p = (module_wheel *)pEl->pGetData();

	/*
	 * Orientazione del terreno nel sistema assoluto
	 */
	Vec3 n = p->pGround->GetRCurr()*p->GroundDirection;
	
	/*
	 * Distanza Wheel Ground nel sistema assoluto
	 */
	Vec3 x = p->pWheel->GetXCurr()-p->pGround->GetXCurr();

	/*
	 * se dDeltaL > 0 c'e' contatto, altrimenti no
	 */
	doublereal dDeltaL = p->dRNP - x*n;

	/*
	 * Reset dati per output
	 */
	p->dDeltaL = dDeltaL;
	p->dInstRadius = p->dRadius-dDeltaL;
	
	p->dSr = 0.;
	p->dAlpha = 0.;

	p->dMuX = 0.;
	p->dMuY = 0.;

	if (dDeltaL < 0.) {
		
		p->F = Zero3;
		p->M = Zero3;
		
		/*
		 * Non assemblo neppure il vettore ;)
		 */
		WorkVec.Resize(0);

		return WorkVec;
	}
	
	/*
	 * Dimensiono il vettore
	 */
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
	WorkVec.Resize(iNumRows);
	WorkVec.Reset(0.);

	integer iGroundFirstMomIndex = p->pGround->iGetFirstMomentumIndex();
	integer iWheelFirstMomIndex = p->pWheel->iGetFirstMomentumIndex();

	/*
	 * Indici equazioni
	 */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.fPutRowIndex(iCnt, iGroundFirstMomIndex+iCnt);
		WorkVec.fPutRowIndex(6+iCnt, iWheelFirstMomIndex+iCnt);
	}

	/*
	 * Velocita' tra Wheel (nell'asse)
	 * e Ground nel sistema assoluto
	 */
	Vec3 va = p->pWheel->GetVCurr()
		-p->pGround->GetVCurr()-(p->pGround->GetWCurr()).Cross(
			p->pGround->GetRCurr()*p->GroundPosition
			);
	/*
	 * Velocita' tra Wheel (nel punto di contatto) 
	 * e Ground nel sistema assoluto
	 */
	Vec3 v = va-(p->pWheel->GetWCurr()).Cross(n*p->dInstRadius);
	
	/*
	 * Componente normale al terreno della velocita'
	 * (positiva se la ruota si allontana dal terreno)
	 */
	p->dVn = n*v;
	
	/*
	 * Stima dell'area di contatto
	 */
	doublereal dA = p->dRefArea*(dDeltaL/p->dInternalRadius);
	
	/*
	 * Stima del volume di compenetrazione tra ruota e terreno
	 */
	doublereal dDeltaV = .5*dA*dDeltaL;

	/*
	 * Stima della pressione del pneumatico (adiabatica)
	 */
	doublereal dP = p->dP0*pow(p->dV0/(p->dV0-dDeltaV), p->dGamma);

	/*
	 * Il punto di applicazione della forza e' xW - R * nG ;
	 * per la forza normale si puo' usare anche la posizione
	 * della ruota, nell'ipotesi che il punto di contatto
	 * stia nell'intersezione della normale al ground passante
	 * per l'asse ruota.
	 *
	 * FIXME: perche' dRadius invece di dInstRadius?
	 */
	Vec3 pc = p->pWheel->GetXCurr()-(n*p->dInstRadius);

	/*
	 * Forza
	 */
	doublereal dFn = (dA*dP*(1.-tanh(p->dVn/p->dHystVRef)));
	p->F = n*dFn;

	if (p->fSlip) {
	
		/*
		 * Direzione di "avanzamento": asse ruota cross normale
		 * al ground
		 */
		Vec3 fwd = (p->pWheel->GetRCurr()*p->WheelAxle).Cross(n);
		doublereal d = fwd.Dot();
		if (d < DBL_EPSILON) {
			std::cerr << "wheel axle is orthogonal to the ground"
				<< std::endl;
			THROW(DataManager::ErrGeneric());
		}
		fwd /= sqrt(d);

		/*
		 * Slip ratio
		 */
		doublereal dvx = fwd.Dot(v);
		doublereal sgn = copysign(1., dvx);
		doublereal dfvx = fabs(dvx);
		doublereal dvax = fwd.Dot(va);
		doublereal dfvax = fabs(dvax);

		/*
		 * FIXME: se la vax va a zero (perche' il velivolo si e'
		 * fermato, ad esempio) lo "sleep" ratio deve essere 
		 * "piccolo", o no?
		 */
		p->dSr = dfvx/(dfvax + p->dvThreshold);
		if (p->dSr > 1.) {
			p->dSr = 1.;
		}

		/*
		 * Direzione laterale: normale cross forward
		 */
		Vec3 lat = n.Cross(fwd);

		/*
		 * Velocita' laterale del mozzo
		 */
		doublereal dvay = lat.Dot(va);

		/*
		 * Angolo di deriva del mozzo
		 * Nota: ristretto al Io-IVo quadrante
		 * questa threshold sul modulo della velocita' fa in modo
		 * che l'angolo vada a zero se il modulo della velocita'
		 * e' troppo piccolo
		 */
		p->dAlpha = atan2(dvay, fabs(dvax));
		if (p->dAlphaThreshold > 0.) {
			doublereal dtmp = tanh(sqrt(dvax*dvax+dvay*dvay)/p->dAlphaThreshold);
			p->dAlpha *= dtmp*dtmp;
		}

		/*
		 * Coefficiente di attrito longitudinale
		 */
		doublereal dMuX0 = p->pMuX0->dGet(p->dSr);
		/*
		 * Nota: alpha/(pi/2) e' compreso tra -1. e 1.
		 */
		p->dMuX = dMuX0*sgn*(1.-fabs(p->dAlpha)/M_PI_2);
		
		/*
		 * Correggo le forze:
		 * uso il coefficiente di attrito longitudinale
		 * con il segno della velocita' del punto di contatto
		 */
		p->F -= fwd*dFn*p->dMuX;

		if (dvay != 0.) {
			doublereal dAlpha = p->dAlpha;

#if 0 /* non serve piu': Alpha e' ristretta al Io-IVo quadrante */
			if (dAlpha > M_PI_2) {
				dAlpha = M_PI-dAlpha;
			} else if (dAlpha < -M_PI_2) {
				dAlpha = -M_PI-dAlpha;
			}
#endif
			
			doublereal dMuY0 = p->pMuY0->dGet(dAlpha);
			doublereal dMuY1 = p->pMuY1->dGet(dAlpha);
			
			p->dMuY = dMuY0+(dMuY1-dMuY0)*p->dSr;

			/*
			 * Correggo le forze
			 */
			p->F -= lat*dFn*p->dMuY;
		}
	}

	/*
	 * Momento
	 */
	p->M = (pc-p->pWheel->GetXCurr()).Cross(p->F);

	WorkVec.Sub(1, p->F);
	WorkVec.Sub(4, (pc-p->pGround->GetXCurr()).Cross(p->F));
	WorkVec.Add(7, p->F);
	WorkVec.Add(10, p->M);

	return WorkVec;
}

static void
before_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

static void
after_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP)
{
	DEBUGCOUTFNAME("after_predict");
}

static void
update(LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP)
{
	DEBUGCOUTFNAME("update");
}

static unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_initial_num_dof");
	return 0;
}

static void
initial_work_space_dim(const LoadableElem* pEl, 
		integer* piNumRows, 
		integer* piNumCols)
{
	DEBUGCOUTFNAME("initial_work_space_dim");
	*piNumRows = 0;
	*piNumCols = 0;   
}

static VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_jac");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(iNumRows, iNumCols, 0.);
#if 0	
	module_wheel* p = (module_wheel *)pEl->pGetData();
#endif /* 0 */
   
	/* set sub-matrix indices and coefs */

	return WorkMat;
}

static SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_res");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.Resize(iNumRows);

#if 0
	module_wheel* p = (module_wheel *)pEl->pGetData(); 
#endif /* 0 */
	
	/* set sub-vector indices and coefs */
   
	return WorkVec;
}

static void
set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
	DEBUGCOUTFNAME("set_value");
}
   
static void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
	DEBUGCOUTFNAME("set_initial_value");
}

static unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_priv_data");
	return 0;
}

static doublereal
d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
	DEBUGCOUTFNAME("d_get_priv_data");
	ASSERT(pEl->iGetNumPrivData() > 0);
	if (i > pEl->iGetNumPrivData()) {
		std::cerr << "Module-template Elem: illegal private data index "
			<< i << std::endl;      
		THROW(ErrGeneric());
	}
	
	/* return i-th priv data */

	return 0.;
}

static void
destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("destroy");
	module_wheel* p = (module_wheel *)pEl->pGetData();
	
	SAFEDELETE(p);
}

static int
i_get_num_connected_nodes(const LoadableElem* /* pEl */ )
{
	DEBUGCOUTFNAME("i_get_num_connected_nodes");
	
	/* wheel + ground */
	return 2;
}

static void
get_connected_nodes(const LoadableElem* pEl, 
		int& NumNodes, Node::Type* NdTyp, unsigned int* NdLabels)
{
	DEBUGCOUTFNAME("get_connected_nodes");
	module_wheel* p = (module_wheel *)pEl->pGetData();
	
	/* wheel + ground */
	NumNodes = 2;

	NdTyp[0] = Node::STRUCTURAL;
	NdTyp[1] = Node::STRUCTURAL;

	NdLabels[0] = p->pWheel->GetLabel();
	NdLabels[1] = p->pGround->GetLabel();
}

static struct
LoadableCalls lc = {
	"wheel2",
	read, /* */
	NULL /* i_get_num_dof */ ,
	NULL /* set_dof */ ,
	output, /* */
	NULL /* restart */ ,
	work_space_dim, /* */
	ass_jac, /* */
	NULL /* ass_mats */ ,
	ass_res, /* */
	NULL /* before_predict */ ,
	NULL /* after_predict */ ,
	NULL /* update */ ,
	NULL /* after_convergence */ ,
	NULL /* i_get_initial_num_dof */ ,
	NULL /* initial_work_space_dim */ ,
	NULL /* initial_ass_jac */ ,
	NULL /* initial_ass_res */ ,
	NULL /* set_value */ ,
	NULL /* set_initial_value */ ,
	NULL /* i_get_num_priv_data */ ,
	NULL /* d_get_priv_data */ ,
	i_get_num_connected_nodes,
	get_connected_nodes,
	destroy
};

extern "C" {
void *calls = &lc;
}

