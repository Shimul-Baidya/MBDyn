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

/* nodi */

#ifndef NODE_H
#define NODE_H

#include "myassert.h"

#include "output.h"
#include "withlab.h"
#include "dofown.h"
#include "simentity.h"

/*
 * Array dei nomi dei nodi.
 * Usato per output
 * @see Node::Type
 */
extern const char* psNodeNames[];

/*
 * Array delle stringhe di identificazione dei tipi di nodi.
 * Usato per input di controllo
 * @see Node::Type
 */
extern const char* psReadControlNodes[];

/*
 * Array delle stringhe di identificazione dei tipi di nodi.
 * Usato per input dei nodi
 * @see Node::Type
 */
extern const char* psReadNodesNodes[];

/* Node - begin */ 

class Node : public WithLabel, public SimulationEntity,
public DofOwnerOwner, public ToBeOutput
{
public:
	/* Enumerazione dei tipi di nodi */
	enum Type {
		UNKNOWN = -1,

		/* Should be Node::SCALAR; keep using Node::ABSTRACT
		 * for backward compatibility */
		ABSTRACT = 0,

		STRUCTURAL,
		ELECTRIC,
		THERMAL,
		PARAMETER,
		HYDRAULIC,

		LASTNODETYPE
	};

public:
	/* Errori: */
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
 	};

	/* Costruttori */

	/*
	 * Costruttore.
	 * @param uL label
	 * @param pDO puntatore al DofOwner relativo, gestito da DofOwnerOwner
	 * @param fOut flag di output, gestito da ToBeOutput
	 * @see DofOwner
	 * @see DofOwnerOwner
	 * @see ToBeOutput
	 */
	Node(unsigned int uL, const DofOwner* pDO, flag fOut);

	/* Distruttore */
	virtual ~Node(void);

	/* Funzioni di servizio */
	const Node *GetNode(void) const { return this; };

	/* Tipo del nodo (usato per debug ecc.) */
	virtual Node::Type GetNodeType(void) const = 0;

	/* Contributo del nodo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const = 0;

	/*
	 * Describe the degrees of freedom
	 */
	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	/* Metodi che operano sui DoF */

	/*
	 * Ritorna il primo indice di riga dei DoF del nodo, in base 0.
	 * Ovvero, l'indice del primo DoF del nodo in un vettore a base zero.
	 * Per avere gli indici in un vettore a base 1 (stile Fortran),
	 * occorre sommare al risultato il numero del DoF.
	 * Vedi iGetFirstColIndex()
	 */
	virtual integer iGetFirstRowIndex(void) const;

	/*
	 * Ritorna gli indici di colonna dei DoF.
	 * Per la numerazione degli indici vedi iGetFirstRowIndex().
	 * Tipicamente gli indici di riga e di colonna sono gli stessi,
	 * tranne in alcuni casi notevoli.
	 */
	virtual integer iGetFirstColIndex(void) const;

	/*
	 * Restituisce il valore del DoF iDof.
	 * Se il nodo e' differenziale, iOrder puo' essere = 1 per avere la derivata
	 */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const = 0;

	/*
	 * Restituisce il valore del DoF iDof al passo precedente.
	 * Se il nodo e' differenziale, iOrder puo' essere = 1 per avere la derivata
	 */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const = 0;

	/*
	 * Setta il valore del DoF iDof a dValue.
	 * Se il nodo e' differenziale, iOrder puo' essere = 1
	 * per operare sulla derivata
	 */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof,
		unsigned int iOrder = 0) = 0;

	/*
	 * priv data and dofs coincide for nodes, unless overridden */

	/*
	 * Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili
	 */
	virtual unsigned int iGetNumPrivData(void) const {
		return iGetNumDof();
	};

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return 0;
	};

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const {
		return dGetDofValue(i);
	};

        /* 
         * Automatic differentiation support:
         * Called on each node before assembly of the Jacobian 
         */
        virtual void UpdateJac(doublereal dCoef);
        /* 
         * Automatic differentiation support:
         * Called on each node before assembly of the Jacobian vector product Jac * Y 
         */
        virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef);
};

Node::Type str2nodetype(const char *const s);

/* Node - end */

#endif /* NODE_H */

