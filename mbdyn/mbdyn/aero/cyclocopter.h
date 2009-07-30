/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

/* Cyclocopter inflow models */

#ifndef CYCLOCOPTER_H
#define CYCLOCOPTER_H

#include "indvel.h"

/* CyclocopterNoInflow - begin */

class CyclocopterNoInflow
: virtual public Elem, public InducedVelocity {
protected:
	const StructNode* pRotor;

	Mat3x3 RRot, RRotTranspose;

public:
	CyclocopterNoInflow(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		flag fOut);
	virtual ~CyclocopterNoInflow(void);

	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(unsigned int uL, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(const Vec3& X) const;

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

/* CyclocopterNoInflow - end */

/* CyclocopterUniform1D - begin */

/*

From Moble: the induced velocity is opposite to the
the force generated by the rotor in the direction 3
of the rotor reference (the direction 1 of this 
reference must be align with the rotation axis of 
the rotor). It should make sense just for a multiblade
rotor (not for the one-blade rotor) when the generated 
force is mainly in one direction.

*/

class CyclocopterUniform1D
: virtual public Elem, public InducedVelocity {
protected:
	const StructNode* pRotor;

	doublereal dOmegaRef;		// Reference rotation speed
	doublereal dRadius;		// Rotor radius
	doublereal dSpan;		// Blade length
	doublereal dArea;		// Cylinder longitudinal area
	
	DriveOwner Weight;
	doublereal dWeight;

	Mat3x3 RRot, RRotTranspose;
	Vec3 RRot3;
	doublereal dUind;
	doublereal dUindPrev;

public:
	CyclocopterUniform1D(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		const doublereal& dOR, const doublereal& dR,
		const doublereal& dL, DriveCaller *pdW,
		flag fOut);
	virtual ~CyclocopterUniform1D(void);

	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(unsigned int uL, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(const Vec3& X) const;

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

/* CyclocopterUnifor1D - end */

/* CyclocopterUniform2D - begin */

/*

Uniform inflow: the indiced velocity is opposite to
the force generated by the rotor in the plane 
perdendicular to the rotor rotation axis. The
rotor reference must have the direction 1 aligned 
with the rotor rotation axis!

*/

class CyclocopterUniform2D
: virtual public Elem, public InducedVelocity {
protected:
	const StructNode* pRotor;

	doublereal dOmegaRef;		// Reference rotation speed
	doublereal dRadius;		// Rotor radius
	doublereal dSpan;		// Blade length
	doublereal dArea;		// Cylinder longitudinal area
	
	DriveOwner Weight;
	doublereal dWeight;

	Mat3x3 RRot;
	Mat3x3 RRotor;
	Vec3 dUind, dUindPrev;
	doublereal dUindMagnitude;

public:
	CyclocopterUniform2D(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		const doublereal& dOR, const doublereal& dR,
		const doublereal& dL, DriveCaller *pdW,
		flag fOut);
	virtual ~CyclocopterUniform2D(void);

	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(unsigned int uL, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(const Vec3& X) const;

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

/* CyclocopterUnifor2D - end */

/* CyclocopterKARI - begin */

/*

From:

"A New VTOL UAV Cyclocopter with Cycloidal Blades System",

Chul Yong Yun, Illkyung Park,
Ho Yong Lee, Jai Sang Jung, In Seong Hwang,
Seung Jo Kim,
Sung Nam Jung

Presented at the American Helicopter Society 60th Annual Forum,
Baltimore, MD, June 7-10, 2004.
*/

class CyclocopterKARI
: virtual public Elem, public InducedVelocity {
protected:
	const StructNode* pRotor;
	Mat3x3 RRot;

public:
	CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, flag fOut);
	virtual ~CyclocopterKARI(void);

	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(unsigned int uL, const Vec3& F, const Vec3& M, const Vec3& X);
#endif

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(const Vec3& X) const;

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

/* CyclocopterKARI - end */

class DataManager;
class MBDynParser;

extern Elem*
ReadCyclocopter(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO, 
	unsigned int uLabel,
	const StructNode* pC,
	const Mat3x3& rrot,
	const StructNode* pR);

#endif /* CYCLOCOPTER_H */

