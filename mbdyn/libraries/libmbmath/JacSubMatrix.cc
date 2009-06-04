/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003
 *
 * Marco Morandini	<morandini@aero.polimi.it>
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

/* here goes Morandini's copyright */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <algorithm>

#include "myassert.h"
#include "JacSubMatrix.h"

ExpandableRowVector::ExpandableRowVector(void) {};
ExpandableRowVector::ExpandableRowVector(const integer n) {
	ReDim(n);
}
ExpandableRowVector::~ExpandableRowVector(void) {};
void ExpandableRowVector::ReDim(const integer n) {
	//we have to accept = 0, some elements do ReDim(0,0) (PointForceElement))
	ASSERTMSGBREAK(n>=0, "ExpandableRowVector:ReDim(), n shold be >= 0");
	x.resize(n,0.);
	xm.resize(n,0);
	idx.resize(n,0); 
}
void ExpandableRowVector::Zero(void) {
	std::fill(x.begin(),x.end(),0.);
}
void ExpandableRowVector::Reset(void) {
	Zero();
	std::fill(xm.begin(),xm.end(),(ExpandableRowVector*)0);
	std::fill(idx.begin(),idx.end(),0);		
}
void ExpandableRowVector::Link(const integer i, const ExpandableRowVector*const xp) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Link() underflow");
	ASSERTMSGBREAK(i <= idx.size(), "ExpandableRowVector::Link() overflow");
	ASSERTMSGBREAK(idx[i-1] == 0, "ExpandableRowVector::Link() fatal error");
	xm[i-1] = xp;
}
void ExpandableRowVector::Set(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Set() underflow");
	ASSERTMSGBREAK(i <= x.size(), "ExpandableRowVector::Set() overflow");
	x[i-1] = xx;
}
void ExpandableRowVector::SetIdx(integer i, integer iidx) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::SetIdx() underflow");
	ASSERTMSGBREAK(i <= idx.size(), "ExpandableRowVector::SetIdx() overflow");
	idx[i-1] = iidx;
}
void ExpandableRowVector::Set(doublereal xx, integer i, integer iidx) {
	Set(xx,i);
	SetIdx(i,iidx);
}
doublereal&
ExpandableRowVector::operator ()(integer i)
{
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
	ASSERTMSGBREAK(i <= x.size(), "ExpandableRowVector::() overflow");
	return x[i-1];
}
const doublereal&
ExpandableRowVector::operator ()(integer i) const
{
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::() underflow");
	ASSERTMSGBREAK(i <= x.size(), "ExpandableRowVector::() overflow");
	return x[i-1];
}
void ExpandableRowVector::Add(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Add() underflow");
	ASSERTMSGBREAK(i <= x.size(), "ExpandableRowVector::Add() overflow");
	x[i-1] += xx;
}
void ExpandableRowVector::Sub(doublereal xx, integer i) {
	ASSERTMSGBREAK(i > 0, "ExpandableRowVector::Sub() underflow");
	ASSERTMSGBREAK(i <= x.size(), "ExpandableRowVector::Sub() overflow");
	x[i-1] -= xx;
}
void ExpandableRowVector::Add(SubVectorHandler& WorkVec, const doublereal c) const {
	for (std::vector<doublereal>::size_type i=0; i<x.size(); i++) {
		if (x[i] == 0.) {
			continue;
		}

		if (idx[i] != 0) {
			WorkVec.Add(idx[i],c*x[i]);
		} else {
			ASSERTMSGBREAK(xm[i] != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
			xm[i]->Add(WorkVec,c*x[i]);
		}
	}
}
void ExpandableRowVector::Sub(SubVectorHandler& WorkVec, const doublereal c) const {
	for (std::vector<doublereal>::size_type i=0; i<x.size(); i++) {
		if (x[i] == 0.) {
			continue;
		}

		if (idx[i] != 0) {
			WorkVec.Sub(idx[i],c*x[i]);
		} else {
			ASSERTMSGBREAK(xm[i] != 0, "ExpandableRowVector::Sub() null pointer to ExpandableRowVector");
			xm[i]->Sub(WorkVec,c*x[i]);
		}
	}
}
void ExpandableRowVector::Add(FullSubMatrixHandler& WM, 
	const integer eq,
	const doublereal c) const {
	for (std::vector<doublereal>::size_type i=0; i<x.size(); i++) {
		if (x[i] == 0.) {
			continue;
		}

		if (idx[i] != 0) {
			WM.IncCoef(eq,idx[i],c*x[i]);
		} else {
			ASSERTMSGBREAK(xm[i] != 0, "ExpandableRowVector::Add() null pointer to ExpandableRowVector");
			xm[i]->Add(WM,eq,c*x[i]);
		}
	}
}
void ExpandableRowVector::Sub(FullSubMatrixHandler& WM,
	const integer eq,
	const doublereal c) const {
	for (std::vector<doublereal>::size_type i=0; i<x.size(); i++) {
		if (x[i] == 0.) {
			continue;
		}

		if (idx[i] != 0) {
			WM.DecCoef(eq,idx[i],c*x[i]);
		} else {
			ASSERTMSGBREAK(xm[i] != 0, "ExpandableRowVector::Sub() null pointer to ExpandableRowVector");
			xm[i]->Sub(WM,eq,c*x[i]);
		}
	}
}

#include<iomanip>
std::ostream & ExpandableRowVector::Write(std::ostream &out, const char *sFill) const {
	out << "LocalDof: ";
	for (std::vector<doublereal>::size_type i=0; i<x.size(); i++) {
		if (idx[i] != 0) {
			out << sFill << std::setw(12) << idx[i];
		} else {
			out << sFill << std::setw(12) << "linked";
		}		
	}
	out << std::endl;
	out << "   Value: ";
	for (std::vector<doublereal>::size_type i=0; i<x.size(); i++) {
		out << sFill << std::setw(12) << x[i];
	}
	out << std::endl;
	return out;
}

std::ostream & operator << (std::ostream & s, const ExpandableRowVector & z) {
	return z.Write(s);	
}

