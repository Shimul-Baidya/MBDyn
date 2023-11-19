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

/* 
 * Copyright 1999-2000 Lamberto Puggelli <puggelli@tiscalinet.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano 
 */

#ifndef PRESNODE_H
#define PRESNODE_H

#include "scalarnode.h"

class PressureNode : virtual public ScalarAlgebraicNode {
public:
     PressureNode(unsigned int uL, const DofOwner* pDO, doublereal dx, DofOrder::Equality eEqualityType, flag fOut); 
   
     virtual ~PressureNode();
   
     virtual Node::Type GetNodeType(void) const;
   
	virtual OutputHandler::OutFiles GetOutputType(void) const { return OutputHandler::PRESNODES; };
	virtual void OutputPrepare(OutputHandler &OH);

     /* returns the dimension of the component */
     const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;

     /* describes the dimension of components of equation */
     virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix = "", bool bInitial = false) const;
};

#endif /* PRESNODE_H */

