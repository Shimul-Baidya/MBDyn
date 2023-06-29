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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "mynewmem.h"
#include "node.h"
#include "nodead.h"
#include "solman.h"


/* Node - begin */

/* Costruttore */
Node::Node(unsigned int uL, const DofOwner* pDO, flag fOut)
: WithLabel(uL), DofOwnerOwner(pDO), ToBeOutput(fOut)
{
	NO_OP;
}

/* Distruttore banale */
Node::~Node(void)
{
	NO_OP;
}

std::ostream&
Node::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

void
Node::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i <= 0);
	desc.resize(0);
}

std::ostream&
Node::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

void
Node::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i <= 0);
	desc.resize(0);
}

/* Ritorna gli indici di riga e colonna. Tipicamente sono gli stessi */
integer
Node::iGetFirstRowIndex(void) const
{
	return iGetFirstIndex();
}

integer
Node::iGetFirstColIndex(void) const
{
	return iGetFirstIndex();
}

Node::Type
str2nodetype(const char *const s)
{
	for (int i = 0; i < Node::LASTNODETYPE; i++) {
		if (strcasecmp(s, psReadNodesNodes[i]) == 0) {
			return Node::Type(i);
		}
	}

	return Node::UNKNOWN;
}

void Node::UpdateJac(doublereal dCoef)
{
}

void Node::UpdateJac(const VectorHandler& Y, doublereal dCoef)
{
}

/* Node - end */

