/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2004

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2004
Leonardo Cassan		<lcassan@tiscalinet.it>

*/

#ifndef VARIABLE_H
#define VARIABLE_H

#include <common.h>

struct s_variable: public s_card {

   double Ic;
   char* Expression;
   double Value;
   
   Boolean _Ic,_Expression,_Value;
   
   s_variable(void);
   ~s_variable(void);
   inline const char* const Gettype(void) const;
   Boolean Test(void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case IC: Set_Param (Ic,value,_Ic); break;
	 case FUNCTION: Set_Param (Value,value,_Value); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif

