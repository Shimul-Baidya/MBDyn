/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <bulk.h>
#include <dataman.h>

/* Legge un elemento bulk */
   
Elem* ReadBulk(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadBulk");
   
   const char* sKeyWords[] = {
      "springsupport"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	SPRINGSUPPORT = 0,	
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
   
   /* lettura del tipo di elemento elettrico */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {      
      cout << "bulk element type: " 
	<< sKeyWords[CurrKeyWord] << endl;
   }   
#endif   

   Elem* pEl = NULL;
   
   switch (CurrKeyWord) {
      /*  */
    case SPRINGSUPPORT: {       
       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       if (SD.pNode->GetNodeType() ==  Node::PARAMETER) {
	  cerr << "Sorry, parameters are not allowed for bulk spring" << endl;
	  THROW(DataManager::ErrGeneric());	      
       }	     	  
       
       DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
       flag fOut = pDM->fReadOutput(HP, Elem::BULK);
       
       SAFENEWWITHCONSTRUCTOR(pEl, 
			      BulkSpringSupport,
			      BulkSpringSupport(uLabel, pDC, SD, fOut));
       
       break;
    }
      
      /* Aggiungere altri elementi elettrici */
      
    default: {
       cerr << endl
	 << "unknown bulk element type in bulk element " << uLabel
	 << " at line " << HP.GetLineData() << endl;       
       THROW(DataManager::ErrGeneric());
    }	
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      cerr << endl
	<< "semicolon expected at line " << HP.GetLineData() << endl;     
      THROW(DataManager::ErrGeneric());
   }   
   
   return pEl;
} /* End of ReadBulk() */
