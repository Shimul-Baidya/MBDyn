/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

/* Classe che consente la manipolazione dei nomi di files
 * e delle estensioni.
 * Uso consigliato:
 * derivare una classe da questa, aggiungendo funzioni
 * che generino direttamente i nomi del file con l'estensione desiderata
 * Ad esempio:


 class MyFile : public FileName {
  public:
    MyFile(char* sFName) : FileName(sFName) { 
       NULL 
    };
    
    char* sOldFile(void) { 
       return _sPutExt(NULL); 
    };
    
    char* sDatFile(void) { 
       return _sPutExt(".dat"); 
    };
 };

*/

#ifndef FILENAME_H
#define FILENAME_H

extern "C" {
#include <stdlib.h>
#include <string.h>
}


const char EXT_SEP = '.';

#ifdef USE_DOS_FILE

const char DIR_SEP = '\\';

class FileName {
 protected:
   char* sName;
   char sExt[4];
   char* sRef;
   
 public:
   FileName(char* n = NULL);       // Acquisisce e seziona il nome del file
   virtual ~FileName(void);        // Dealloca le stringhe usate
   int iInit(char* n);             // Acquisisce e seziona il nome del file
   char* _sPutExt(char* n = NULL); // Aggiunge una nuova estensione (di default attacca la vecchia)
   char* sGet(void);               // Restituisce il nome del file con la vecchia estensione          
};

#endif /* USE_DOS_FILE */





#ifdef USE_UNIX_FILE

const char DIR_SEP = '/';

class FileName {
 private:
   char* sName;
   char* sExt;
   char* sRef;
   unsigned int iMaxSize;
   unsigned int iCurSize;
   
 public:
   FileName(const char* n = NULL, int i = 0);
   virtual ~FileName(void);
   int iInit(const char* n, int i = 0);
   char* _sPutExt(const char* n);
   char* sGet(void) const;
};

#endif /* USE_UNIX_FILE */

#endif /* FILENAME_H */

