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

#ifndef ___GINAC_H__INCLUDED___
#define ___GINAC_H__INCLUDED___

#include "myassert.h"

#ifdef USE_MULTITHREAD
#include <ac/pthread.h>
#endif

class GiNaCEntity {
protected:
     GiNaCEntity();
     ~GiNaCEntity();

#ifdef USE_MULTITHREAD
     friend class GiNaCGuard;
     class GiNaCGuard {
     public:
          GiNaCGuard() {
               ASSERT(iInitMutex > 0);
               pthread_mutex_lock(&GiNaCEntity::GiNaCMutex);
          }
          ~GiNaCGuard() {
               ASSERT(iInitMutex > 0);
               pthread_mutex_unlock(&GiNaCEntity::GiNaCMutex);
          }
     };
#endif
private:
#ifdef USE_MULTITHREAD
     static int iInitMutex;
     static pthread_mutex_t GiNaCMutex;
#endif
};

#endif
