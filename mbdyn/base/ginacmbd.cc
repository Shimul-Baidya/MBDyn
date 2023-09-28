/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Paolo Mantegazza     <paolo.mantegazza@polimi.it>
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

#include "ginacmbd.h"

#ifdef USE_MULTITHREAD
int GiNaCEntity::iInitMutex;
pthread_mutex_t GiNaCEntity::GiNaCMutex;
#endif

GiNaCEntity::GiNaCEntity()
{
#ifdef USE_MULTITHREAD
     ASSERT(iInitMutex >= 0);

     if (!iInitMutex) {
          if (pthread_mutex_init(&GiNaCMutex, nullptr)) {
               silent_cerr("GiNaCEntity::GiNaCEntity():"
                           "mutex init failed\n");
               throw ErrGeneric(MBDYN_EXCEPT_ARGS);
          }
     }

     ++iInitMutex;
#endif
}

GiNaCEntity::~GiNaCEntity()
{
#ifdef USE_MULTITHREAD
     ASSERT(iInitMutex > 0);

     --iInitMutex;

     if (!iInitMutex) {
          if (pthread_mutex_destroy(&GiNaCMutex)) {
               silent_cerr("GiNaCEntity::~GiNaCEntity(): "
                           "mutex init failed\n");
               ASSERT(0);
          }
     }
#endif
}
