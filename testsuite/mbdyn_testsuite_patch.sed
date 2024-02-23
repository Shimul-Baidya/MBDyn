# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2023
#
# Pierangelo Masarati	<pierangelo.masarati@polimi.it>
# Paolo Mantegazza	<paolo.mantegazza@polimi.it>
#
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
#
# Changing this copyright notice is forbidden.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
#
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
# Copyright (C) 2023(-2023) all rights reserved.

# The copyright of this code is transferred
# to Pierangelo Masarati and Paolo Mantegazza
# for use in the software MBDyn as described
# in the GNU Public License version 2.1

# Use mbdyn_testsuite_patch.sed to patch input files (e.g. in order to test different solver settings)

## Add a new line after "begin: initial value;"
/^[[:space:]]*\<begin\>:[[:space:]]*\<initial[[:space:]]*value\>[[:space:]]*;[[:space:]]*$/a include: "${MBD_TESTSUITE_INITIAL_VALUE_BEGIN}";

## Add a new line before "end: initial value;"
/^[[:space:]]*\<end\>:[[:space:]]*\<initial[[:space:]]*value\>[[:space:]]*;[[:space:]]*$/i include: "${MBD_TESTSUITE_INITIAL_VALUE_END}";

## Add a new line after "begin: control data;"
/^[[:space:]]*\<begin\>:[[:space:]]*control[[:space:]]*data[[:space:]]*;[[:space:]]*$/a include: "${MBD_TESTSUITE_CONTROL_DATA_BEGIN}";

## Add a new line before "end: control data;"
/^[[:space:]]*\<end\>:[[:space:]]*control[[:space:]]*data[[:space:]]*;[[:space:]]*$/i include: "${MBD_TESTSUITE_CONTROL_DATA_END}";

## Remove "method: <method_data>;"
s/^[[:space:]]*\<method\>[[:space:]]*:.*;.*//g

## Remove "tolerance: <tolerance_data>;"
s/^[[:space:]]*\<tolerance\>[[:space:]]*:.*;.*//g

## Remove "derivatives tolerance: <tolerance_data>;"
s/^[[:space:]]*\<derivatives[[:space:]]*tolerance\>[[:space:]]*:.*;.*//g

## Remove "linear solver: <solver_data>;"
s/^[[:space:]]*(\<linear|\<nonlinear)[[:space:]]*solver[[:space:]]*:([[:space:]]|[[:alnum:]]|[,]|[\.]|[+]|[-])*;[[:space:]]*$//g

## Remove "threads: assembly;" and "threads: solver;"
s/^[[:space:]]*threads[[:space:]]*:[[:space:]]*(assembly|solver)[[:space:]]*,([[:space:]]|[[:alnum:]]|[,]|[\.]|[\n])*;[[:space:]]*$//g

## Remove "threads: disable;"
s/^[[:space:]]*threads[[:space:]]*:[[:space:]]*disable[[:space:]]*;[[:space:]]*$//g

## Remove "use automatic differentiation;"
s/^[[:space:]]*use[[:space:]]*automatic[[:space:]]*differentiation[[:space:]]*;[[:space:]]$//g
