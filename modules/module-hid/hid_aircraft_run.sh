#!/bin/sh -f

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
# 

echo "warning: hid_aircraft requires user interaction and cannot be executed by MBDyn's testsuite!"

input_file=hid_aircraft
output_file=hid_aircraft

while ! test -z "$1"; do
    case "$1" in
        -f)
            input_file="$2"
            shift
            ;;
        -o)
            output_file="$2"
            shift
            ;;
    esac
    shift
done

rm -f aircraft.sock
exec mbdyn -f "${input_file}" -o "${output_file}"
