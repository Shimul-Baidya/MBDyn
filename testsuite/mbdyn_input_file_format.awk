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

## Use mbdyn_input_file_format.awk in order to check if a regular file is an input file for MBDyn.

BEGIN {
    MBD_SEC_NONE       = 0x00;
    MBD_SEC_DATA       = 0x01;
    MBD_SEC_INIT_VAL   = 0x02;
    MBD_SEC_INV_DYN    = 0x04;
    MBD_SEC_CONTR_DATA = 0x08;
    MBD_SEC_NODES      = 0x10;
    MBD_SEC_ELEM       = 0x20;
}

BEGINFILE {
    MBD_SEC_ALL_INIT_VAL = or(MBD_SEC_DATA, MBD_SEC_INIT_VAL, MBD_SEC_CONTR_DATA, MBD_SEC_NODES, MBD_SEC_ELEM);
    MBD_SEC_ALL_INV_DYN = or(MBD_SEC_DATA, MBD_SEC_INV_DYN, MBD_SEC_CONTR_DATA, MBD_SEC_NODES, MBD_SEC_ELEM);
    current_section = 0;
    sections_found = 0;
}

/^[[:space:]]*begin:[[:space:]]*data[[:space:]]*;[[:space:]]*$/ {
    # print "begin data"

    current_section = MBD_SEC_DATA;
}

/^[[:space:]]*end:[[:space:]]*data[[:space:]]*;[[:space:]]*$/ {
    # print "end data"

    if (current_section == MBD_SEC_DATA) {
        sections_found = or(sections_found, current_section);
    }

    current_section = MBD_SEC_NONE;
}

/^[[:space:]]*begin:[[:space:]]*initial[[:space:]]*value[[:space:]]*;[[:space:]]*$/ {
    # print "begin initial value"

    current_section = MBD_SEC_INIT_VAL;
}

/^[[:space:]]*end:[[:space:]]*initial[[:space:]]*value[[:space:]]*;[[:space:]]*$/ {
    # print "end initial value"

    if (current_section == MBD_SEC_INIT_VAL) {
        sections_found = or(sections_found, current_section);
    }

    current_section = MBD_SEC_NONE;
}

/^[[:space:]]*begin:[[:space:]]*inverse[[:space:]]*dynamics[[:space:]]*;[[:space:]]*$/ {
    # print "begin inverse dynamics"

    current_section = MBD_SEC_INV_DYN;
}

/^[[:space:]]*end:[[:space:]]*inverse[[:space:]]*dynamics[[:space:]]*;[[:space:]]*$/ {
    # print "end inverse dynamics"

    if (current_section == MBD_SEC_INV_DYN) {
        sections_found = or(sections_found, current_section);
    }

    current_section = MBD_SEC_NONE;
}

/^[[:space:]]*begin:[[:space:]]*control[[:space:]]*data[[:space:]]*;[[:space:]]*$/ {
    # print "begin control data"

    current_section = MBD_SEC_CONTR_DATA;
}

/^[[:space:]]*end:[[:space:]]*control[[:space:]]*data[[:space:]]*;[[:space:]]*$/ {
    # print "end control data"

    if (current_section == MBD_SEC_CONTR_DATA) {
        sections_found = or(sections_found, current_section);
    }

    current_section = MBD_SEC_NONE;
}

/^[[:space:]]*begin:[[:space:]]*nodes[[:space:]]*;[[:space:]]*$/ {
    # print "begin nodes"

    current_section = MBD_SEC_NODES;
}

/^[[:space:]]*end:[[:space:]]*nodes[[:space:]]*;[[:space:]]*$/ {
    # print "end nodes"

    if (current_section == MBD_SEC_NODES) {
        sections_found = or(sections_found, current_section);
    }

    current_section = MBD_SEC_NONE;
}

/^[[:space:]]*begin:[[:space:]]*elements[[:space:]]*;[[:space:]]*$/ {
    # print "begin elements"

    current_section = MBD_SEC_ELEM;
}

/^[[:space:]]*end:[[:space:]]*elements[[:space:]]*;[[:space:]]*$/ {
    # print "end elements"

    if (current_section == MBD_SEC_ELEM) {
        sections_found = or(sections_found, current_section);
    }

    current_section = MBD_SEC_NONE;
}

ENDFILE {
    if ((!exclude_initial_value && sections_found == MBD_SEC_ALL_INIT_VAL) || (!exclude_inverse_dynamics && sections_found == MBD_SEC_ALL_INV_DYN)) {
        print FILENAME;
    }
}
