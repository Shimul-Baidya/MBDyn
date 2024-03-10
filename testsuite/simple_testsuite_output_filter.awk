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

# AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
#        Copyright (C) 2023(-2023) all rights reserved.
#
#        The copyright of this code is transferred
#        to Pierangelo Masarati and Paolo Mantegazza
#        for use in the software MBDyn as described
#        in the GNU Public License version 2.1

# purpose:
## Currently the amount of output is limited to 4MB in GitLab-CI.
## Some of our tests based on simple_testsuite.sh are exceeding this limit.
## In such a situation you may get a message like this one:

## Job's log exceeded limit of 4194304 bytes.
## Job execution will continue but no more output will be collected.

## As a workaround, we may run simple_testsuite.sh as follows:
## simple_testsuite.sh ... | awk -f simple_testsuite_output_filter.awk

BEGINFILE {
    output = 0;
    reports = 0;
}

/^@END_SIMPLE_TESTSUITE_REPORT@$/ {
    if (output) {
        ++reports;
    }
    output = 0;
}

output != 0 {
    print
}

/^@BEGIN_SIMPLE_TESTSUITE_REPORT@$/ {
    output = 1;
}

ENDFILE {
    printf("%d reports found\n", reports);
}
