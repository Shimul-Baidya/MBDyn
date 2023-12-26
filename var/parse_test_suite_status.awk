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

# parse_test_suite_status.awk: helper to parse the output from Octave's __run_test_suite__ function

BEGIN {
    failed = -1;
    passed = -1;
}

/^  PASS\>/ {
    passed = $2;
}

/^  FAIL\>/ {
    failed = $2;
}

END {
    if (failed < 0 || passed < 0) {
        printf("failed to parse \"%s\"\n", FILENAME);
        exit 1;
    }

    printf("%d/%d tests failed\n", failed, failed + passed);
    
    if (failed > 0) {
        printf("not all tests passed\n");
        exit 1;
    }

    if (passed == 0) {
        printf("no tests were executed\n");
        exit 1;        
    }

    printf("__run_test_suite__ passed\n");
}
