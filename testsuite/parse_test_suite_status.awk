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
    ## Output from Octave's function "__run_test_suite__"
    passed = $2;
}

/^  FAIL\>/ {
    ## Output from Octave's function "__run_test_suite__"
    failed = $2;
}

/^PASSES\>/ {
    ## Output from Octave's function "test"
    passed = strtonum($2);
    total = strtonum($5);
    failed = total - passed;
}

/^FAILED\>/ {
    ## Output from Octave's function "test"
    failed = strtonum($2);
}

/\?\?\?\?\? .+has no tests available$|\?\?\?\?\? .+source code with tests for dynamically linked function not found$/ {
    ## Output from Octave's function "test"
    ## This is not considered as a failure.
    passed = 0;
    failed = 0;
}

/^!!!!! test failed$/ {
    if (failed < 0) {
        failed = 0;
        passed = 0;
    }

    ++failed;
}

END {
    if (failed < 0 || passed < 0) {
        printf("Failed to parse output file \"%s\"!\n", FILENAME);
        exit 1;
    }

    printf("%d/%d tests passed!\n", passed, failed + passed);

    if (failed > 0) {
        printf("%d/%d tests failed!\n", failed, failed + passed);
        exit 1;
    }

    if (passed == 0) {
        printf("No tests were executed\n");
        exit 0;
    }
}
