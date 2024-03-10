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

# Purpose:
# Perform a classification of error messages.
# So, syntax errors and socket errors may be suppressed for mbdyn-tests-private-test-job.yml.
# For example:
# awk -v suppressed_errors="syntax|module|loadable|feature|socket" parse_mbdyn_error_message.awk ${log_file}
# if test $? -eq 0; then
#   echo "suppressed"
# else
#   echo "failed"
# fi

function set_error(err, msg)
{
    global::error_found = err;
    global::message = msg;
}

BEGIN {
    status = 0;
}

BEGINFILE {
    set_error("unknown", "");
}

/^\[.*\/libraries\/libmbutil\/mathp\.cc:[0-9]*,func=ExpressionElement\* MathParser::stmt\(\)\] \(unable to find variable ".*"\) at line [0-9]*$/ {
    set_error("syntax", "variable");
}

/^HighParser::GetValue\(\): error return from MathParser at line [0-9]+, file <.*>$/ {
    set_error("syntax", "getvalue");
}

/^Parser error in HighParser::GetDescription, invalid call to GetDescription at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "getdescription");
}

/^Cannot stat file <.+> at line [[:digit:]]+, file <.+>/ {
    set_error("syntax", "filenotfound");
}

/^hydraulic fluid expected at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "hydraulicfluid");
}

/^unknown hydraulic element type for hydraulic element [[:digit:]]+ at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "hydraulicfluid");
}

/^need at least [[:digit:]]+ Gauss integration point at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "aerodynamic");
}

/^C81Data\([[:digit:]]+\): unable to open file '.+' at line [[:digit:]]+, file <.+>$/ {
    if (split($0, files, /[\<\>']/) > 2) {
        msg = files[2];
    } else {
        msg = "";
    }

    set_error("syntax", sprintf("c81data:%s", msg));
}

/^FixedStepFileDrive\([[:digit:]]+\): can't open file ".+"$/ {
    set_error("syntax", "fixedstepfiledrive");
}

/^Unknown integration method at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "method");
}

/^line [[:digit:]]+, file <.+>: illegal viscosity/ {
    set_error("syntax", "hydraulicfluid");
}

/^semicolon expected .+ line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "semicolon");
}

/^Force\([[:digit:]]+\): "position" keyword expected at line [[:digit:]]+, file <.+>; still using deprecated syntax\?$/ {
    set_error("syntax", "force");
}

/^DataManager::ReadNode: Structural\([[:digit:]]+\) not defined at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "undefined-node");
}

/^ReadStructNode\([[:digit:]]+\): semicolon expected at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "structnode");
}

/^AerodynamicBody\([[:digit:]]+\): unknown output mode at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "aerodynamicbody");
}

/^"output file name" no longer supported at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "outputfilename");
}

/^\[.+libraries\/libmbutil\/mathp.cc:[[:digit:]]+,func=ExpressionElement\* MathParser::GetExpr\(\)] \(statement separator expected) at line [[:digit:]]+$/ {
    set_error("syntax", "separator");
}

/^StringDriveCaller::dGet\(\): \[.+libraries\/libmbutil\/mathp.cc:[[:digit:]]+,func=int mp_stop\(const MathArgs&\)\]$/ {
    set_error("syntax", "stringdrivecaller");
}

/^SocketStreamDrive\([[:digit:]]+, ".+"\): unable to read local path at line [[:digit:]]+, file <.+>$/ {
    set_error("socket", "path");
}

/^ExtSocketHandler\([[:digit:]]+\): unable to read local path at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "extsockedhandler");
}

/^ReadAuthMethod: line [[:digit:]]+, file <.+>: no working crypt\([[:digit:]]+\)$/ {
    set_error("syntax", "crypt");
}

/^error - unknown element type at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "element");
}

/^template drive caller [[:digit:]]+ is undefined at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "drivecaller");
}

/^Parser error in HighParser::IsKeyWord\(\), missing separator at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "separator");
}

/^Unknown description at line [[:digit:]]+, file <.+>; aborting\.\.\.$/ ||
/^unknown description at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "description");
}

/^Error while expecting "end: elements;" at line [[:digit:]]+, file <.+>$/ {
    if (!match(global::error_found, /module/)) {
        set_error("syntax", "elements");
    }
}

/^PlaneDispJoint\([[:digit:]]\): unsupported; use an InPlane and a RevoluteRotation$/ {
    set_error("syntax", "inplanejoint");
}

/^ExtForce\([[:digit:]]+\): unrecognised communicator type at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "extforce");
}

/^Joint\([[:digit:]]+\): "kinematic" obolete; replace with a "total \[pin\] joint" at line [[:digit:]]+, file <.+>$/ {
    set_error("syntax", "joint");
}

/^error in allocation of Joint Regularization\([[:digit:]]+\)$/ {
    set_error("element", "jointregularization");
}

/^RTPOSIXSolver: sched_setscheduler failed \(1: .*\)$/ {
    set_error("feature", "RTPOSIX");
}

/^ReadRTSolver: need to configure --with-rtai to use default RTAI real-time solver at line [[:digit:]]+, file <.+>$/ {
    set_error("feature", "RTAI");
}

/^ReadAuthMethod: line [[:digit:]]+, file <.+>: no PAM support$/ {
    set_error("feature", "PAM");
}

/^solid\([[:digit:]]+\): initial assembly with incompressible constitutive laws is not yet implemented!$/ {
    set_error("feature", "solid-upc-assembly");
}

/^ModuleLoad_int: unable to open module\>/ {

    if (split($0, modules, /[\<\>]/) > 2) {
        msg = modules[2];
    } else {
        msg = "";
    }

    set_error("module", msg);
}

/^ParseUserDefinedElem\([0-9]+\): unknown user-defined element type at line [0-9]+, file <.+>$/ {
    if (split($0, elements, /[\<\>]/) > 2) {
        msg = elements[2];
    } else {
        msg = "";
    }

    set_error("loadable", msg);
}

/^Loadable\([[:digit:]]+\): unable to open module <.+> \(file not found\) at line [[:digit:]]+, file <.+>$/ {
    if (split($0, modules, /[\<\>]/) > 2) {
        msg = modules[2];
    } else {
        msg = "";
    }

    set_error("loadable", msg);
}

/^(Harwell|Meschach) solver not available; requested at line [[:digit:]]+, file <.+>$/ {
    set_error("feature", "linearsolver");
}

/^"use jdqz" needs to configure --with-jdqz at line [[:digit:]]+, file <.+>$/ {
    set_error("feature", "jdqz");
}

/^UseLocalSocket\(".+"\): bind\(\) failed/ {
    set_error("socket", "bind");
}

/^UseSocket\(\): connection timed out$/ {
    set_error("socket", "timeout");
}

/^SocketDrive\([[:digit:]]+\): bind failed \([[:digit:]]+: Address already in use\)[[:space:]]*$/ {
    set_error("socket", "addressinuse");
}

/^Initial derivatives calculation [[:digit:]]+ does not converge; aborting\.\.\.$/ ||
/^Max iterations number [[:digit:]]+ has been reached during Step=[[:digit:]]+, Time=.+; TimeStep=.+ cannot be reduced further; aborting\.\.\.$/ ||
/^Initial assembly iterations reached maximum number [[:digit:]]; aborting\.\.\.$/ ||
/^Max iterations number [[:digit:]]+ has been reached during first step, Time=.+; TimeStep=.+ cannot be reduced further; aborting\.\.\.$/ {
    set_error("solver", "convergence");
}

/^MBDyn was interrupted$/ {
    set_error("interrupted", "");
}

/^An error occurred during the execution of MBDyn; aborting\.\.\. $/ {
    if (match(global::error_found, /\<unknown\>/)) {
        set_error("generic", "");
    }
}

ENDFILE {
    if (length(suppressed_errors) && global::error_found ~ suppressed_errors) {
        ignore = "suppressed";
    } else {
        ignore = "failed";
        status = 1;

        if (print_filename) {
            printf("%s:%s:%s:%s\n", FILENAME, global::error_found, global::message, ignore);
        } else {
            printf("%s:%s:%s\n", global::error_found, global::message, ignore);
        }
    }
}

END {
    exit status;
}
