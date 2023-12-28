#!/bin/bash -f

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

## Simple testsuite just to see if a model can run

program_name="$0"
mbdyn_module_test_timeout="unlimited"
mbdyn_testsuite_prefix_output=""
mbdyn_testsuite_prefix_input=""

while ! test -z "$1"; do
    case "$1" in
        --prefix-output|-o)
            mbdyn_testsuite_prefix_output="$2"
            shift
            ;;
        --prefix-input|-i)
            mbdyn_testsuite_prefix_input="$2"
            shift
            ;;
        --timeout|-t)
            mbdyn_module_test_timeout="$2"
            shift
            ;;
        --help|-h)
            printf "%s\n  --prefix-output <output-dir>\n  --prefix-input <input-dir>\n  --timeout <timeout-seconds>\n  --help\n" "${program_name}" > /dev/stderr
            exit 1;
            ;;
        *)
            printf "%s: invalid argument \"%s\"\n" "${program_name}" "$1" > /dev/stderr
            exit 1
            ;;
    esac
    shift
done

if test -z "${mbdyn_module_test_timeout}"; then
    mbdyn_module_test_timeout="unlimited"
fi

if test -z "${mbdyn_testsuite_prefix_output}"; then
    printf "%s: missing argument --prefix-output <output-dir>\n" "${program_name}" > /dev/stderr
    exit 1
fi

if test -z "${mbdyn_testsuite_prefix_input}"; then
    printf "%s: missing argument --prefix-input <input-dir>\n" "${program_name}" > /dev/stderr
    exit 1
fi

if ! test -d "${mbdyn_testsuite_prefix_input}"; then
    printf "%s: argument --prefix-input \"%s\" is not valid\n" "${program_name}" "${mbdyn_testsuite_prefix_input}" > /dev/stderr
    exit 1
fi

passed_tests=""
failed_tests=""
killed_tests=""

for mbd_filename in `find "${mbdyn_testsuite_prefix_input}" -type f -print0 | xargs -0 awk "/begin: initial value;/{print FILENAME}"`; do
    mbd_basename=`basename -s ".mbdyn" "${mbd_filename}"`
    mbd_basename=`basename -s ".mbd" "${mbd_basename}"`
    if test -f "${mbd_filename}"; then
        mbd_time_file="${mbdyn_testsuite_prefix_output}${mbd_basename}_time.log"
        mbd_output_file="${mbdyn_testsuite_prefix_output}${mbd_basename}_mbdyn_output"
        mbd_log_file="${mbd_output_file}.stdout"
        echo ${mbd_log_file}
        mbd_command="ulimit -t ${mbdyn_module_test_timeout}; /usr/bin/time --verbose --output \"${mbd_time_file}\" mbdyn -C -f \"${mbd_filename}\" -o \"${mbd_output_file}\""
        printf "%s:%s\n" "${mbd_basename}" "${mbd_command}"
        status="failed"
        sh -c "${mbd_command}" >& "${mbd_log_file}"

        rc=$?

        case ${rc} in
            0)
                status="passed"
                ;;
            137)
                status="killed"
                ;;
            *)
                status="failed"
                ;;
        esac
        printf "test %s:%d\n" "${status}" "${rc}"
        rm -f "${mbdyn_testsuite_prefix_output}${mbd_basename}*"
        case "${status}" in
            passed)
                passed_tests="${passed_tests} ${mbd_filename}"
                ;;
            killed)
                killed_tests="${killed_tests} ${mbd_filename}"
                ;;
            *)
                failed_tests="${failed_tests} ${mbd_filename}"
                ;;
        esac
    fi
done

if test -z "${passed_tests}"; then
    echo "No tests passed"
else
    echo "The following tests passed:"
    for mbd_filename in ${passed_tests}; do
        echo "  ${mbd_filename}"
    done
fi

if test -z "${killed_tests}"; then
    echo "No tests were killed because of timeout"
else
    echo "The following tests were killed because of timeout:"
    for mbd_filename in ${killed_tests}; do
        echo " ${mbd_filename}"
    done
fi

if test -z "${failed_tests}"; then
    echo "No tests failed"
else
    echo "The following tests failed:"
    for mbd_filename in ${failed_tests}; do
        echo " ${mbd_filename}"
    done
    exit 1
fi
