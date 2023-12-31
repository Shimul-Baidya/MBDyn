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

## Use the simple testsuite in order to check if a model can run or not

program_name="$0"
mbdyn_testsuite_timeout="unlimited"
mbdyn_testsuite_prefix_output=""
mbdyn_testsuite_prefix_input=""
mbdyn_input_filter=""
mbdyn_verbose_output="no"
OCTAVE_EXEC="${OCTAVE_EXEC:-octave}"

## Do not use multithreaded BLAS by default, because this could cause performance issues!
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

while ! test -z "$1"; do
    case "$1" in
        --prefix-output)
            mbdyn_testsuite_prefix_output="$2"
            shift
            ;;
        --prefix-input)
            mbdyn_testsuite_prefix_input="$2"
            shift
            ;;
        --timeout)
            mbdyn_testsuite_timeout="$2"
            shift
            ;;
        --filter)
            mbdyn_input_filter="$2"
            shift
            ;;
        --verbose)
            mbdyn_verbose_output="$2"
            shift
            ;;
        --help)
            printf "%s\n  --prefix-output <output-dir>\n  --prefix-input <input-dir>\n  --timeout <timeout-seconds>\n  --help\n" "${program_name}"
            exit 1;
            ;;
        *)
            printf "%s: invalid argument \"%s\"\n" "${program_name}" "$1"
            exit 1
            ;;
    esac
    shift
done

if test -z "${mbdyn_testsuite_timeout}"; then
    mbdyn_testsuite_timeout="unlimited"
fi

if test -z "${mbdyn_testsuite_prefix_output}"; then
    printf "%s: missing argument --prefix-output <output-dir>\n" "${program_name}"
    exit 1
fi

mbdyn_testsuite_prefix_output="$(realpath ${mbdyn_testsuite_prefix_output})"

if ! test -d "${mbdyn_testsuite_prefix_output}"; then
    if ! mkdir -p "${mbdyn_testsuite_prefix_output}"; then
        echo "Failed to create directory \"${mbdyn_testsuite_prefix_output}\""
        exit 1
    fi
fi

if test -z "${mbdyn_testsuite_prefix_input}"; then
    printf "%s: missing argument --prefix-input <input-dir>\n" "${program_name}"
    exit 1
fi

if ! test -d "${mbdyn_testsuite_prefix_input}"; then
    printf "%s: argument --prefix-input \"%s\" is not valid\n" "${program_name}" "${mbdyn_testsuite_prefix_input}"
    exit 1
fi

mbdyn_testsuite_prefix_input=`realpath ${mbdyn_testsuite_prefix_input}`

passed_tests=""
failed_tests=""
timeout_tests=""

search_expression="-type f"

if ! test -z "${mbdyn_input_filter}"; then
    search_expression=`printf -- "-name %s -and %s" "${mbdyn_input_filter}" "${search_expression}"`
fi

## Octave allows us to set TMPDIR in order to store all the temporary files in a single folder.
## This will make it easier to delete those files.
export TMPDIR="${mbdyn_testsuite_prefix_output}"

for mbd_filename in `find ${mbdyn_testsuite_prefix_input} '(' ${search_expression} ')' -print0 | xargs -0 awk "/begin: initial value;/{print FILENAME}"`; do
    mbd_basename=`basename -s ".mbdyn" "${mbd_filename}"`
    mbd_basename=`basename -s ".mbd" "${mbd_basename}"`

    if test -f "${mbd_filename}"; then
        mbd_time_file="${mbdyn_testsuite_prefix_output}/${mbd_basename}_time.log"
        mbd_output_file="${mbdyn_testsuite_prefix_output}/${mbd_basename}_mbdyn_output"
        mbd_log_file="${mbd_output_file}.stdout"
        echo ${mbd_log_file}

        mbd_script_name=`basename ${mbd_filename}`
        mbd_script_name=`basename -s .mbd ${mbd_script_name}`
        mbd_script_name=`basename -s .mbdyn ${mbd_script_name}`
        mbd_dir_name=`dirname ${mbd_filename}`
        mbd_script_name_sh="${mbd_dir_name}/${mbd_script_name}_run.sh"
        mbd_script_name_m1="${mbd_dir_name}/${mbd_script_name}.m"
        mbd_script_name_m2="${mbd_dir_name}/${mbd_script_name}_gen.m"
        mbd_command=""

        for mbd_script_name in "${mbd_script_name_m2}" "${mbd_script_name_m1}" "${mbd_script_name_sh}"; do
            if ! test -z "${mbd_command}"; then
                break
            fi
            if test -f "${mbd_script_name}"; then
                echo "A custom test script ${mbd_script_name} was found for input file ${mbd_filename}; It will be used to run the model"
                case "${mbd_script_name}" in
                    *_gen.m)
                        mbd_command="${OCTAVE_EXEC} -q -f ${mbd_script_name} -f ${mbd_filename} -o ${mbd_output_file}; mbdyn -C -f ${mbd_filename} -o ${mbd_output_file}"
                        ;;
                    *.m)
                        mbd_command="${OCTAVE_EXEC} -q -f ${mbd_script_name} -f ${mbd_filename} -o ${mbd_output_file}"
                        ;;
                    *.sh)
                        chmod +x "${mbd_script_name}"
                        mbd_command="${mbd_script_name} -f ${mbd_filename} -o ${mbd_output_file}"
                        ;;
                esac
            fi
        done

        if test -z "${mbd_command}"; then
            echo "No custom test script was found for input file ${mbd_filename}; The default command will be used to run the model"
            mbd_command="mbdyn -C -f ${mbd_filename} -o ${mbd_output_file}"
        fi

        mbd_command="/usr/bin/time --verbose --output ${mbd_time_file} ${mbd_command}"
        status="failed"

        if test -z "${mbdyn_testsuite_timeout}"; then
            mbdyn_testsuite_timeout="unlimited"
        fi

        case "${mbdyn_testsuite_timeout}" in
            unlimited)
                echo "no timeout is applied"
                ;;
            *)
                echo "timeout after ${mbdyn_testsuite_timeout}"
                ## Octave will not dump a core file on SIGINT
                mbd_command="timeout --signal=SIGINT ${mbdyn_testsuite_timeout} ${mbd_command}"
            ;;
        esac

        printf "%s:%s\n" "${mbd_basename}" "${mbd_command}"

        curr_dir="$(pwd)"

        if ! cd "${mbd_dir_name}"; then
            exit 1
        fi

        if test "${mbdyn_verbose_output}" = "yes"; then
            ${mbd_command}
        else
            ${mbd_command} >& "${mbd_log_file}"
        fi

        rc=$?

        if ! cd "${curr_dir}"; then
            exit 1
        fi

        case ${rc} in
            0)
                status="passed"
                ;;
            124)
                status="timeout"
                ;;
            130)
                echo "interrupted"
                exit ${rc}
                ;;
            *)
                status="failed"
                ;;
        esac
        printf "Test \"%s\" %s with status %d\n" "${mbd_basename}" "${status}" "${rc}"
        echo rm -f "${mbdyn_testsuite_prefix_output}/${mbd_basename}*"
        case "${status}" in
            passed)
                passed_tests="${passed_tests} ${mbd_filename}"
                ;;
            timeout)
                timeout_tests="${timeout_tests} ${mbd_filename}"
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

if test -z "${timeout_tests}"; then
    echo "No tests were killed because of timeout"
else
    echo "The following tests were killed because of timeout:"
    for mbd_filename in ${timeout_tests}; do
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
