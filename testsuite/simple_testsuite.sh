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

set -o pipefail ## Needed for commands like "mbdyn -f input_file |& tee logfile"

program_name="$0"
mbdyn_testsuite_timeout="unlimited"
mbdyn_testsuite_prefix_output=""
mbdyn_testsuite_prefix_input=""
mbdyn_input_filter=""
mbdyn_verbose_output="no"
mbdyn_keep_output="unexpected"
declare -i mbd_exit_status_mask=0 ## Define the errors codes which should not cause the pipeline to fail
declare -i mbd_test_idx_start=1
declare -i mbd_test_idx_offset=1
OCTAVE_EXEC="${OCTAVE_EXEC:-octave}"

## Disable multithreaded BLAS by default
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

## Might be used for Octave scripts (e.g. via mboct-mbdyn-pkg)
export MBD_NUM_THREADS=${MBD_NUM_THREADS:-`awk -F ':' 'BEGIN{cores=-1;}/^cpu cores\>/{cores=strtonum($2);}END {print cores;}' /proc/cpuinfo`}

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
        --test-index-start)
            mbd_test_idx_start="$2"
            shift
            ;;
        --test-index-offset)
            mbd_test_idx_offset="$2"
            shift
            ;;
        --verbose)
            mbdyn_verbose_output="$2"
            shift
            ;;
        --keep-output)
            mbdyn_keep_output="$2"
            shift
            ;;
        --help)
            printf "%s\n  --prefix-output <output-dir>\n  --prefix-input <input-dir>\n  --timeout <timeout-seconds>\n  --help\n" "${program_name}"
            exit 1;
            ;;
        --exit-status-mask)
            ((mbd_exit_status_mask=$2))
            shift
            ;;
        *)
            printf "%s: invalid argument \"%s\"\n" "${program_name}" "$1"
            exit 1
            ;;
    esac
    shift
done

if test ${mbd_test_idx_start} -lt 1; then
    printf "%s: invalid argument --test-index-start %d\n" "${program_name}" ${mbd_test_idx_start}
    exit 1
fi

if test ${mbd_test_idx_offset} -lt 1; then
    printf "%s: invalid argument --test-index-offset %d\n" "${program_name}" ${mbd_test_idx_offset}
    exit 1
fi

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
modules_not_found=""
loadables_not_found=""
unexpected_faults=""

search_expression="-type f"

if ! test -z "${mbdyn_input_filter}"; then
    search_expression=`printf -- "-name %s -and %s" "${mbdyn_input_filter}" "${search_expression}"`
fi

## Octave allows us to set TMPDIR in order to store all the temporary files in a single folder.
## This will make it easier to delete those files.
export TMPDIR="${mbdyn_testsuite_prefix_output}"
declare -i idx_test=0

for mbd_filename in `find ${mbdyn_testsuite_prefix_input} '(' ${search_expression} ')' -print0 | xargs -0 awk "/begin: initial value;/{print FILENAME}"`; do
    ((++idx_test))

    printf "%4d: \"%s\"\n" $((idx_test)) "${mbd_filename}"

    if test $((idx_test)) -lt $((mbd_test_idx_start)); then
        echo "  skipped ..."
        continue
    fi

    ((mbd_test_idx_start+=mbd_test_idx_offset))

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
        mbd_script_name_m1="${mbd_dir_name}/${mbd_script_name}_run.m"
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
                    *_run.m)
                        mbd_command="${OCTAVE_EXEC} -q -f ${mbd_script_name} -f ${mbd_filename} -o ${mbd_output_file}"
                        ;;
                    *_run.sh)
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
                mbd_command="timeout --signal=SIGTERM ${mbdyn_testsuite_timeout} ${mbd_command}"
            ;;
        esac

        printf "%s:%s\n" "${mbd_basename}" "${mbd_command}"

        curr_dir="$(pwd)"

        if ! cd "${mbd_dir_name}"; then
            exit 1
        fi

        rm -f "${mbd_time_file}"
        rm -f "${mbd_log_file}"

        if test "${mbdyn_verbose_output}" = "yes"; then
            ${mbd_command} |& tee "${mbd_log_file}"
        else
            ${mbd_command} >& "${mbd_log_file}"
        fi

        rc=$?

        if ! cd "${curr_dir}"; then
            exit 1
        fi

        mbd_module_not_found=""
        mbd_loadable_not_found=""

        case ${rc} in
            0)
                num_steps=`awk 'BEGIN{num_steps=0}/^End of simulation at time [0-9.-]+ after [0-9]+ steps;$/{num_steps=$8} END{print num_steps}' "${mbd_log_file}"`
                status=$(printf 'passed{Steps=%d}' "${num_steps}")
                ;;
            1)
                mbd_module_not_found=`awk -v FPAT='[^:<>]+' '/^ModuleLoad_int: unable to open module\>/{print $3;}' "${mbd_log_file}"`
                mbd_loadable_not_found=`awk -v FPAT='[^<>]+' '/^ParseUserDefinedElem\([0-9]+\): unknown user-defined element type at line [0-9]+, file <.+>$/{ print $2 }' "${mbd_log_file}"`

                if ! test -z "${mbd_module_not_found}"; then
                    status="module"
                else
                    if ! test -z "${mbd_loadable_not_found}"; then
                        status="loadable"
                    else
                        status="failed"
                    fi
                fi
                ;;
            124)
                status="timeout"
                echo "Timeout"
                ;;
            130)
                echo "Interrupted"
                exit ${rc}
                ;;
            143)
                echo "Terminated"
                exit ${rc}
                ;;
            137)
                echo "Killed"
                exit ${rc}
                ;;
            *)
                status="unexpected"
                ;;
        esac

        printf "Test \"%s\" %s with status %d\n" "${mbd_basename}" "${status}" "${rc}"

        if test -f "${mbd_time_file}"; then
            cat "${mbd_time_file}"
        fi

        keep_output_flag="no"

        case "${mbdyn_keep_output}" in
            all)
                keep_output_flag="yes"
                ;;
            failed)
                case "${status}" in
                    failed|unexpected)
                        keep_output_flag="yes"
                        ;;
                esac
                ;;
            unexpected)
                case "${status}" in
                    unexpected)
                        keep_output_flag="yes"
                        ;;
                esac
                ;;
        esac

        mbd_output_file=$(awk -F '"' '/^output in file\>/{print $2}' "${mbd_log_file}")

        if test "${keep_output_flag}" = "no"; then
            if test -f "${mbd_output_file}.log"; then
                find "${mbdyn_testsuite_prefix_output}" '(' -type f -and -wholename $(printf '%s.*' "${mbd_output_file}") ')' -delete
            fi

            rm -f "${mbd_log_file}"
            rm -f "${mbd_time_file}"
        fi

        case "${status}" in
            passed*)
                passed_tests="${passed_tests} ${mbd_filename}:${status}(${rc})"
                ;;
            timeout)
                timeout_tests="${timeout_tests} ${mbd_filename}:${status}(${rc})"
                ;;
            module)
                modules_not_found="${modules_not_found} ${mbd_filename}[${mbd_module_not_found}]:${status}(${rc})"
                ;;
            loadable)
                loadables_not_found="${loadables_not_found} ${mbd_filename}:${status}(${rc})"
                ;;
            failed)
                failed_tests="${failed_tests} ${mbd_filename}:${status}(${rc})"
                ;;
            *)
                unexpected_faults="${unexpected_faults} ${mbd_filename}:${status}(${rc})"
                ;;
        esac
    fi
done

declare -i exit_status=0x0

function print_files()
{
    declare -i idx=0
    format="$1"
    shift

    printf "${format}" "$(echo ${*} | wc -w)"

    for mbd_filename in $*; do
        printf "%4d:%s\n" "$((++idx))" "${mbd_filename}"
    done
}

if test -z "${passed_tests}"; then
    echo "No tests passed"
    ((exit_status|=0x1))
else
    print_files "PASSED:The following %d tests passed with zero exit status:\n" ${passed_tests}
fi

if test -z "${timeout_tests}"; then
    echo "No tests were killed because of timeout"
else
    print_files "TIMEOUT:The following %d tests were killed because of timeout:\n" ${timeout_tests}
    ((exit_status|=0x2))
fi

if test -z "${modules_not_found}"; then
    echo "All modules were found"
else
    print_files "FAILED-MODULE:The following %d tests failed because a loadable module was not found:\n" ${modules_not_found}
    ((exit_status|=0x4))
fi

if test -z "${loadables_not_found}"; then
    echo "All loadables were found"
else
    print_files "FAILED-LOADABLE:The following %d tests failed because a loadable element was not found:\n" ${loadables_not_found}
    ((exit_status|=0x8))
fi

if test -z "${failed_tests}"; then
    echo "No tests failed with status 1"
else
    print_files "FAILED:The following %d tests failed with status 1:\n" ${failed_tests}
    ((exit_status|=0x10))
fi

if test -z "${unexpected_faults}"; then
    echo "No tests returned with unexpected exit status"
else
    print_files "FAILED-UNEXPECTED:The following %d tests failed with unexpected exit status:\n" ${unexpected_faults}
    ((exit_status|=0x40))
fi

if test $((exit_status&~mbd_exit_status_mask)) == 0 && test $((exit_status)) != 0; then
    printf "Minor failures (0x%X) were detected during this test, but they will be ignored!\n" $((exit_status&mbd_exit_status_mask))
fi

## Mask all the error codes which should be ignored for the current test
((exit_status&=~mbd_exit_status_mask))

printf "${program_name} exit status 0x%X\n" $((exit_status))

exit $((exit_status))
