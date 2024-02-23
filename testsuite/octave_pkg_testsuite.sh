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

## Octave package testsuite based on Octave's __run_test_suite__ function
## At the moment "mboct-mbdyn-pkg" is the only package which will perform unit tests on MBDyn

set -o pipefail ## Needed for commands like "octave --eval ${cmd} |& tee logfile"

program_name="$0"
OCT_PKG_LIST="${OCT_PKG_LIST:-mboct-mbdyn-pkg:no:master:yes:unlimited}"
OCT_PKG_TEST_DIR="${OCT_PKG_TEST_DIR:-octave-pkg-testsuite}"
OCTAVE_EXEC="${OCTAVE_EXEC:-octave}"
OCT_PKG_TESTS_VERBOSE="${OCT_PKG_TESTS_VERBOSE:-no}"
OCT_PKG_PRINT_RES="${OCT_PKG_PRINT_RES:-no}"
OCT_PKG_TEST_MODE="${OCT_PKG_TEST_MODE:-pkg}"
OCT_PKG_INSTALL_PREFIX="${OCT_PKG_INSTALL_PREFIX:-}"
OCTAVE_CMD_ARGS="-qfH"
OCT_PKG_FUNCTION_FILTER='/.+\.(tst|m)\>/'
## Do not print any output from Octave which does not pass through this filter, even if "--verbose yes" is used!
## This is strictly required because the amount of output is limited to 4194304 bytes by GitLab
OCT_GREP_FILTER_EXPR='^command: "mbdyn|^!!!!! test failed$|/^PASSES\>/|[[:alnum:]]+/[[:alnum:]]+/[[:alnum:]]+\.m\>|\<PASS\>|\<FAIL\>|\<pass\>|\<fail\>|^Summary|^Integrated test scripts|\.m files have no tests\.$'

## Disable multithreaded BLAS by default
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

## Will be used only by mboct-mbdyn-pkg
MBD_NUM_TASKS=${MBD_NUM_TASKS:-$(( $(lscpu | awk '/^Socket\(s\)/{ print $2 }') * $(lscpu | awk '/^Core\(s\) per socket/{ print $4 }') ))}
MBD_NUM_THREADS=${MBD_NUM_THREADS:-1}

echo $program_name

if test "$(basename ${program_name})" = "${program_name}" && ! test -z "$(which ${program_name})"; then
    ## Path of script was inside the environment variable PATH
    program_name="$(which ${program_name})"
    echo "${program_name} found on PATH"
fi

program_name=`realpath ${program_name}`

if ! test -z "${program_name}"; then
    export AWKPATH="$(dirname "${program_name}"):${AWKPATH}" ## Needed for parse_test_suite_status.awk
fi

while ! test -z "$1"; do
    case "$1" in
        --octave-pkg-list)
            OCT_PKG_LIST="$2"
            shift
            ;;
        --octave-pkg-test-dir)
            OCT_PKG_TEST_DIR="$2"
            shift
            ;;
        --octave-exec)
            OCTAVE_EXEC="$2"
            shift
            ;;
        --octave-pkg-test-mode)
            OCT_PKG_TEST_MODE="$2"
            shift
            ;;
        --octave-pkg-prefix)
            OCT_PKG_INSTALL_PREFIX="$2"
            shift
            ;;
        --octave-function-filter)
            if ! test -z "$2"; then
                OCT_PKG_FUNCTION_FILTER="$2"
            fi
            shift
            ;;
        --tasks)
            MBD_NUM_TASKS="$2"
            shift
            ;;
        --threads)
            MBD_NUM_THREADS="$2"
            shift
            ;;
        --verbose)
            OCT_PKG_TESTS_VERBOSE="$2"
            shift
            ;;
        --print-resources)
            OCT_PKG_PRINT_RES="$2"
            shift
            ;;
        --help|-h)
            printf "%s\n --octave-pkg-list <list-of-packages-and-flags>\n --octave-pkg-test-dir <output-dir>\n --print-resources {all|time|disc|profile}\n --octave-pkg-prefix <pkg-install-dir>\n --octave-exec <octave-executable>\n --octave-pkg-test-mode {pkg|single}\n --help\n" "${program_name}"
            exit 1;
            ;;
        *)
            printf "%s: invalid argument \"%s\"\n" "${program_name}" "$1"
            exit 1
            ;;
    esac
    shift
done

if ! test -d "${OCT_PKG_TEST_DIR}"; then
    if ! mkdir -p "${OCT_PKG_TEST_DIR}"; then
        echo "Failed to create directory \"${OCT_PKG_TEST_DIR}\"!"
        exit 1
    fi
fi

if ! test -d "${OCT_PKG_TEST_DIR}"; then
    echo "Directory \"${OCT_PKG_TEST_DIR}\" does not exist!"
    exit 1
fi

if ! test -z "${OCT_PKG_INSTALL_PREFIX}"; then
    if ! test -f "${OCT_PKG_INSTALL_PREFIX}/octave_packages"; then
        echo "Directory \"${OCT_PKG_INSTALL_PREFIX}\" is not valid!"
        exit 1
    fi
fi

case "${OCT_PKG_TEST_MODE}" in
    pkg|single)
        echo "Test mode: ${OCT_PKG_TEST_MODE}"
        ;;
    *)
        echo "Invalid value for --octave-pkg-test-mode \"${OCT_PKG_TEST_MODE}\""
        exit 1
esac

OCT_PKG_TEST_DIR="$(realpath ${OCT_PKG_TEST_DIR})"

test_status="passed"
failed_packages=""
octave_pkg_testsuite_pid=$$

case "${OCT_PKG_PRINT_RES}" in
    all|*profile*)
        oct_pkg_profile_on_cmd="profile('on');"
        oct_pkg_profile_off_fmt="profile('off');prof=profile('info');save('-binary','%s','prof');"
        oct_pkg_profile_post_fmt="load('%s');profshow(prof);"
        ;;
    *)
        oct_pkg_profile_on_cmd=""
        oct_pkg_profile_off_fmt=""
        oct_pkg_profile_post_fmt=""
        ;;
esac

function octave_pkg_testsuite_run()
{
    octave_code_cmd=""
    octave_status_file=""
    octave_pkg_name=""
    octave_pkg_task_id=""

    while ! test -z "$1"; do
        case "$1" in
            --exec)
                octave_code_cmd="$2"
                shift
                ;;
            --pkg)
                octave_pkg_name="$2"
                shift
                ;;
            --status)
                octave_status_file="$2"
                shift
                ;;
            --task-id)
                octave_pkg_task_id="$2"
                shift
                ;;
        esac
        shift
    done

    if test -z "${octave_code_cmd}"; then
        echo "Invalid argument --exec"
        return 1
    fi

    if test -z "${octave_status_file}"; then
        echo "Invalid argument --status"
        return 1
    fi

    if test -z "${octave_pkg_name}"; then
        echo "Invalid argument --pkg"
        return 1
    fi

    if test -z "${octave_pkg_task_id}"; then
        echo "Invalid argument --task-id"
        return 1
    fi

    rm -f "${octave_status_file}"

    ## Octave allows us to set TMPDIR in order to store all the temporary files in a single folder.
    ## This will make it easier to delete those files.
    export TMPDIR="${OCT_PKG_TEST_DIR}/${octave_pkg_name}/${octave_pkg_task_id}"

    if ! test -d "${TMPDIR}"; then
        if ! mkdir "${TMPDIR}"; then
            echo "Failed to create directory ${TMPDIR}"
            return 1
        fi
    fi

    octave_pkg_timing_file="${TMPDIR}/fntests.tm"

    case "${OCT_PKG_PRINT_RES}" in
        all|*time*)
            octave_pkg_timing_cmd="/usr/bin/time --verbose --output ${octave_pkg_timing_file}"
            ;;
        *)
            octave_pkg_timing_cmd=""
            ;;
    esac

    OCTAVE_CMD=$(printf '%s%s %s %s --eval %s' "${TIMEOUT_CMD}" "${octave_pkg_timing_cmd}" "${OCTAVE_EXEC}" "${OCTAVE_CMD_ARGS}" "${octave_code_cmd}")

    case "${OCT_PKG_PRINT_RES}" in
        all|*disk*)
            echo "Memory usage before test:"
            vmstat -S M
            echo "Temporary files before test:"
            ls -lhF "${TMPDIR}"
            echo "Disk usage before test:"
            df -h
            ;;
    esac

    curr_test_status="failed"

    pkg_test_output_file="${TMPDIR}/fntests.out"
    pkg_test_log_file="${TMPDIR}/fntests.log" ## created by __run_test_suite__

    ## Make sure that we do not read any old stuff ...
    rm -f "${pkg_test_output_file}"
    rm -f "${pkg_test_log_file}"

    if ! test -z "${octave_pkg_timing_file}"; then
        rm -f "${octave_pkg_timing_file}"
    fi

    echo "${OCTAVE_CMD}"

    curr_dir="`pwd`"

    ## Octave's __run_test_suite__ function will create the file "fntests.log" inside the current directory.
    if ! cd "${TMPDIR}"; then
        echo "Invalid directory ${TMPDIR}"
        return 1
    fi

    case "${OCT_PKG_TESTS_VERBOSE}" in
        yes)
            echo "Verbose output enabled:"
            ## If grep returns a nonzero status, this is not considered as an error!
            ${OCTAVE_CMD} 2>&1 | tee "${pkg_test_output_file}" | (grep -i -E "${OCT_GREP_FILTER_EXPR}" || true)
            rc=$?
            ;;
        *)
            echo "Verbose output disabled:"
            ${OCTAVE_CMD} >& "${pkg_test_output_file}"
            rc=$?
            ;;
    esac

    if ! cd "${curr_dir}"; then
        return 1
    fi

    case ${rc} in
        0)
            echo "${OCTAVE_CMD} completed with status 0"
            case ${OCT_PKG_TEST_MODE} in
                pkg|single)
                    if test -f "${pkg_test_output_file}"; then
                        pkg_test_log_parse="${pkg_test_output_file}"
                    else
                        pkg_test_log_parse="${pkg_test_log_file}"
                    fi

                    if test -f "${pkg_test_log_parse}"; then
                        if awk -f parse_test_suite_status.awk "${pkg_test_log_parse}"; then
                            curr_test_status="passed"
                        else
                            cat "${pkg_test_log_parse}"
                        fi
                    else
                        echo "File ${pkg_test_log_parse} not found"
                        curr_test_status="unexpected"
                    fi
                    ;;
            esac
            ;;
        124)
            echo "${OCTAVE_CMD} failed with timeout"
            ;;
    esac

    if test -f "${octave_pkg_timing_file}"; then
        echo "Resources used by ${OCTAVE_CMD}"
        cat "${octave_pkg_timing_file}"
        rm -f "${octave_pkg_timing_file}"
    fi

    case "${OCT_PKG_PRINT_RES}" in
        all|*disk*)
            echo "Memory usage after test:"
            vmstat -S M
            echo "Temporary files after test:"
            ls -lhF ${OCT_PKG_TEST_DIR}
            echo "Disk usage after test:"
            df -h
            ;;
    esac

    case "${curr_test_status}" in
        passed)
            printf "octave testsuite for package \"%s\" passed\n" "${octave_pkg_name}"
            ;;
        *)
            printf "octave testsuite for package \"%s\" failed\n" "${octave_pkg_name}"
            if test -f "${pkg_test_log_file}"; then
                cat "${pkg_test_log_file}";
            else
                echo "${pkg_test_log_file} not found";
            fi
            ;;
    esac

    printf "%d:%s:%s:%s\n" "${rc}" "${curr_test_status}" "${octave_pkg_name}" "${octave_code_cmd}" > "${octave_status_file}"

    echo "Remove all temporary files after the test:"
    ## "oct-" is the default prefix of Octave's tempname() function
    export -n TMPDIR
    find "${TMPDIR}" '(' -type f -and '(' -name 'oct-*' -or -name 'fntests.*' ')' ')' -delete
}

for pkgname_and_flags in ${OCT_PKG_LIST}; do
    pkgname=$(echo ${pkgname_and_flags} | awk -F ":" "{print \$1}")
    pkg_test_flag=$(echo ${pkgname_and_flags} | awk -F ":" "{print \$4}")
    pkg_test_timeout=$(echo ${pkgname_and_flags} | awk -F ":" "{print \$5}")

    if test -z "${pkg_test_timeout}"; then
        pkg_test_timeout="unlimited"
    fi

    case "${pkg_test_flag}" in
        yes)
            printf "Test package \"%s\"\n" "${pkgname}"
            ;;
        *)
            printf "Tests of package \"%s\" skipped\n" "${pkgname}"
            continue
            ;;
    esac

    echo "Temporary directory: ${OCT_PKG_TEST_DIR}/${pkgname}"

    if ! test -d "${OCT_PKG_TEST_DIR}/${pkgname}"; then
        if ! mkdir "${OCT_PKG_TEST_DIR}/${pkgname}"; then
            echo "Failed to create directory ${OCT_PKG_TEST_DIR}/${pkgname}"
            exit 1
        fi
    fi

    if test -z "${OCT_PKG_INSTALL_PREFIX}"; then
        oct_pkg_prefix_cmd=""
    else
        oct_pkg_prefix_cmd=$(printf "pkg('local_list','%s');" "${OCT_PKG_INSTALL_PREFIX}/octave_packages")
    fi

    oct_pkg_sigterm_dumps_core="sigterm_dumps_octave_core(false);"
    oct_pkg_load_cmd=$(printf "pkg('load','%s');" "${pkgname}")
    oct_pkg_list_cmd=$(printf "p=pkg('list','%s');" "${pkgname}")
    oct_pkg_run_test_suite_cmd="__run_test_suite__({p{1}.dir},{p{1}.dir});"

    case "${OCT_PKG_TEST_MODE}" in
        pkg)
            oct_pkg_profile_data="${OCT_PKG_TEST_DIR}/oct_pkg_profile_data_${octave_pkg_testsuite_pid}_${pkgname}.mat"
            oct_pkg_profile_off_cmd=$(printf "${oct_pkg_profile_off_fmt}" "${oct_pkg_profile_data}")
            OCTAVE_CODE="${oct_pkg_sigterm_dumps_core}${oct_pkg_prefix_cmd}${oct_pkg_load_cmd}${oct_pkg_list_cmd}${oct_pkg_profile_on_cmd}${oct_pkg_run_test_suite_cmd}${oct_pkg_profile_off_cmd}"
            ;;
        single)
            OCTAVE_CMD_FUNCTIONS=$(printf "p=pkg('list','-verbose','%s');dir(fullfile(p{1}.dir,'*.m'));dir(fullfile(p{1}.dir,'*.tst'));" "${pkgname}")
            OCTAVE_PKG_FUNCTIONS=`${OCTAVE_EXEC} ${OCTAVE_CMD_ARGS} --eval "${oct_pkg_prefix_cmd}${OCTAVE_CMD_FUNCTIONS}" | awk "${OCT_PKG_FUNCTION_FILTER}"`
            rc=$?
            if test ${rc} != 0; then
                echo "${OCTAVE_CMD_FUNCTIONS} failed with status ${rc}"
                exit 1
            fi
            OCTAVE_CODE=""
            ((oct_pkg_func_index=0))
            for pkg_function_name in ${OCTAVE_PKG_FUNCTIONS}; do
                ((++oct_pkg_func_index))
                oct_pkg_test_function_cmd=$(printf "test('%s');" "${pkg_function_name}")
                oct_pkg_profile_data=`printf '%s/oct_pkg_profile_data_%d_%s_%03d.mat' "${OCT_PKG_TEST_DIR}" ${octave_pkg_testsuite_pid} "${pkgname}" $((oct_pkg_func_index))`
                oct_pkg_profile_off_cmd=$(printf "${oct_pkg_profile_off_fmt}" "${oct_pkg_profile_data}")
                OCTAVE_CODE="${OCTAVE_CODE} ${oct_pkg_sigterm_dumps_core}${oct_pkg_prefix_cmd}${oct_pkg_load_cmd}${oct_pkg_profile_on_cmd}${oct_pkg_test_function_cmd}${oct_pkg_profile_off_cmd}"
            done
            ;;
        *)
            exit 1
            ;;
    esac

    case "${pkg_test_timeout}" in
        unlimited)
            TIMEOUT_CMD=""
            ;;
        *)
            ## Need to call Octave's "sigterm_dumps_octave_core" in order to avoid any core dumps!
            TIMEOUT_CMD="timeout --signal=SIGTERM ${pkg_test_timeout} "
            ;;
    esac

    octave_status_file_format=`printf '%s/octave_pkg_testsuite_run_%08X_%%s.status' "${OCT_PKG_TEST_DIR}/${pkgname}" "${octave_pkg_testsuite_pid}"`

    oct_pkg_num_tests=`echo ${OCTAVE_CODE} | wc -w`

    if test ${MBD_NUM_TASKS} -gt 1 && test ${oct_pkg_num_tests} -gt 1; then
        echo "Parallel execution using ${MBD_NUM_TASKS} tasks:"
        export MBD_NUM_THREADS
        export TIMEOUT_CMD
        export OCTAVE_EXEC
        export OCTAVE_CMD_ARGS
        export OCT_PKG_PRINT_RES
        export OCT_PKG_TEST_DIR
        export OCT_GREP_FILTER_EXPR
        export OCT_PKG_TEST_MODE
        export OCT_PKG_TESTS_VERBOSE
        export -f octave_pkg_testsuite_run
        octave_status_file=`printf "${octave_status_file_format}" '{#}'`
        octave_parallel_args="-j${MBD_NUM_TASKS} -n1 octave_pkg_testsuite_run --status ${octave_status_file} --exec {} --pkg ${pkgname} --task-id {#}"
        printf '%s\n' ${OCTAVE_CODE} | parallel ${octave_parallel_args}
    else
        echo "Serial execution:"
        ((idx_test=0))
        for octave_code_cmd in ${OCTAVE_CODE}; do
            ((++idx_test))
            status_file=`printf ${octave_status_file_format} $((idx_test))`
            octave_pkg_testsuite_run --exec "${octave_code_cmd}" --pkg "${pkgname}" --status "${status_file}" --task-id $((idx_test))
        done
    fi

    ((idx_test=0))
    ((tests_passed=0))
    ((tests_failed=0))

    for octave_code_cmd in ${OCTAVE_CODE}; do
        ((++idx_test))
        status_file=`printf ${octave_status_file_format} $((idx_test))`

        if ! test -f "${status_file}"; then
            echo "Status file ${status_file} not found"
            ((++tests_failed))
            continue
        fi

        curr_test_status=`awk -F ':' '{print $2}' ${status_file}`

        if ! test "${curr_test_status}" = "passed"; then
            ((++tests_failed))
        fi

        rm -rf "${status_file}"

        ((++tests_passed))
    done

    printf "%s:%d/%d tests passed\n" "${pkgname}" $((tests_passed)) $((tests_passed+tests_failed))
    printf "%s:%d/%d tests failed\n" "${pkgname}" $((tests_failed)) $((tests_passed+tests_failed))

    if test $((tests_failed)) -gt 0 || test $((tests_passed)) -le 0; then
        test_status="failed"
        failed_packages="${failed_packages} ${pkgname}"
    fi
done

if ! test -z "${oct_pkg_profile_post_fmt}"; then
    oct_pkg_profile_files=`find "${OCT_PKG_TEST_DIR}" '(' -type f -and -name "oct_pkg_profile_data_${octave_pkg_testsuite_pid}_*.mat" ')'`

    for oct_pkg_profile_file in ${oct_pkg_profile_files}; do
        oct_pkg_profile_post_cmd=$(printf "${oct_pkg_profile_post_fmt}" "${oct_pkg_profile_file}")
        echo "Profile information for \"${oct_pkg_profile_file}\":"
        ${OCTAVE_EXEC} ${OCTAVE_CMD_ARGS} --eval "${oct_pkg_profile_post_cmd}"
        rm -f "${oct_pkg_profile_file}"
    done
fi

case "${test_status}" in
    passed)
        echo "All tests passed"
        exit 0
        ;;
    *)
        echo "The following packages did not pass:"
        for pkgname in ${failed_packages}; do
            printf "%s\n" "${pkgname}"
        done
        exit 1
        ;;
esac
