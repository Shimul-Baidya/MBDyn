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

program_name="$0"
program_dir=$(realpath $(dirname "${program_name}"))

if ! test -f "${program_dir}/simple_testsuite.sh"; then
    program_dir=$(realpath $(which "${program_name}"))
fi

if test -f "${program_dir}/mbdyn_input_file_format.awk"; then
    export PATH=${program_dir}:${PATH}
fi

mbdyn_testsuite_prefix_output=""
mbdyn_keep_output="unexpected"
declare -i mbd_exit_status_mask=0
other_arguments=""

while ! test -z "$1"; do
    case "$1" in
        --prefix-output)
            mbdyn_testsuite_prefix_output="$2"
            shift
            ;;
        --patch-input)
            shift
            ;;
        --keep-output)
            other_arguments="${other_arguments} $1 $2"
            mbdyn_keep_output="$2"
            shift
            ;;
        --exit-status-mask)
            ((mbd_exit_status_mask=$2))
            shift
            ;;
        *)
            other_arguments="${other_arguments} $1 $2"
            shift
            ;;
    esac
    shift
done

((mbd_exit_status_mask|=0x1))
failed_tests=""

for mbd_linear_solver in 'naive' 'umfpack' 'klu' 'pardiso' 'pardiso_64' 'y12' 'spqr' 'qr' 'lapack'; do
    for mbd_mh_type in 'map' 'cc' 'dir' 'grad'; do
        for mbd_mat_scale in 'rowmaxcolumnmax' 'iterative' 'lapack' 'rowmax' 'columnmax' 'rowsum' 'columnsum'; do
            for mbd_mat_scale_when in 'never' 'always' 'once'; do
                case "${mbd_mat_scale_when}" in
                    never)
                        case "${mbd_mat_scale}" in
                            rowmaxcolumnmax)
                                ;;
                            *)
                                continue
                                ;;
                        esac
                        ;;
                esac
                case "${mbd_linear_solver}" in
                    naive|lapack|qr)
                        case "${mbd_mh_type}" in
                            map)
                            ;;
                            *)
                                continue
                                ;;
                        esac
                        ;;
                    y12)
                        case "${mbd_mh_type}" in
                            map|cc|dir)
                            ;;
                            *)
                                continue
                                ;;
                        esac
                        ;;
                    pardiso|pardiso_64|spqr)
                        case "${mbd_mh_type}" in
                            map|grad)
                            ;;
                            *)
                                continue
                                ;;
                        esac
                        case "${mbd_mat_scale_when}" in
                            never)
                                ;;
                            *)
                                continue
                                ;;
                        esac
                esac
                for mbd_use_autodiff in 'autodiff' 'noautodiff'; do
                    case "${mbd_mh_type}" in
                        grad)
                            case "${mbd_use_autodiff}" in
                                noautodiff)
                                    continue
                                    ;;
                            esac
                            ;;
                    esac
                    for mbd_nonlin_solver in 'newtonraphson' 'linesearch' 'nox' 'nox-newton-krylov' 'nox-direct' 'nox-broyden' 'mcpnewtonminfb' 'mcpnewtonfb' 'bfgs'; do
                        case "${mbd_nonlin_solver}" in
                            mcpnewtonminfb)
                                case "${mbd_use_autodiff}" in
                                    noautodiff)
                                        continue
                                        ;;
                                esac
                                ;;
                            nox|nox-direct|nox-broyden)
                                case "${mbd_linear_solver}" in
                                    naive|qr|lapack)
                                        continue
                                        ;;
                                esac
                                ;;
                            bfgs)
                                case "${mbd_linear_solver}" in
                                    spqr|qr)
                                    ;;
                                    *)
                                        continue
                                        ;;
                                esac
                        esac

                        case "${mbd_linear_solver}" in
                            naive)
                                mbd_linear_solver_flags=",colamd"
                                ;;
                            *)
                                mbd_linear_solver_flags=""
                                ;;
                        esac

                        case "${mbd_nonlin_solver}" in
                            nox)
                                mbd_nonlin_solver_flags="nox, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            nox-newton-krylov)
                                mbd_nonlin_solver_flags="nox, modified, 10, jacobian operator, newton krylov, use preconditioner as solver, no, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            nox-direct)
                                mbd_nonlin_solver_flags="nox, use preconditioner as solver, yes, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            nox-broyden)
                                mbd_nonlin_solver_flags="nox, modified, 10, direction, broyden, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            *)
                                mbd_nonlin_solver_flags="${mbd_nonlin_solver}"
                                ;;
                        esac

                        mbd_output_dir="${mbdyn_testsuite_prefix_output}/${mbd_linear_solver}/${mbd_mh_type}/${mbd_mat_scale}/${mbd_mat_scale_when}/${mbd_use_autodiff}/${mbd_nonlin_solver}"
                        mkdir -p "${mbd_output_dir}"
                        export MBD_TESTSUITE_INITIAL_VALUE_BEGIN="${mbd_output_dir}/mbd_init_val_begin.set"
                        export MBD_TESTSUITE_INITIAL_VALUE_END="${mbd_output_dir}/mbd_init_val_end.set"
                        export MBD_TESTSUITE_CONTROL_DATA_BEGIN="${mbd_output_dir}/mbd_control_data_begin.set"
                        export MBD_TESTSUITE_CONTROL_DATA_END="${mbd_output_dir}/mbd_control_data_end.set"
                        printf '    # mbd_init_val_begin.set currently not used!\n' > "${MBD_TESTSUITE_INITIAL_VALUE_BEGIN}"
                        printf '    linear solver: %s,%s%s,scale,%s,%s;\n' "${mbd_linear_solver}" "${mbd_mh_type}" "${mbd_linear_solver_flags}" "${mbd_mat_scale}" "${mbd_mat_scale_when}" > "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                        printf '    abort after: derivatives;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                        printf '    threads: disable;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                        printf '    nonlinear solver: %s;\n' "${mbd_nonlin_solver_flags}" >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                        printf '    output: iterations;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                        printf '    derivatives max iterations: 10;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                        printf '    derivatives coefficient: auto;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                        case ${mbd_use_autodiff} in
                            autodiff)
                                mbd_use_autodiff_cmd="use automatic differentiation;"
                                ;;
                            *)
                                mbd_use_autodiff_cmd="# automatic differentiation disabled"
                                ;;
                        esac

                        printf '    # mbd_control_data_begin.set currently not used!\n' > "${MBD_TESTSUITE_CONTROL_DATA_BEGIN}"
                        printf '    %s\n' "${mbd_use_autodiff_cmd}" > "${MBD_TESTSUITE_CONTROL_DATA_END}"

                        simple_testsuite_log_file="${mbd_output_dir}/mbdyn-testsuite-patched.log"
                        cat "${MBD_TESTSUITE_INITIAL_VALUE_END}" "${MBD_TESTSUITE_CONTROL_DATA_END}"

                        simple_testsuite.sh --patch-input "yes" --prefix-output "${mbd_output_dir}" ${other_arguments} --exit-status-mask $((mbd_exit_status_mask)) >& "${simple_testsuite_log_file}"

                        rc=$?

                        test_status="FAILED"

                        if test ${rc} -eq 0; then
                            test_status="PASSED"
                        else
                            case ${rc} in
                                130)
                                    echo "Interrupted"
                                    exit 1;
                                    ;;
                                143)
                                    echo "Terminated"
                                    exit 1;
                                    ;;
                                137)
                                    echo "Killed"
                                    exit 1
                                    ;;
                            esac
                            failed_tests="${failed_tests} ${mbd_output_dir}"
                        fi

                        keep_output_flag="no"

                        case "${mbdyn_keep_output}" in
                            all)
                                keep_output_flag="yes"
                                ;;
                            failed)
                                case "${test_status}" in
                                    FAILED)
                                        keep_output_flag="yes"
                                        ;;
                                esac
                                ;;
                        esac

                        if test "${keep_output_flag}" = "no"; then
                            rm -f "${simple_testsuite_log_file}" "${MBD_TESTSUITE_INITIAL_VALUE_BEGIN}" "${MBD_TESTSUITE_INITIAL_VALUE_END}" "${MBD_TESTSUITE_CONTROL_DATA_BEGIN}" "${MBD_TESTSUITE_CONTROL_DATA_END}"
                        fi

                        printf 'TEST %s\n' "${test_status}"
                    done
                done
            done
        done
    done
done
if test -z "${failed_tests}"; then
    echo "All tests passed"
    exit 0
else
    printf '%d tests failed\n' "$(echo ${failed_tests} | wc -w)"
    for failed_test in ${failed_tests}; do
        printf "Test %s failed\n" "${failed_test}"
    done
    exit 1
fi
