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
    program_dir=$(dirname $(realpath $(which "${program_name}")))
fi

if test -f "${program_dir}/mbdyn_input_file_format.awk"; then
    export PATH=${program_dir}:${PATH}
fi

mbdyn_testsuite_prefix_output=""
mbdyn_keep_output="unexpected"
## FIXME: fourbar_int will fail with aztecoo and amesos
## mbdyn_linear_solvers="aztecoo amesos naive umfpack klu pardiso pardiso_64 y12 spqr qr lapack"
mbdyn_linear_solvers="naive umfpack klu pardiso pardiso_64 y12 spqr qr lapack"
mbdyn_matrix_handlers="map cc dir grad"
mbdyn_matrix_scale_methods="rowmaxcolumnmax iterative lapack rowmax columnmax rowsum columnsum"
mbdyn_matrix_scale_when="never always once"
mbdyn_nonlinear_solvers="newtonraphson linesearch linesearch-modified nox nox-newton-krylov nox-direct nox-broyden-linesearch nox-broyden-trust-region nox-broyden-inexact-trust-region mcpnewtonminfb mcpnewtonfb bfgs"
mbdyn_autodiff_options="autodiff noautodiff"
mbdyn_method="impliciteuler cranknicolson ms2,0.6 ms3,0.6 ms4,0.6 ss2,0.6 ss3,0.6 ss4,0.6 hope,0.6 Bathe,0.6 msstc3,0.6 msstc4,0.6 msstc5,0.6 mssth3,0.6 mssth4,0.6 mssth5,0.6 DIRK33 DIRK43 DIRK54 hybrid,ms,0.6"
mbdyn_output="netcdf-text"
mbdyn_abort_after="input assembly derivatives regularstep,2"
mbdyn_skip_initial_joint_assembly="not-skip skip"
mbdyn_initial_assembly_of_deformable_and_force_elements="exclude include"
declare -i mbd_exit_status_mask=0
other_arguments=""

while ! test -z "$1"; do
    case "$1" in
        --prefix-output)
            mbdyn_testsuite_prefix_output="$2"
            shift
            ;;
        --linear-solvers)
            mbdyn_linear_solvers="$2"
            shift
            ;;
        --matrix-handlers)
            mbdyn_matrix_handlers="$2"
            shift
            ;;
        --scale-methods)
            mbdyn_matrix_scale_methods="$2"
            shift
            ;;
        --scale-when)
            mbdyn_matrix_scale_when="$2"
            shift
            ;;
        --nonlinear-solvers)
            mbdyn_nonlinear_solvers="$2"
            shift
            ;;
        --output-format)
            mbdyn_output="$2"
            shift
            ;;
        --abort-after)
            mbdyn_abort_after="$2"
            shift
            ;;
        --skip-initial-joint-assembly)
            mbdyn_skip_initial_joint_assembly="$2"
            shift
            ;;
        --initial-assembly-of-deformable-and-force-elements)
            mbdyn_initial_assembly_of_deformable_and_force_elements="$2"
            shift
            ;;
        --patch-input)
            shift
            ;;
        --autodiff)
            mbdyn_autodiff_options="$2"
            shift
            ;;
        --method)
            mbdyn_method="$2"
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
        --help)
            printf "%s\n  --prefix-output <output_dir>\n" "${program_name}"
            printf "  --prefix-input <input_dir>\n"
            printf "  --linear-solvers \"{naive|umfpack|klu|pardiso|pardiso_64|y12|spqr|qr|lapack} {...}\"\n"
            printf "  --matrix-handlers \"{map|cc|dir|grad} {...}\"\n"
            printf "  --scale-methods \"{rowmaxcolumnmax|iterative|lapack|rowmax|columnmax|rowsum|columnsum} {...}\"\n"
            printf "  --scale-when \"{never|always|once} {...}\"\n"
            printf "  --nonlinear-solvers \"{newtonraphson|linesearch|linesearch-modified|nox|nox-newton-krylov|nox-direct|nox-broyden-linesearch|nox-broyden-trust-region|nox-broyden-inexact-trust-region|mcpnewtonminfb|mcpnewtonfb|bfgs} {...}\"\n"
            printf "  --autodiff {autodiff|noautodiff}\n"
            printf "  --method \"{impliciteuler|cranknicolson|ms2,0.6|ms3,0.6|ms4,0.6|ss2,0.6|ss3,0.6|ss4,0.6|hope,0.6|Bathe,0.6|msstc3,0.6|msstc4,0.6|msstc5,0.6|mssth3,0.6|mssth4,0.6|mssth5,0.6|DIRK33|DIRK43|DIRK54|hybrid,ms,0.6} {...}\"\n"
            printf "  --timeout <timeout_seconds>\n"
            printf "  --regex-filter <input_file_path_filter>\n"
            printf "  --exclude-inverse-dynamics {0|1}\n"
            printf "  --exclude-initial-value {0|1}\n"
            printf "  --threads <number_of_threads_per_task>\n"
            printf "  --tasks <number_of_tasks>\n"
            printf "  --verbose {yes|no}\n"
            printf "  --keep-output {all|failed|unexpected}\n"
            printf "  --mbdyn-args-add \"<arg1> <arg2> ... <argN>\"\n"
            printf "  --exec-gen {yes|no}\n"
            printf "  --exec-solver {yes|no}\n"
            printf "  --exec-status-mask <mask_errors_to_be_ignored>\n"
            printf "  --print-resources {no|all|time}\n"
            printf "  --suppressed-errors {syntax|element|feature|module|loadable|socked|interrupted}\n"
            printf "  --help\n"
            exit 1;
            ;;
        *)
            other_arguments="${other_arguments} $1 $2"
            shift
            ;;
    esac
    shift
done

((mbd_exit_status_mask|=0x1))

if ! test -d "${mbdyn_testsuite_prefix_output}"; then
    if ! mkdir -p "${mbdyn_testsuite_prefix_output}"; then
        exit 1
    fi
fi

simple_testsuite_log_file="${mbdyn_testsuite_prefix_output}/mbdyn-testsuite-patched.log"

simple_testsuite.sh --prefix-output "${mbdyn_testsuite_prefix_output}" ${other_arguments} --exec-solver no --exit-status-mask $((mbd_exit_status_mask)) >& "${simple_testsuite_log_file}"

rc=$?

if test "${mbdyn_keep_output}" = "no"; then
    rm -f "${simple_testsuite_log_file}"
fi

failed_tests=""

for mbd_linear_solver in ${mbdyn_linear_solvers}; do
    for mbd_mh_type in ${mbdyn_matrix_handlers}; do
        case "${mbd_linear_solver}" in
            naive|lapack|qr|aztecoo|amesos)
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
                ;;
        esac
        for mbd_mat_scale in ${mbdyn_matrix_scale_methods}; do
            for mbd_mat_scale_when in ${mbdyn_matrix_scale_when}; do
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
                    aztecoo|amesos)
                        case "${mbd_mat_scale_when}" in
                            never)
                            ;;
                            *)
                                continue
                                ;;
                        esac
                        ;;
                    pardiso|pardiso_64|qr|spqr|y12)
                        case "${mbd_mat_scale_when}" in
                            never)
                            ;;
                            *)
                                continue
                                ;;
                        esac
                esac
                for mbd_use_autodiff in ${mbdyn_autodiff_options}; do
                    case "${mbd_mh_type}" in
                        grad)
                            case "${mbd_use_autodiff}" in
                                noautodiff)
                                    continue
                                    ;;
                            esac
                            ;;
                    esac
                    for mbd_nonlin_solver in ${mbdyn_nonlinear_solvers}; do
                        case "${mbd_nonlin_solver}" in
                            mcpnewtonminfb)
                                case "${mbd_use_autodiff}" in
                                    noautodiff)
                                        continue
                                        ;;
                                esac
                                ;;
                            nox|nox-direct|nox-broyden*)
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
                                mbd_linear_solver_flags_pre=",colamd"
                                mbd_linear_solver_flags_post=""
                                ;;
                            amesos|aztecoo)
                                mbd_linear_solver_flags_pre=""
                                mbd_linear_solver_flags_post=", tolerance, 1e-8, max iterations, 100, preconditioner, klu,verbose,3"
                                ;;
                            umfpack|pardiso|pardiso_64)
                                mbd_linear_solver_flags_pre=""
                                mbd_linear_solver_flags_post=",max iterations, 10"
                                ;;
                            *)
                                mbd_linear_solver_flags_pre=""
                                mbd_linear_solver_flags_post=""
                                ;;
                        esac

                        case "${mbd_nonlin_solver}" in
                            nox)
                                mbd_nonlin_solver_flags="nox, minimum step, 1e-12, recovery step, 1e-4"
                                ;;
                            nox-newton-krylov)
                                mbd_nonlin_solver_flags="nox, modified, 10, jacobian operator, newton krylov, use preconditioner as solver, no, minimum step, 1e-12, recovery step, 1e-4"
                                ;;
                            nox-direct)
                                mbd_nonlin_solver_flags="nox, use preconditioner as solver, yes, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            nox-broyden-linesearch)
                                mbd_nonlin_solver_flags="nox, modified, 10, direction, broyden, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            nox-broyden-trust-region)
                                mbd_nonlin_solver_flags="nox, modified, 10, solver, trust region based, direction, broyden, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            nox-broyden-inexact-trust-region)
                                mbd_nonlin_solver_flags="nox, modified, 10, solver, inexact trust region based, direction, broyden, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            nox-broyden-tensor)
                                mbd_nonlin_solver_flags="nox, modified, 10, solver, tensor based, direction, broyden, minimum step, 1e-12, recovery step, 1e-12"
                                ;;
                            linesearch)
                                mbd_nonlin_solver_flags="linesearch, default solver options, heavy nonlinear, divergence check, no, lambda min, 1, print convergence info, yes, verbose, yes"
                                ;;
                            linesearch-heavy-nonlinear)
                                mbd_nonlin_solver_flags="linesearch, default solver options, heavy nonlinear, divergence check, no, lambda min, 1e-12, print convergence info, yes, verbose, yes"
                                ;;
                            linesearch-modified)
                                mbd_nonlin_solver_flags="linesearch, modified, 0, default solver options, heavy nonlinear, divergence check, no, lambda min, 1, print convergence info, yes, verbose, yes"
                                ;;
                            linesearch-modified-heavy-nonlinear)
                                mbd_nonlin_solver_flags="linesearch, modified, 0, default solver options, heavy nonlinear, divergence check, no, lambda min, 1e-12, print convergence info, yes, verbose, yes"
                                ;;
                            *)
                                mbd_nonlin_solver_flags="${mbd_nonlin_solver}"
                                ;;
                        esac

                        for mbd_method in ${mbdyn_method}; do
                            for mbd_output in ${mbdyn_output}; do
                                for mbd_abort_after in ${mbdyn_abort_after}; do
                                    for mbd_skip_initial_joint_assembly in ${mbdyn_skip_initial_joint_assembly}; do
                                        for mbd_initial_assembly_of_deformable_and_force_elements in ${mbdyn_initial_assembly_of_deformable_and_force_elements}; do
                                            case "${mbd_abort_after}" in
                                                derivatives)
                                                    case "${mbd_skip_initial_joint_assembly}" in
                                                        skip)
                                                            case "${mbd_linear_solver}" in
                                                                umfpack)
                                                                ;;
                                                                *)
                                                                    continue
                                                                    ;;
                                                            esac
                                                            case "${mbd_mh_type}" in
                                                                map)
                                                                ;;
                                                                *)
                                                                    continue
                                                                    ;;
                                                            esac
                                                            case "${mbd_mat_scale}" in
                                                                rowmaxcolumnmax)
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
                                                            case "${mbd_nonlin_solver}" in
                                                                newtonraphson)
                                                                ;;
                                                                *)
                                                                    continue
                                                                    ;;
                                                            esac
                                                            case "${mbd_method}" in
                                                                impliciteuler)
                                                                ;;
                                                                *)
                                                                    continue
                                                                    ;;
                                                            esac
                                                            ;;
                                                    esac
                                                    ;;
                                                *)
                                                    case "${mbd_skip_initial_joint_assembly}" in
                                                        skip)
                                                            continue
                                                            ;;
                                                    esac
                                                    ;;
                                            esac

                                            case "${mbd_initial_assembly_of_deformable_and_force_elements}" in
                                                include)
                                                    case "${mbd_abort_after}" in
                                                        assembly)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_skip_initial_joint_assembly}" in
                                                        skip)
                                                            continue
                                                            ;;
                                                    esac
                                                    ;;
                                                *)
                                                    ;;
                                            esac

                                            case "${mbd_abort_after}" in
                                                input)
                                                    case "${mbd_linear_solver}" in
                                                        umfpack)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_mh_type}" in
                                                        map)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_mat_scale}" in
                                                        rowmaxcolumnmax)
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
                                                    case "${mbd_nonlin_solver}" in
                                                        newtonraphson)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_method}" in
                                                        impliciteuler)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    ;;
                                                assembly)
                                                    case "${mbd_nonlin_solver}" in
                                                        newtonraphson)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_method}" in
                                                        impliciteuler)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    ;;
                                                derivatives)
                                                    case "${mbd_method}" in
                                                        impliciteuler)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    ;;
                                                regularstep,*)
                                                    case "${mbd_nonlin_solver}" in
                                                        newtonraphson)
                                                        ;;

                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_linear_solver}" in
                                                        umfpack)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_mh_type}" in
                                                        map)
                                                        ;;
                                                        *)
                                                            continue
                                                            ;;
                                                    esac
                                                    case "${mbd_mat_scale}" in
                                                        rowmaxcolumnmax)
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
                                                    ;;
                                            esac

                                            mbd_output_dir="${mbdyn_testsuite_prefix_output}/${mbd_linear_solver}/${mbd_mh_type}/${mbd_mat_scale}/${mbd_mat_scale_when}/${mbd_use_autodiff}/${mbd_nonlin_solver}/${mbd_method}/${mbd_output}/${mbd_abort_after}/${mbd_skip_initial_joint_assembly}/${mbd_initial_assembly_of_deformable_and_force_elements}"

                                            mkdir -p "${mbd_output_dir}"

                                            export MBD_TESTSUITE_INITIAL_VALUE_BEGIN="${mbd_output_dir}/mbd_init_val_begin.set"
                                            export MBD_TESTSUITE_INITIAL_VALUE_END="${mbd_output_dir}/mbd_init_val_end.set"
                                            export MBD_TESTSUITE_CONTROL_DATA_BEGIN="${mbd_output_dir}/mbd_control_data_begin.set"
                                            export MBD_TESTSUITE_CONTROL_DATA_END="${mbd_output_dir}/mbd_control_data_end.set"
                                            printf '    # mbd_init_val_begin.set currently not used!\n' > "${MBD_TESTSUITE_INITIAL_VALUE_BEGIN}"
                                            printf '    linear solver: %s,%s%s,scale,%s,%s%s;\n' "${mbd_linear_solver}" "${mbd_mh_type}" "${mbd_linear_solver_flags_pre}" "${mbd_mat_scale}" "${mbd_mat_scale_when}" "${mbd_linear_solver_flags_post}" > "${MBD_TESTSUITE_INITIAL_VALUE_END}"

                                            case ${mbd_abort_after} in
                                                none)
                                                ;;
                                                *)
                                                    printf '    abort after: %s;\n' "${mbd_abort_after}" >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            esac

                                            printf '    threads: disable;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            printf '    nonlinear solver: %s;\n' "${mbd_nonlin_solver_flags}" >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            printf '    output: iterations, cpu time, solver condition number, stat, yes;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            printf '    tolerance: 1e-4;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            printf '    derivatives tolerance: 1e-4;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            printf '    derivatives max iterations: 10;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            printf '    derivatives coefficient: auto;\n' >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"
                                            printf '    method: %s;\n' "${mbd_method}" >> "${MBD_TESTSUITE_INITIAL_VALUE_END}"


                                            case ${mbd_use_autodiff} in
                                                autodiff)
                                                    mbd_use_autodiff_cmd="use automatic differentiation;"
                                                    ;;
                                                *)
                                                    mbd_use_autodiff_cmd="# automatic differentiation disabled"
                                                    ;;
                                            esac

                                            printf '    print: all, to file;\n' > "${MBD_TESTSUITE_CONTROL_DATA_BEGIN}"
                                            printf '    %s\n' "${mbd_use_autodiff_cmd}" > "${MBD_TESTSUITE_CONTROL_DATA_END}"

                                            case "${mbd_initial_assembly_of_deformable_and_force_elements}" in
                                                include)
                                                    mbd_init_ass_prefix=""
                                                    ;;
                                                *)
                                                    mbd_init_ass_prefix="##"
                                                    ;;
                                            esac

                                            printf '    %s%s\n' "${mbd_init_ass_prefix}" "initial assembly of deformable and force elements;" >> "${MBD_TESTSUITE_CONTROL_DATA_END}"
                                            printf '    %s%s\n' "${mbd_init_ass_prefix}" "initial stiffness: 1e6, 1e6;" >> "${MBD_TESTSUITE_CONTROL_DATA_END}"
                                            printf '    %s%s\n' "${mbd_init_ass_prefix}" "max iterations: 10;" >> "${MBD_TESTSUITE_CONTROL_DATA_END}"

                                            case "${mbd_skip_initial_joint_assembly}" in
                                                skip)
                                                    mbd_skip_init_ass_cmd="skip initial joint assembly"
                                                    ;;
                                                *)
                                                    mbd_skip_init_ass_cmd="## skip initial joint assembly"
                                                    ;;
                                            esac

                                            printf '    %s;\n' "${mbd_skip_init_ass_cmd}" >> "${MBD_TESTSUITE_CONTROL_DATA_END}"

                                            case ${mbd_output} in
                                                netcdf)
                                                    printf '    output results: netcdf, no text;\n' >> "${MBD_TESTSUITE_CONTROL_DATA_END}"
                                                    ;;
                                                netcdf-text)
                                                    printf '    output results: netcdf, text;\n' >> "${MBD_TESTSUITE_CONTROL_DATA_END}"
                                                    ;;
                                                text)
                                                    printf '    # output results: text;\n' >> "${MBD_TESTSUITE_CONTROL_DATA_END}"
                                                    ;;
                                                none)
                                                    printf '    default output: none;\n' >> "${MBD_TESTSUITE_CONTROL_DATA_END}"
                                                    ;;
                                            esac

                                            simple_testsuite_log_file="${mbd_output_dir}/mbdyn-testsuite-patched.log"

                                            echo "${mbd_output_dir}" > "${simple_testsuite_log_file}"

                                            cat "${MBD_TESTSUITE_INITIAL_VALUE_BEGIN}" >> "${simple_testsuite_log_file}"
                                            cat "${MBD_TESTSUITE_INITIAL_VALUE_END}" >> "${simple_testsuite_log_file}"
                                            cat "${MBD_TESTSUITE_CONTROL_DATA_BEGIN}" >> "${simple_testsuite_log_file}"
                                            cat "${MBD_TESTSUITE_CONTROL_DATA_END}" >> "${simple_testsuite_log_file}"

                                            simple_testsuite.sh --exec-gen "no" --patch-input "yes" --prefix-output "${mbd_output_dir}" --exit-status-mask $((mbd_exit_status_mask)) ${other_arguments} 2>&1 >> "${simple_testsuite_log_file}"

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

                                            printf 'TEST \"%s\" %s\n' "${mbd_output_dir}" "${test_status}"
                                        done
                                    done
                                done
                            done
                        done
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
