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

if ! test -f "${program_dir}/octave_pkg_build_job.sh"; then
    program_dir=$(dirname $(realpath $(which "${program_name}")))
fi

if ! test -d "${program_dir}"; then
    echo "Failed to detect program directory"
    exit 1
fi

OCTAVE_EXEC="${OCTAVE_EXEC:-octave}"
OCT_PKG_LIST="${OCT_PKG_LIST:-nurbs:yes:master:no:unlimited netcdf:yes:master:no:unlimited mboct-octave-pkg:yes:master:no:unlimited mboct-numerical-pkg:yes:master:no:unlimited mboct-mbdyn-pkg:yes:master:no:unlimited mboct-fem-pkg:yes:master:no:unlimited}"
OCT_PKG_BUILD_DIR="${OCT_PKG_BUILD_DIR:-${program_dir}/var/cache/tmp/build/octave-pkg}"
OCT_PKG_BINARY_DIR="${OCT_PKG_BUILD_DIR}/binary-packages"
MBD_COMPILER_FLAGS="${MBD_COMPILER_FLAGS:--Ofast -Wall -march=native -mtune=native -Wno-unused-variable}"
MBD_CLEAN_ALL="${MBD_CLEAN_ALL:-no}"
MKL_INSTALL_PREFIX="${MKL_INSTALL_PREFIX:-/usr/lib}"
MKL_PKG_CONFIG="${MKL_PKG_CONFIG:-mkl-dynamic-lp64-gomp}"
NC_CONFIG="${NC_CONFIG:-nc-config}"

while ! test -z "$1"; do
    case "$1" in
        --mkl-install-prefix)
            MKL_INSTALL_PREFIX="$2"
            shift
            ;;
        --mkl-pkg-config)
            MKL_PKG_CONFIG="$2"
            shift
            ;;
        --netcdf-config)
            NC_CONFIG="$2"
            shift
            ;;
        --nlopt-install-prefix)
            NL_INSTALL_PREFIX="$2"
            shift
            ;;
        --compiler-flags)
            MBD_COMPILER_FLAGS="$2"
            shift
            ;;
        --linker-flags)
            export LDFLAGS="$2"
            shift
            ;;
        --octave-pkg-list)
            OCT_PKG_LIST="$2"
            shift
            ;;
        --octave-exec)
            OCTAVE_EXEC="$2"
            shift
            ;;
        --octave-pkg-prefix)
            OCT_PKG_INSTALL_PREFIX="$2"
            shift
            ;;
        --help|-h)
            printf "%s\n --octave-pkg-list <list-of-packages-and-flags>\n --octave-pkg-prefix <pkg-install-dir>\n --octave-exec <octave-executable>\n --help\n" "${program_name}"
            exit 1;
            ;;
        *)
            printf "%s: invalid argument \"%s\"\n" "${program_name}" "$1"
            exit 1
            ;;
    esac
    shift
done

if ! test -z "${NC_INSTALL_PREFIX}" && test -x "${NC_INSTALL_PREFIX}/bin/nc-config"; then
    NC_CONFIG="${NC_INSTALL_PREFIX}/bin/nc-config"
else
    NC_INSTALL_PREFIX=`${NC_CONFIG} --prefix`
fi

if ! "${OCTAVE_EXEC}" --version; then
    echo "Command ${OCTAVE_EXEC} --version failed"
    exit 1
fi

if ! test -z "${OCT_PKG_INSTALL_PREFIX}" && ! mkdir -p "${OCT_PKG_INSTALL_PREFIX}"; then
    echo "Failed to create directory ${OCT_PKG_INSTALL_PREFIX}"
    exit 1
fi

if ! mkdir -p "${OCT_PKG_BUILD_DIR}"; then
    echo "Failed to create directory ${OCT_PKG_BUILD_DIR}"
    exit 1
fi

if ! mkdir -p "${OCT_PKG_BINARY_DIR}"; then
    echo "Failed to create directory ${OCT_PKG_BINARY_DIR}"
    exit 1
fi

echo "octave packages installation job"
echo "Detecting MKL ..."

if test -d "${MKL_INSTALL_PREFIX}"; then
    echo "Search for ${MKL_PKG_CONFIG}.pc ..."

    MKL_PKG_CONFIG_FILE=`find "${MKL_INSTALL_PREFIX}" '(' -type f -and -name "${MKL_PKG_CONFIG}.pc" ')'`

    if test -f "${MKL_PKG_CONFIG_FILE}"; then
        MKL_PKG_CONFIG_PATH=`dirname "${MKL_PKG_CONFIG_FILE}"`
        echo "MKL_PKG_CONFIG_PATH=${MKL_PKG_CONFIG_PATH}"
    else
        echo "File ${MKL_PKG_CONFIG}.pc not found"
    fi
fi

if test -d "${MKL_PKG_CONFIG_PATH}"; then
    export PKG_CONFIG_PATH="${MKL_PKG_CONFIG_PATH}:${PKG_CONFIG_PATH}"
else
    echo "Warning: MKL_PKG_CONFIG_PATH does not exist"
fi

if pkg-config --validate "${MKL_PKG_CONFIG}"; then
    if test -z "${PARDISO_INC}"; then
        export PARDISO_INC=`pkg-config --cflags ${MKL_PKG_CONFIG}`
    fi
    if test -z "${PARDISO_LIBS}"; then
        export PARDISO_LIBS=`pkg-config --libs ${MKL_PKG_CONFIG}`
        export PARDISO_LIBS="${PARDISO_LIBS} `pkg-config --libs-only-L ${MKL_PKG_CONFIG} | sed  's/^-L\//-Wl,-rpath=\//g'`"
    fi
else
    echo "Warning: MKL pkg-config ${MKL_PKG_CONFIG} is not valid"
fi

NL_PKG_CONFIG_PATH=""

if ! test -z "${NL_INSTALL_PREFIX}"; then
    NL_PKG_CONFIG_PATH=`find ${NL_INSTALL_PREFIX} '(' -type d -and -name pkgconfig ')'`
fi

if ! test -z "${NL_PKG_CONFIG_PATH}"; then
    export PKG_CONFIG_PATH="${NL_PKG_CONFIG_PATH}:${PKG_CONFIG_PATH}"
fi

export NLOPT_LIBS=${NLOPT_LIBS:-`pkg-config --libs nlopt`}
export NLOPT_INC=${NLOPT_INC:-`pkg-config --cflags nlopt`}

NL_LIBDIR=`pkg-config --variable=libdir nlopt`

if ! test -z "${NL_LIBDIR}"; then
    NLOPT_LIBS="${NLOPT_LIBS} -Wl,-rpath=${NL_LIBDIR}"
fi

if ! test -z "${NC_INSTALL_PREFIX}"; then
    NC_PKG_CONFIG_PATH=`find ${NC_INSTALL_PREFIX} '(' -type d -and -name pkgconfig ')'`

    if ! test -z "${NC_PKG_CONFIG_PATH}"; then
        export PKG_CONFIG_PATH="${NC_PKG_CONFIG_PATH}:${PKG_CONFIG_PATH}"
    fi

    export PATH="${NC_INSTALL_PREFIX}/bin:${PATH}"
fi

echo PKG_CONFIG_PATH="${PKG_CONFIG_PATH}"
echo NLOPT_LIBS="${NLOPT_LIBS}"
echo NLOPT_INC="${NLOPT_INC}"
echo PARDISO_LIBS="${PARDISO_LIBS}"
echo PARDISO_INC="${PARDISO_INC}"

if ! cd "${OCT_PKG_BUILD_DIR}"; then
    "Cannot change directory to ${OCT_PKG_BUILD_DIR}"
    exit 1
fi

if ! test -z "${OCT_PKG_INSTALL_PREFIX}"; then
    OCTAVE_LOCAL_LIST=`printf 'pkg("local_list","%s");' "${OCT_PKG_INSTALL_PREFIX}/octave_packages"`
else
    OCTAVE_LOCAL_LIST=""
fi

if ! test -z "${OCT_PKG_INSTALL_PREFIX}"; then
    OCTAVE_PREFIX_CMD=`printf 'pkg("prefix","%s","%s");%s' "${OCT_PKG_INSTALL_PREFIX}" "${OCT_PKG_INSTALL_PREFIX}" "${OCTAVE_LOCAL_LIST}"`
else
    OCTAVE_PREFIX_CMD=""
fi

for pkgname_and_flags in ${OCT_PKG_LIST}; do
    pkgname=$(echo ${pkgname_and_flags} | awk -F ":" "{print \$1}")
    pkg_rebuild_flag=$(echo ${pkgname_and_flags} | awk -F ":" "{print \$2}")
    pkg_branch=$(echo ${pkgname_and_flags} | awk -F ":" "{print \$3}")
    printf "build package \"%s\" branch \"%s\"\n" "${pkgname}" "${pkg_branch}"

    case "${pkg_rebuild_flag}" in
        yes|no)
        ;;
        *)
            pkg_rebuild_flag="no"
            ;;
    esac

    if test -z "${pkg_branch}"; then
        pkg_branch="master"
    fi

    OCTAVE_CMD=`printf '%spkg("load","%s");' "${OCTAVE_LOCAL_LIST}" "${pkgname}"`

    if test "${pkg_rebuild_flag}" = "no" -a "${MBD_CLEAN_ALL}" = "no"; then
        echo ${OCTAVE_EXEC} -qfH --eval ${OCTAVE_CMD}

        if ${OCTAVE_EXEC} -qfH --eval "${OCTAVE_CMD}"; then
            printf "installation of \"%s\" will be skipped, because the package is already installed\n" "${pkgname}"
            continue
        else
            printf "package \"%s\" is not installed\n" "${pkgname}"
        fi
    fi

    OCTAVE_CMD=`printf '%spkg("uninstall","-nodeps","-local","-verbose","%s");' "${OCTAVE_LOCAL_LIST}" "${pkgname}"`

    echo ${OCTAVE_EXEC} -qfH --eval "${OCTAVE_CMD}"

    # Force a new installation
    if ! ${OCTAVE_EXEC} -qfH --eval "${OCTAVE_CMD}"; then
        printf "Warning: failed to uninstall package \"%s\"\n" "${pkgname}"
    fi

    case "${pkgname}" in
        mboct-*-pkg)
            if ! test -d "${pkgname}"; then
                if ! git clone -b "${pkg_branch}" -- "https://github.com/octave-user/${pkgname}.git" "${pkgname}"; then
                    printf "failed to clone package \"%s\" from github.com\n" "${pkgname}"
                    exit 1
                fi
            fi

            if ! test -d "${pkgname}"; then
                printf "package \"%s\" not found!\n" "${pkgname}"
                exit 1
            fi

            if ! cd "${pkgname}"; then
                exit 1
            fi

            if ! git fetch --all; then
                exit 1
            fi

            if ! git checkout "${pkg_branch}" --force; then
                exit 1
            fi

            if ! git pull --force; then
                exit 1
            fi

            if ! cd ..; then
                exit 1
            fi

            pkg_tar_file_binary=$(find "${OCT_PKG_BINARY_DIR}" -name "${pkgname}*.tar.gz" | tail -n1)

            if ! test -f "${pkg_tar_file_binary}" || test "${pkg_rebuild_flag}" == "yes"; then
                if ! make -C "${pkgname}" dist; then
                    printf "failed to build package \"%s\"\n" "${pkgname}"
                    exit 1
                fi

                if ! test -d "${pkgname}"; then
                    exit 1
                fi

                pkg_tar_file=$(find "${pkgname}" '(' -type f -and -path "*/target/${pkgname}*.tar.gz" ')')

                if ! test -f "${pkg_tar_file}"; then
                    printf "failed to create the package file \"%s\"\n" "${pkgname}"
                    exit 1
                fi

                OCTAVE_CMD=`printf '%spkg("build","-nodeps","-verbose","%s","%s");' "${OCTAVE_PREFIX_CMD}" "${OCT_PKG_BINARY_DIR}" "${pkg_tar_file}"`

                echo ${OCTAVE_EXEC} --eval "${OCTAVE_CMD}"

                CXXFLAGS="${MBD_COMPILER_FLAGS}" ${OCTAVE_EXEC} --eval "${OCTAVE_CMD}"

                if test $? != 0; then
                    printf "failed to build binary package \"%s\"\n" "${pkgname}"
                    exit 1
                fi

                pkg_tar_file_binary=`find "${OCT_PKG_BINARY_DIR}" -name "${pkgname}*.tar.gz" | tail -n1`

                if ! test -f "${pkg_tar_file_binary}"; then
                    printf "binary package \"%s\" not found\n" "${pkgname}"
                    exit 1
                fi
            fi

            OCTAVE_CMD=`printf '%spkg("install","-local","-verbose","%s");' "${OCTAVE_PREFIX_CMD}" "${pkg_tar_file_binary}"`

            echo "Installing binary package ..."
            echo ${OCTAVE_EXEC} --eval "${OCTAVE_CMD}"

            if ! ${OCTAVE_EXEC} --eval "${OCTAVE_CMD}"; then
                printf "failed to install binary package \"%s\"\n" "${pkg_tar_file_binary}"
                exit 1
            fi

            case "${pkgname}" in
                mboct-fem-pkg)
                    INSTALL_PREFIX_FEM_PRE_MESH_SIZE="${OCT_PKG_INSTALL_PREFIX:-/usr/local}"

                    ## FIXME: Need to install fem_pre_mesh_size since octave's package manager does not run "make install"
                    mkdir -p "${INSTALL_PREFIX_FEM_PRE_MESH_SIZE}/bin"

                    if ! cd "${pkgname}/src"; then
                        echo "Directory src not found"
                        exit 1
                    fi

                    if ! test -x configure; then
                        chmod u+x bootstrap
                        ./bootstrap
                    fi

                    if ! test -f Makefile; then
                        if ! ./configure CXXFLAGS="${MBD_COMPILER_FLAGS}" --prefix="${INSTALL_PREFIX_FEM_PRE_MESH_SIZE}"; then
                            echo "configure failed"
                            exit 1
                        fi
                    fi

                    ## No need to recompile the whole package!
                    if ! make fem_pre_mesh_size; then
                        echo "Failed to build fem_pre_mesh_size"
                        exit 1
                    fi

                    if ! install fem_pre_mesh_size "${INSTALL_PREFIX_FEM_PRE_MESH_SIZE}/bin"; then
                        echo "Failed to install fem_pre_mesh_size"
                        exit 1
                    fi
                    if ! cd ..; then
                        echo "Directory is not valid"
                        exit 1
                    fi
                    ;;
            esac
            ;;
        *)
            case "${pkgname}" in
                netcdf)
                    LDFLAGS_SAVE="${LDFLAGS}"
                    export LDFLAGS="-Wl,-rpath=`${NC_INSTALL_PREFIX}/bin/nc-config --libdir` ${LDFLAGS}"
                    ;;
            esac

            # Assume that the package is hosted at octave-forge
            OCTAVE_CMD=`printf '%spkg("install","-local","-verbose","-forge","%s");' "${OCTAVE_PREFIX_CMD}" "${pkgname}"`

            echo ${OCTAVE_EXEC} --eval "${OCTAVE_CMD}"

            if ! ${OCTAVE_EXEC} -qfH --eval "${OCTAVE_CMD}"; then
                printf "failed to install package \"%s\" from octave-forge\n" "${pkgname}"
                exit 1
            fi

            export LDFLAGS="${LDFLAGS_SAVE}"
            ;;
    esac

    OCTAVE_CMD=`printf '%spkg("list");pkg("load","%s");' "${OCTAVE_LOCAL_LIST}" "${pkgname}"`

    echo ${OCTAVE_EXEC} -qfH --eval "${OCTAVE_CMD}"

    if ${OCTAVE_EXEC} -qfH --eval "${OCTAVE_CMD}"; then
        printf "package \"%s\" was loaded successfully\n" "${pkgname}"
    else
        printf "failed to load package \"%s\"\n" "${pkgname}"
        exit 1
    fi
done
