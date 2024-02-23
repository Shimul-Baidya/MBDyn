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

# Purpose:
# Make sure that we are not starting a new pipeline until all remaining processes from a previous pipeline were terminated.
# Although this should be ensured by the gitlab runner, actually it is not always the case if a pipeline was cancelled.
# For a better solution see also https://gitlab.com/gitlab-org/gitlab-runner/-/issues/27443#note_1197449654

MBD_FORCE_KILL=${MBD_FORCE_KILL:-no}
program_name=$0

while ! test -z "$1"; do
    case "$1" in
        --force-kill)
            MBD_FORCE_KILL="$2"
            shift
            ;;
    esac
    shift
done

OUTPUT_FORMAT="pid,comm,tty,state,cputime,vsize,uname,pcpu,pmem,tty"
export AWKPATH=$(dirname ${program_name}):${AWKPATH}

function list_pid()
{
    ps -u $(whoami) -o "pid,comm,tty,state" | awk -f mbdyn_cleanup_filter.awk
}

function list_process()
{
    pid=$1

    ps -h -p ${pid} -o "${OUTPUT_FORMAT}"
}

success="yes"

declare -i killed_pids=0

pid_list=$(list_pid)

if test -z "${pid_list}"; then
    echo "There are no matching processes"
    exit 0
fi

echo "Matching processes for user $(whoami)"
for pid in ${pid_list}; do
    list_process ${pid}

    if test "${MBD_FORCE_KILL}" = "yes"; then
        echo "Process ${pid} will be killed"
        if ! kill -9 ${pid}; then
            echo "Failed to kill ${pid}"
            success="no"
        fi
        ((++killed_pids))
    else
        echo "Process ${pid} is still running: use MBD_FORCE_KILL=yes to kill this process!"
        success="no"
    fi
done

if test $((killed_pids)) -gt 0; then
    pid_list=$(list_pid)

    if ! test -z "${pid_list}"; then
        echo "Matching processes for user $(whoami)"
    fi

    for pid in $(list_pid); do
        list_process ${pid}
    done
fi

if test "${MBD_FORCE_KILL}" = "yes" && ! test "${success}" = "yes"; then
    exit 1;
fi

if test "${MBD_FORCE_KILL}" = "no" && ! test "${success}" = "yes"; then
    exit 1;
fi
