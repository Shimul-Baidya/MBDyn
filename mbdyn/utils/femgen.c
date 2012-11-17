/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "ac/f2c.h"

extern int __FC_DECL__(femgen)(char outname[73]);

static void
usage(FILE *outf, int rc)
{
	fprintf(outf,
"usage: femgen [-h] [[-f] <outfile>]\n"
"\t-h\t\tthis message\n"
"\t-f <outfile>\toutput file name (up to 72 characters long)\n"
		);
	exit(rc);

}

int
main(int argc, char *argv[])
{
	char outname[73] = { ' ' };

	for (;;) {
		int opt = getopt(argc, argv, "f:h");
		if (opt == -1) {
			break;
		}

		switch (opt) {
		case 'f': {
			size_t len = strlen(optarg);
			if (len >= sizeof(outname)) {
				fprintf(stderr, "femgen: output file name '%s' too long; trim to 72 bytes or less\n", optarg);
				exit(EXIT_FAILURE);
			}

			strncpy(outname, optarg, len + 1);
			} break;

		case '?':
		case 'h':
			usage(stdout, EXIT_SUCCESS);

		default:
			fprintf(stderr, "femgen: unhandled option '%c'\n", opt);
			usage(stderr, EXIT_FAILURE);
		}
	}

	if (optind < argc) {
		if (outname[0] != ' ') {
			fprintf(stderr, "femgen: output file name already set using '-f' option\n");
			exit(EXIT_FAILURE);

		} else {
			size_t len = strlen(argv[optind]);
			if (len >= sizeof(outname)) {
				fprintf(stderr, "femgen: output file name '%s' too long; trim to 72 bytes or less\n", argv[optind]);
				exit(EXIT_FAILURE);
			}

			strncpy(outname, argv[optind], len + 1);
			optind++;
		}
	}

	if (optind < argc) {
		fprintf(stderr, "femgen: extra args ignored\n");
	}

	return __FC_DECL__(femgen)(outname);
}
