/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

/* socket driver */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

/* include del programma */

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>

#include <rtai_out_elem.h>
#include <mbrtai_utils.h>
#include <dataman.h>

/* RTMBDynOutElem - begin */

RTMBDynOutElem::RTMBDynOutElem(unsigned int uL, const std::string& m,
	const std::string& host, unsigned long n, bool c, StreamContent *pSC,
	bool bNonBlocking)
: Elem(uL, flag(0)),
StreamOutElem(uL, m, 1),
host(host), node(n), create(c), port(-1),
bNonBlocking(bNonBlocking),
pSC(pSC),
mbx(0),
f_send(0)
{
	/* RATIONALE:
	 *
	 * if host/node is present, the mailbox is "remote";
	 * if it not, we may need to create it
	 */

	if (create) {
		ASSERT(node == 0);

		if (rtmbdyn_rt_mbx_init(name.c_str(), size, &mbx)) {
			silent_cerr("RTMBDyn mailbox(" << name << ") "
				"init failed" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		if (node) {
			/* get port ... */
			port = rtmbdyn_rt_request_port(node);
			/* FIXME: what in case of failure? */
		}

		if (rtmbdyn_RT_get_adr(node, port, name.c_str(), &mbx)) {
			silent_cerr("RTMBDyn mailbox(" << name << ") "
				"get_adr failed" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	if (bNonBlocking) {
		f_send = rtmbdyn_RT_mbx_send_if;
	} else {
		f_send = rtmbdyn_RT_mbx_send;
	}
}

RTMBDynOutElem::~RTMBDynOutElem(void)
{
	if (mbx) {
		rtmbdyn_rt_mbx_delete(&mbx);
	}

	if (pSC) {
		SAFEDELETE(pSC);
	}
}

std::ostream&
RTMBDynOutElem::Restart(std::ostream& out) const
{
	return out << "# RTMBDynOutElem(" << GetLabel() << ") "
		"not implemented yet" << std::endl;
}

void
RTMBDynOutElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{
	pSC->Prepare();

	int rc = f_send(node, -port, mbx, pSC->GetBuf(), pSC->GetSize());
	if (rc != pSC->GetSize()) {
		/* FIXME: error */
	}
}

void
RTMBDynOutElem::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
}

Elem *
ReadRTMBDynOutElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel, StreamContent::Type type)
{
	unsigned long node = 0;
	std::string host;
	std::string name;
	bool create = false;

	if (HP.IsKeyWord("name") || HP.IsKeyWord("stream" "name")) {
		const char *m = HP.GetStringWithDelims();
		if (m == NULL) {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"unable to read mailbox name "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (strlen(m) != 6) {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"illegal mailbox name \"" << m << "\" "
				"(must be exactly 6 chars) "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		name = m;

	} else {
		silent_cerr("RTMBDynOutElem(" << uLabel << "): "
			"missing mailbox name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (HP.IsKeyWord("create")) {
		if (HP.IsKeyWord("yes")) {
			create = true;
		} else if (HP.IsKeyWord("no")) {
			create = false;
		} else {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"\"create\" must be \"yes\" or \"no\" "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
	
	if (HP.IsKeyWord("local") || HP.IsKeyWord("path")) {
		const char *m = HP.GetStringWithDelims();
		silent_cout("RTMBDynOutElem(" << uLabel << "): "
			"local path \"" << m << "\" silently ignored "
			"at line " << HP.GetLineData() << std::endl);

	}

	if (HP.IsKeyWord("port")){
		int p = HP.GetInt();
		
		silent_cout ("RTMBDynOutElem(" << uLabel << "): "
			"port " << p << " silently ignored "
			"at line " << HP.GetLineData() << std::endl);
	}

	if (HP.IsKeyWord("host")) {
#if defined(HAVE_GETHOSTBYNAME) || defined(HAVE_INET_ATON)
		const char *h;
		
		h = HP.GetStringWithDelims();
		if (h == NULL) {
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"unable to read host "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (create) {
			silent_cout("RTMBDynOutElem(" << uLabel << "): "
				"host name \"" << h << "\" silently ignored "
				"at line " << HP.GetLineData() << std::endl);			
		} else {

			host = h;

			/* resolve host
		 	* FIXME: non-reentrant ... */
#if defined(HAVE_GETHOSTBYNAME)
			struct hostent *he = gethostbyname(host.c_str());
			if (he != NULL)
			{
				node = ((unsigned long *)he->h_addr_list[0])[0];
			} 
#elif defined(HAVE_INET_ATON)
			struct in_addr addr;
			if (inet_aton(host.c_str(), &addr)) {
				node = addr.s_addr;
			}
#endif /* ! HAVE_GETHOSTBYNAME && HAVE_INET_ATON */
			else {
				silent_cerr("RTMBDynOutElem(" << uLabel << "): "
					"unable to convert host "
					"\"" << host << "\" "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
#else /* ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON */
			silent_cerr("RTMBDynOutElem(" << uLabel << "): "
				"host (RTAI RPC) not supported "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
#endif /* ! HAVE_GETHOSTBYNAME && ! HAVE_INET_ATON */
		}
	}

	bool bNonBlocking = true;
	while (HP.IsArg()) {
		if (HP.IsKeyWord("blocking")) {
			bNonBlocking = false;

		} else if (HP.IsKeyWord("non" "blocking")) {
			bNonBlocking = true;

		} else {
			break;
		}
	}

	StreamContent *pSC = ReadStreamContent(pDM, HP, type);

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("RTMBDynOutElem(" << uLabel << "): "
			"semicolon expected at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
      
	/* costruzione dell'elemento */
	Elem *pEl = NULL;
	SAFENEWWITHCONSTRUCTOR(pEl, RTMBDynOutElem,
			RTMBDynOutElem(uLabel,
				host, name, node, create, pSC, bNonBlocking));
	return pEl;
}

