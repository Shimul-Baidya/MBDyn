/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2002
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

#ifndef FRICTION_H
#define FRICTION_H

doublereal 
norm(const doublereal &v)
{
	return fabs(v);
}

doublereal
dir(const doublereal &v)
{
	return copysign(1., v);
}

doublereal
norm(const Vec3& v)
{
	return v.Norm();
}

Vec3
dir(const Vec3& v)
{
	doublereal d = v.Norm();
	if (d < DBL_EPSILON) {
		return Zero3;
	}
	return v/d;
}

Vec3
tanh(const Vec3& v)
{
	doublereal d = v.Norm();
	if (d < DBL_EPSILON) {
		return Zero3;
	}
	return v*(tanh(d)/d);
}

template <class T>
class Friction {
public:
	enum UpdateType { FIRST, ANY };

protected:
	T m_F;

public:
	virtual ~Friction(void) { NO_OP; };

	virtual void Update(const T &force, 
			const T &position, 
			const T &velocity, 
			Friction::UpdateType = ANY) = 0;

	virtual const T &F(void) const {
		return m_F;
	};
};

template <class T>
class StepFriction : public Friction<T> {
private:
	doublereal m_mu0;
public:
	StepFriction(const doublereal& mu0) 
	: m_mu0(mu0) {
		NO_OP;
	};
	
	virtual ~StepFriction(void) { NO_OP; };

	virtual void Update(const T &force,
			const T & /* position */, 
			const T & velocity, 
			Friction::UpdateType = ANY) {
		m_F = dir(velocity)*(-norm(force)*m_mu0);
	};
};

template <class T>
class TanhFriction : public Friction<T> {
private:
	doublereal m_mu0;
	doublereal m_vRef;
public:
	TanhFriction(const doublereal& mu0, const doublereal vRef) 
	: m_mu0(mu0), m_vRef(vRef) {
		NO_OP;
	};
	
	virtual ~TanhFriction(void) { NO_OP; };

	virtual void Update(const T &force,
			const T & /* position */, 
			const T &velocity, 
			Friction::UpdateType = ANY) {
#if 0
		m_F = dir(velocity)*(-norm(force)*m_mu0*tanh(norm(velocity)/m_vRef));
#else
		m_F = tanh(velocity/m_vRef)*(-norm(force)*m_mu0);
#endif
	};
};

#endif /* FRICTION_H */

