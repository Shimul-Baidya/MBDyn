! $Header$
! MBDyn (C) is a multibody analysis code. 
! http://www.mbdyn.org
! 
! Copyright (C) 1996-2017
! 
! Pierangelo Masarati	<masarati@aero.polimi.it>
! Paolo Mantegazza	<mantegazza@aero.polimi.it>
! 
! Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
! via La Masa, 34 - 20156 Milano, Italy
! http://www.aero.polimi.it
! 
! Changing this copyright notice is forbidden.
! 
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation (version 2 of the License).
! 
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
! AUTHOR: Reinhard Resch <R.RESCH@secop.com>
!        Copyright (C) 2011(-2017) all rights reserved.
!
!        The copyright of this code is transferred
!        to Pierangelo Masarati and Paolo Mantegazza
!        for use in the software MBDyn as described
!        in the GNU Public License version 2.1

!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
! 

MODULE HYDRODYNAMIC_PLAIN_BEARING_DV
  USE ISO_C_BINDING
  USE DIFFSIZES
  USE HYDRODYNAMIC_PLAIN_BEARING
  IMPLICIT NONE
CONTAINS
!  Differentiation of hydrodynamic_plain_bearing_force in forward (tangent) mode:
!   variations   of useful results: k
!   with respect to varying inputs: e omega e_dot
!   RW status of diff variables: e:in omega:in k:out e_dot:in
  SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE_DV(bdat, omega, omegad, e&
&   , ed, e_dot, e_dotd, k, kd, eps, eps_dot, delta, sod, sov, mu, beta&
&   , nbdirs) BIND(C)
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
!---------------
!       PRINT OUT OF RESULTS
!--------------------------------------------------------------------------------------------------------------------------------
!---------------
!~ PRINT *,'HYDRODYNAMIC_PLAIN_BEARING_FORCE'
!~ PRINT *,'b=',bdat%b
!~ PRINT *,'d=',bdat%d
!~ PRINT *,'eta=',eta
!~ PRINT *,'Psi=',Psi
!~ PRINT *,'e=',e
!~ PRINT *,'e_dot=',e_dot
!~ PRINT *,'omega_w=',omega_w
!~ PRINT *,'delta=',delta
!~ PRINT *,'phi=',phi
!~ PRINT *,'kappa=',kappa
!~ PRINT *,'abs_e=',abs_e
!~ PRINT *,'abs_e_dot=',abs_e_dot
!~ PRINT *,'eps_dot=',eps_dot
!~ PRINT *,'eps=',eps
!~ PRINT *,'delta_dot=',delta_dot
!~ PRINT *,'a(1)=',bdat%a(1)
!~ PRINT *,'a(2)=',bdat%a(2)
!~ PRINT *,'a(3)=',bdat%a(3)
!~ PRINT *,'a(4)=',bdat%a(4)
!~ PRINT *,'a(5)=',bdat%a(5)
!~ PRINT *,'a(6)=',bdat%a(6)
!~ PRINT *,'a(7)=',bdat%a(7)
!~ PRINT *,'a(8)=',bdat%a(8)
!~ PRINT *,'a(9)=',bdat%a(9)
!~ PRINT *,'SoD=',SoD
!~ PRINT *,'SoV=',SoV
!~ PRINT *,'beta=',beta
!~ PRINT *,'mu=',mu
!~ PRINT *,'alpha=',alpha
!~ PRINT *,'omega_res=',omega_res
!~ PRINT *,'abs_FD=',abs_FD
!~ PRINT *,'abs_FV=',abs_FV
!~ PRINT *,'abs_MR=',abs_MR
!~ PRINT *, 'k=',k
!-------------------------------------------------------------------------------------------------------------------------
!       hydrodynamic plain bearing calculation according to Butenschoen's theory
!-------------------------------------------------------------------------------------------------------------------------
!       COORDINATE SYSTEM:
!       x ... axial direction
!       y, z ... radial direction
!-------------------------------------------------------------------------------------------------------------------------
!        INPUT PARAMETERS
!-------------------------------------------------------------------------------------------------------------------------
!       b   ... bearing width [m]
!       d   ... shaft diameter [m]
!       Psi ... relative radial clearance Psi = ( D - d ) / D [1]
!       eta ... dynamic oil viscosity [Pa*s]
!       omega(1) ... angular velocity of the shaft [rad/s]
!	    omega(2) ... angular velocity of the bearing [rad/s]
!       e ...  radial eccentricity of the shaft 
!               e(1) = ey  [m]
!               e(2) = ez  [m]
!       e_dot(1) ... velocity of the shaft relative to the bearing 
!               e_dot(1) = ey_dot [m/s]
!               e_dot(2) = ez_dot [m/s]
!
!--------------------------------------------------------------------------------------------------------------------------
!       OUTPUT PARAMETERS
!--------------------------------------------------------------------------------------------------------------------------
!       k ... force on the bearing
!               k(1) = Fy [N]
!               k(2) = Fz [N]
!               k(3) = Mx [Nm]
!       eps ...             relative eccentricity of the shaft [1]
!       eps_dot ...      time derivative of the relative eccentricity [1/s]
!       delta ...          angle of minimum clearance between shaft and bearing [rad]
!       SoD ....           Sommerfeld number for rotation [1]
!       SoV ...            Sommerfeld number for displacement [1]
!       mu ...             friction coefficient [N/N]
!       beta  ...          angle between force for rotation and minimum clearance [rad]
    TYPE(BEARING_DATA), INTENT(IN) :: bdat
    DOUBLE PRECISION, INTENT(IN) :: omega(2), e(2), e_dot(2)
    DOUBLE PRECISION, INTENT(IN) :: omegad(nbdirsmax, 2), ed(nbdirsmax, &
&   2), e_dotd(nbdirsmax, 2)
    DOUBLE PRECISION, INTENT(OUT) :: k(3), eps, eps_dot, delta, sod, sov&
&   , mu, beta
  DOUBLE PRECISION, INTENT(OUT) :: kd(nbdirsmax, 3)
    DOUBLE PRECISION :: abs_e, abs_e_dot
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: abs_ed, abs_e_dotd
  DOUBLE PRECISION :: delta_dot, alpha, kappa, phi
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: delta_dotd, alphad, kappad&
&   , phid
  DOUBLE PRECISION :: omega_res, abs_fd, abs_fv, abs_mr
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: omega_resd, abs_fdd, &
&   abs_fvd, abs_mrd
    INTRINSIC ATAN2, SQRT, SUM, COS, SIN, ABS, SIGN
  INTRINSIC HUGE
  DOUBLE PRECISION, DIMENSION(2) :: arg1
  DOUBLE PRECISION, DIMENSION(nbdirsmax, 2) :: arg1d
  DOUBLE PRECISION :: arg2
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg2d
  DOUBLE PRECISION :: result1
    DOUBLE PRECISION :: arg10
  INTEGER :: nd
  INTEGER :: nbdirs
  DOUBLE PRECISION :: eps_dotd(nbdirsmax)
  DOUBLE PRECISION :: sovd(nbdirsmax)
    DOUBLE PRECISION :: abs1d(nbdirsmax)
  DOUBLE PRECISION :: epsd(nbdirsmax)
    DOUBLE PRECISION :: mud(nbdirsmax)
    DOUBLE PRECISION :: abs0d(nbdirsmax)
  DOUBLE PRECISION :: deltad(nbdirsmax)
  DOUBLE PRECISION :: sodd(nbdirsmax)
  DOUBLE PRECISION :: betad(nbdirsmax)
    DOUBLE PRECISION :: abs1
    DOUBLE PRECISION :: abs0
  arg1(:) = e(:)**2
  arg2 = SUM(arg1(:))
  DO nd=1,nbdirs
! angle of the position with minimum clearance between shaft and bearing
    deltad(nd) = (ed(nd, 2)*e(1)-ed(nd, 1)*e(2))/(e(2)**2+e(1)**2)
! angle of the velocity vector of the shaft relative to the bearing
      phid(nd) = (e_dotd(nd, 2)*e_dot(1)-e_dotd(nd, 1)*e_dot(2))/(e_dot(&
&       2)**2+e_dot(1)**2)
! angle between velocity vector and minimum clearance
    kappad(nd) = phid(nd) - deltad(nd)
! absolute value of the eccentricity of the shaft inside the bearing
    arg1d(nd, :) = 2*e(:)*ed(nd, :)
    arg2d(nd) = SUM(arg1d(nd, :))
    IF (arg2 .EQ. 0.0) THEN
      abs_ed(nd) = 0.D0
    ELSE
        abs_ed(nd) = arg2d(nd)/(2.0*SQRT(arg2))
    END IF
! absolute value of the velocity of the shaft relative to the bearing
    arg1d(nd, :) = 2*e_dot(:)*e_dotd(nd, :)
    arg2d(nd) = SUM(arg1d(nd, :))
! relative eccentricity of the shaft
      epsd(nd) = 2d0*abs_ed(nd)/bdat%s
  END DO
    delta = ATAN2(e(2), e(1))
    phi = ATAN2(e_dot(2), e_dot(1))
  kappa = phi - delta
    abs_e = SQRT(arg2)
  arg1(:) = e_dot(:)**2
  arg2 = SUM(arg1(:))
    abs_e_dot = SQRT(arg2)
  DO nd=1,nbdirs
    IF (arg2 .EQ. 0.0) THEN
      abs_e_dotd(nd) = 0.D0
    ELSE
        abs_e_dotd(nd) = arg2d(nd)/(2.0*SQRT(arg2))
    END IF
! time derivative of the relative eccentricity of the shaft
      eps_dotd(nd) = 2d0*(COS(kappa)*abs_e_dotd(nd)-kappad(nd)*SIN(kappa&
&       )*abs_e_dot)/bdat%s
  END DO
    eps_dot = 2d0*COS(kappa)*abs_e_dot/bdat%s
    eps = 2d0*abs_e/bdat%s
  IF (eps_dot .NE. 0d0) THEN
    DO nd=1,nbdirs
! eps is positive if it's time derivative is positive too
!       attention the signum function is zero if eps_dot is zero
!       but eps must not be zero in this case
        epsd(nd) = SIGN(1d0, eps_dot)*epsd(nd)
    END DO
      eps = SIGN(1d0, eps_dot)*eps
  END IF
! time derivative of angle of minimum clearance
  IF (abs_e .EQ. 0d0) THEN
! avoid division by zero
    delta_dot = 0d0
    DO nd=1,nbdirs
      delta_dotd(nd) = 0.D0
    END DO
  ELSE
    DO nd=1,nbdirs
        delta_dotd(nd) = ((ed(nd, 1)*e_dot(2)+e(1)*e_dotd(nd, 2)-ed(nd, &
&         2)*e_dot(1)-e(2)*e_dotd(nd, 1))*(e(2)**2+e(1)**2)-(e(1)*e_dot(&
&         2)-e(2)*e_dot(1))*(2*e(2)*ed(nd, 2)+2*e(1)*ed(nd, 1)))/(e(2)**&
&         2+e(1)**2)**2
    END DO
    delta_dot = (e(1)*e_dot(2)-e(2)*e_dot(1))/(e(2)**2+e(1)**2)
  END IF
    CALL SOMMERFELD_NUMBERS_EXT_DV(bdat, eps, epsd, omega, omegad, &
&                            delta_dot, delta_dotd, sod, sodd, sov, sovd&
&                            , beta, betad, mu, mud, nbdirs)
    omega_res = omega(1) + omega(2) - 2d0*delta_dot
    DO nd=1,nbdirs
! effective hydrodynamic angular velocity according to Butenschoen 1976
      omega_resd(nd) = omegad(nd, 1) + omegad(nd, 2) - 2d0*delta_dotd(nd&
&       )
! angle of the force for rotation
      alphad(nd) = deltad(nd) - SIGN(1d0, omega_res)*betad(nd)
    END DO
    alpha = delta - beta*SIGN(1d0, omega_res)
    IF (omega_res .GE. 0.) THEN
      DO nd=1,nbdirs
        abs0d(nd) = omega_resd(nd)
      END DO
      abs0 = omega_res
    ELSE
      DO nd=1,nbdirs
        abs0d(nd) = -omega_resd(nd)
      END DO
      abs0 = -omega_res
    END IF
    DO nd=1,nbdirs
! absolute value of the force for rotation
      abs_fdd(nd) = bdat%b*bdat%d*bdat%eta*(sodd(nd)*abs0+sod*abs0d(nd))&
&       /bdat%psi**2
! absolute value of the force for displacement
      abs_fvd(nd) = bdat%b*bdat%d*bdat%eta*(sovd(nd)*eps_dot+sov*&
&       eps_dotd(nd))/bdat%psi**2
    END DO
    abs_fd = sod*(bdat%b*bdat%d*bdat%eta*abs0)/bdat%psi**2
    abs_fv = sov*(bdat%b*bdat%d*bdat%eta*eps_dot)/bdat%psi**2
    result1 = HUGE(1d0)
    IF (mu .GE. result1) THEN
      IF (omega(1) - omega(2) .GE. 0.) THEN
        DO nd=1,nbdirs
          abs1d(nd) = omegad(nd, 1) - omegad(nd, 2)
        END DO
        abs1 = omega(1) - omega(2)
      ELSE
        DO nd=1,nbdirs
          abs1d(nd) = -(omegad(nd, 1)-omegad(nd, 2))
        END DO
        abs1 = -(omega(1)-omega(2))
      END IF
      DO nd=1,nbdirs
        abs_mrd(nd) = pi*bdat%b*bdat%d**2*bdat%eta*abs1d(nd)/bdat%psi/&
&         2d0
      END DO
      abs_mr = pi*bdat%b*bdat%d**2*bdat%eta*abs1/bdat%psi/2d0
    ELSE
      DO nd=1,nbdirs
        abs_mrd(nd) = bdat%d*(mud(nd)*abs_fd+mu*abs_fdd(nd))/2d0
      END DO
      abs_mr = mu*abs_fd*bdat%d/2d0
    END IF
! friction torque
    arg10 = omega(1) - omega(2)
    DO nd=1,nbdirs
! sum of force for rotation and force for displacement
      kd(nd, :) = 0.D0
      kd(nd, 1) = abs_fdd(nd)*COS(alpha) - abs_fd*alphad(nd)*SIN(alpha) &
&       + abs_fvd(nd)*COS(delta) - abs_fv*deltad(nd)*SIN(delta)
      kd(nd, 2) = abs_fdd(nd)*SIN(alpha) + abs_fd*alphad(nd)*COS(alpha) &
&       + abs_fvd(nd)*SIN(delta) + abs_fv*deltad(nd)*COS(delta)
      kd(nd, 3) = SIGN(1d0, arg10)*abs_mrd(nd)
    END DO
    k(1) = abs_fd*COS(alpha) + abs_fv*COS(delta)
    k(2) = abs_fd*SIN(alpha) + abs_fv*SIN(delta)
    k(3) = abs_mr*SIGN(1d0, arg10)
  END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE_DV
!  Differentiation of sommerfeld_numbers_ext in forward (tangent) mode:
!   variations   of useful results: sod beta sov mu
!   with respect to varying inputs: omega delta_dot eps
  SUBROUTINE SOMMERFELD_NUMBERS_EXT_DV(bdat, eps, epsd0, omega, omegad, &
&   delta_dot, delta_dotd, sod, sodd0, sov, sovd0, beta, betad, mu, mud&
&   , nbdirs)
    USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
    IMPLICIT NONE
    TYPE(BEARING_DATA), INTENT(IN) :: bdat
    DOUBLE PRECISION, INTENT(IN) :: eps, delta_dot, omega(2)
    DOUBLE PRECISION, INTENT(IN) :: epsd0(nbdirsmax), delta_dotd(&
&   nbdirsmax), omegad(nbdirsmax, 2)
    DOUBLE PRECISION, INTENT(OUT) :: sod, sov, beta, mu
    DOUBLE PRECISION, DIMENSION(nbdirsmax), INTENT(OUT) :: sodd0, sovd0&
&   , betad, mud
    DOUBLE PRECISION :: epsd, sodd, sovd, eps_max
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: epsdd, soddd, sovdd
    INTRINSIC ABS
    INTRINSIC SIGN
    INTEGER :: nd
    INTEGER :: nbdirs
    DOUBLE PRECISION :: abs0
    IF (eps .GE. 0.) THEN
      abs0 = eps
    ELSE
      abs0 = -eps
    END IF
    IF (abs0 .LT. bdat%eps_max) THEN
! According to the thesis of Butenschoen those approximations are valid until eps = 0.999
      CALL SOMMERFELD_NUMBERS_DV(bdat, eps, epsd0, omega, omegad, &
&                          delta_dot, delta_dotd, sod, sodd0, sov, sovd0&
&                          , beta, betad, mu, mud, nbdirs)
    ELSE
! Do a linear extrapolation above eps_max
      eps_max = SIGN(bdat%eps_max, eps)
      DO nd=1,nbdirs
        epsdd(nd) = epsd0(nd)
      END DO
      epsd = eps - eps_max
      CALL SOMMERFELD_NUMBERS_D_DV(bdat, eps_max, epsd, epsdd, omega, &
&                            omegad, delta_dot, delta_dotd, sod, sodd, &
&                            soddd, sov, sovd, sovdd, beta, mu, mud, &
&                            nbdirs)
      DO nd=1,nbdirs
        sodd0(nd) = soddd(nd)
        sovd0(nd) = sovdd(nd)
      END DO
      sod = sod + sodd
      sov = sov + sovd
      DO nd=1,nbdirs
        betad(nd) = 0.D0
      END DO
    END IF
  END SUBROUTINE SOMMERFELD_NUMBERS_EXT_DV
!  Differentiation of sommerfeld_numbers in forward (tangent) mode:
!   variations   of useful results: sod beta sov mu
!   with respect to varying inputs: omega delta_dot eps
  SUBROUTINE SOMMERFELD_NUMBERS_DV(bdat, eps, epsd, omega, omegad, &
&   delta_dot, delta_dotd, sod, sodd, sov, sovd, beta, betad, mu, mud, &
&   nbdirs)
    USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
    IMPLICIT NONE
    TYPE(BEARING_DATA), INTENT(IN) :: bdat
    DOUBLE PRECISION, INTENT(IN) :: eps, delta_dot, omega(2)
    DOUBLE PRECISION, INTENT(IN) :: epsd(nbdirsmax), delta_dotd(&
&   nbdirsmax), omegad(nbdirsmax, 2)
    DOUBLE PRECISION, INTENT(OUT) :: sod, sov, beta, mu
    DOUBLE PRECISION, DIMENSION(nbdirsmax), INTENT(OUT) :: sodd, sovd, &
&   betad, mud
    DOUBLE PRECISION, PARAMETER :: eps_min=1d-6
    INTRINSIC ABS
    INTRINSIC SQRT
    INTRINSIC ACOS
    INTRINSIC ATAN2
    INTRINSIC HUGE
    INTRINSIC SIN
    DOUBLE PRECISION :: arg1
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg1d
    DOUBLE PRECISION :: result1
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: result1d
    DOUBLE PRECISION :: pwx1
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: pwx1d
    DOUBLE PRECISION :: pwr1
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: pwr1d
    DOUBLE PRECISION :: result2
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: result2d
    DOUBLE PRECISION :: arg2
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg2d
    DOUBLE PRECISION :: arg3
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg3d
    INTEGER :: nd
    INTEGER :: nbdirs
    DOUBLE PRECISION :: abs1d(nbdirsmax)
    DOUBLE PRECISION :: abs4d(nbdirsmax)
    DOUBLE PRECISION :: abs7d(nbdirsmax)
    DOUBLE PRECISION :: abs0d(nbdirsmax)
    DOUBLE PRECISION :: abs3d(nbdirsmax)
    DOUBLE PRECISION :: abs6d(nbdirsmax)
    DOUBLE PRECISION :: abs8
    DOUBLE PRECISION :: abs7
    DOUBLE PRECISION :: abs6
    DOUBLE PRECISION :: abs5
    DOUBLE PRECISION :: abs4
    DOUBLE PRECISION :: abs3
    DOUBLE PRECISION :: abs2
    DOUBLE PRECISION :: abs1
    DOUBLE PRECISION :: abs0
    DOUBLE PRECISION :: abs5d(nbdirsmax)
    DOUBLE PRECISION :: abs8d(nbdirsmax)
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
        abs0d(nd) = epsd(nd)
    END DO
      abs0 = eps
  ELSE
    DO nd=1,nbdirs
        abs0d(nd) = -epsd(nd)
    END DO
      abs0 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
        abs4d(nd) = epsd(nd)
    END DO
      abs4 = eps
  ELSE
    DO nd=1,nbdirs
        abs4d(nd) = -epsd(nd)
    END DO
      abs4 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
        abs7d(nd) = epsd(nd)
    END DO
      abs7 = eps
  ELSE
    DO nd=1,nbdirs
        abs7d(nd) = -epsd(nd)
    END DO
      abs7 = -eps
  END IF
    arg1 = pi**2*(1d0-eps**2) + 16d0*eps**2
    result1 = SQRT(arg1)
  pwx1 = 1d0 - eps**2
  DO nd=1,nbdirs
! Sommerfeld number for rotation according to Butenschoen 1976
      arg1d(nd) = 16d0*2*eps*epsd(nd) - pi**2*2*eps*epsd(nd)
      IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.D0
    ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
      sodd(nd) = (bdat%a(1)*((bdat%b**2*abs0d(nd)*2d0*(1d0-eps**2)**2/&
&       bdat%d**2+bdat%b**2*abs0*2d0*2**2*(1d0-eps**2)*eps*epsd(nd)/bdat&
&       %d**2)*result1*(abs4-1d0)/(2d0**2*(1d0-eps**2)**4)+bdat%b**2*&
&       abs0*(result1d(nd)*(abs4-1d0)+result1*abs4d(nd))/(bdat%d**2*2d0*&
&       (1d0-eps**2)**2))*(bdat%a(2)+abs7)-bdat%b**2*abs0*result1*bdat%a&
&       (1)*(abs4-1d0)*abs7d(nd)/(bdat%d**2*2d0*(1d0-eps**2)**2))/(bdat%&
&       a(2)+abs7)**2
! Sommerfeld number for displacement according to Butenschoen 1976
    pwx1d(nd) = -(2*eps*epsd(nd))
    IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. -(5d0/2d0) .EQ. INT(-(&
&       5d0/2d0)))) THEN
      pwr1d(nd) = -(5d0*pwx1**((-1)-5d0/2d0)*pwx1d(nd)/2d0)
    ELSE IF (pwx1 .EQ. 0.0 .AND. -(5d0/2d0) .EQ. 1.0) THEN
      pwr1d(nd) = pwx1d(nd)
    ELSE
      pwr1d(nd) = 0.0
    END IF
    IF (eps .EQ. 1.0 .OR. eps .EQ. (-1.0)) THEN
      result1d(nd) = 0.D0
    ELSE
        result1d(nd) = -(epsd(nd)/SQRT(1.0-eps**2))
    END IF
      arg1d(nd) = -(2*eps*epsd(nd))
  END DO
    sod = (bdat%b/bdat%d)**2*abs0/(2d0*(1d0-eps**2)**2)*result1*bdat%a(1&
&     )*(abs4-1d0)/(bdat%a(2)+abs7)
  pwr1 = pwx1**(-(5d0/2d0))
    result1 = ACOS(eps)
    arg1 = 1d0 - eps**2
    result2 = SQRT(arg1)
  DO nd=1,nbdirs
      IF (arg1 .EQ. 0.0) THEN
      result2d(nd) = 0.D0
    ELSE
        result2d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
      sovd(nd) = (4d0*bdat%b**2*bdat%a(8)*((pwr1d(nd)*(1d0-eps)-pwr1*&
&       epsd(nd))*((pi/2d0-1d0/2d0*result1)*(1d0+2d0*eps**2)+3d0/2d0*eps&
&       *result2)+pwr1*(1d0-eps)*((pi/2d0-1d0/2d0*result1)*2d0*2*eps*&
&       epsd(nd)-result1d(nd)*(1d0+2d0*eps**2)/2d0+3d0*(epsd(nd)*result2&
&       +eps*result2d(nd))/2d0))*(-bdat%a(9)-eps)/bdat%d**2+4d0*bdat%b**&
&       2*pwr1*((pi/2d0-1d0/2d0*result1)*(1d0+2d0*eps**2)+3d0/2d0*eps*&
&       result2)*bdat%a(8)*(1d0-eps)*epsd(nd)/bdat%d**2)/(-bdat%a(9)-eps&
&       )**2
  END DO
    sov = 4d0*(bdat%b/bdat%d)**2*pwr1*((pi/2d0-1d0/2d0*result1)*(1d0+2d0&
&     *eps**2)+3d0/2d0*eps*result2)*bdat%a(8)*(1d0-eps)/(-bdat%a(9)-eps)
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
        abs1d(nd) = epsd(nd)
    END DO
      abs1 = eps
  ELSE
    DO nd=1,nbdirs
        abs1d(nd) = -epsd(nd)
    END DO
      abs1 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
        abs5d(nd) = epsd(nd)
    END DO
      abs5 = eps
  ELSE
    DO nd=1,nbdirs
        abs5d(nd) = -epsd(nd)
    END DO
      abs5 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
        abs8d(nd) = epsd(nd)
    END DO
      abs8 = eps
  ELSE
    DO nd=1,nbdirs
        abs8d(nd) = -epsd(nd)
    END DO
      abs8 = -eps
  END IF
    arg1 = 1d0 - eps**2
    result1 = SQRT(arg1)
  arg2 = pi*result1
    arg3 = 2d0*abs1
  DO nd=1,nbdirs
! angle between force for rotation and minimum clearance according to Butenschoen 1976
      arg1d(nd) = -(2*eps*epsd(nd))
      IF (arg1 .EQ. 0.0) THEN
      result1d(nd) = 0.D0
    ELSE
        result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
    END IF
    arg2d(nd) = pi*result1d(nd)
      arg3d(nd) = 2d0*abs1d(nd)
      betad(nd) = (arg2d(nd)*arg3-arg3d(nd)*arg2)*(bdat%a(3)+bdat%a(4)*&
&       abs5+bdat%a(5)*eps**2+bdat%a(6)*abs8**3+bdat%a(7)*eps**4)/(arg2&
&       **2+arg3**2) + ATAN2(arg2, arg3)*(bdat%a(4)*abs5d(nd)+bdat%a(5)*&
&       2*eps*epsd(nd)+bdat%a(6)*3*abs8**2*abs8d(nd)+bdat%a(7)*4*eps**3*&
&       epsd(nd))
  END DO
    beta = ATAN2(arg2, arg3)*(bdat%a(3)+bdat%a(4)*abs5+bdat%a(5)*eps**2+&
&     bdat%a(6)*abs8**3+bdat%a(7)*eps**4)
  IF (eps .GE. 0.) THEN
      abs2 = eps
  ELSE
      abs2 = -eps
  END IF
    IF (abs2 .LT. eps_min) THEN
! avoid division infinite by infinite in case of zero relative eccentricity
! use analytical limit of abs_MR for eps going to zero
      mu = HUGE(1d0)
    DO nd=1,nbdirs
        mud(nd) = 0.D0
    END DO
  ELSE
      IF ((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot) .GE. 0.&
&     ) THEN
      DO nd=1,nbdirs
          abs3d(nd) = ((omegad(nd, 1)-omegad(nd, 2))*(omega(2)+omega(1)-&
&           2d0*delta_dot)-(omega(1)-omega(2))*(omegad(nd, 2)+omegad(nd&
&           , 1)-2d0*delta_dotd(nd)))/(omega(2)+omega(1)-2d0*delta_dot)&
&           **2
      END DO
        abs3 = (omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot)
    ELSE
      DO nd=1,nbdirs
          abs3d(nd) = -(((omegad(nd, 1)-omegad(nd, 2))*(omega(2)+omega(1&
&           )-2d0*delta_dot)-(omega(1)-omega(2))*(omegad(nd, 2)+omegad(&
&           nd, 1)-2d0*delta_dotd(nd)))/(omega(2)+omega(1)-2d0*delta_dot&
&           )**2)
      END DO
        abs3 = -((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot))
    END IF
    IF (eps .GE. 0.) THEN
      DO nd=1,nbdirs
          abs6d(nd) = epsd(nd)
      END DO
        abs6 = eps
    ELSE
      DO nd=1,nbdirs
          abs6d(nd) = -epsd(nd)
      END DO
        abs6 = -eps
    END IF
      arg1 = 1d0 - eps**2
      result1 = SQRT(arg1)
    DO nd=1,nbdirs
! friction coefficient according to Butenschoen
        arg1d(nd) = -(2*eps*epsd(nd))
        IF (arg1 .EQ. 0.0) THEN
        result1d(nd) = 0.D0
      ELSE
          result1d(nd) = arg1d(nd)/(2.0*SQRT(arg1))
      END IF
        mud(nd) = bdat%psi*((pi*abs3d(nd)*result1*sod-abs3*pi*(result1d(&
&         nd)*sod+result1*sodd(nd)))/(result1*sod)**2+(betad(nd)*COS(&
&         beta)*abs6+SIN(beta)*abs6d(nd))/2d0)
    END DO
      mu = bdat%psi*(abs3*pi/(result1*sod)+SIN(beta)*abs6/2d0)
  END IF
  END SUBROUTINE SOMMERFELD_NUMBERS_DV

  !  Differentiation of sommerfeld_numbers_d in forward (tangent) mode:
!   variations   of useful results: sodd sovd mu
!   with respect to varying inputs: omega delta_dot epsd
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
!  Differentiation of sommerfeld_numbers in forward (tangent) mode:
!   variations   of useful results: sod beta sov mu
!   with respect to varying inputs: eps
!   RW status of diff variables: eps:in sod:out beta:out sov:out
!                mu:out
  SUBROUTINE SOMMERFELD_NUMBERS_D_DV(bdat, eps, epsd, epsdd, omega, &
&   omegad, delta_dot, delta_dotd, sod, sodd, soddd, sov, sovd, sovdd, &
&   beta, mu, mud, nbdirs)
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
    TYPE(BEARING_DATA), INTENT(IN) :: bdat
    DOUBLE PRECISION, INTENT(IN) :: eps, delta_dot, omega(2)
    DOUBLE PRECISION, INTENT(IN) :: delta_dotd(nbdirsmax), omegad(&
&   nbdirsmax, 2)
    DOUBLE PRECISION, INTENT(IN) :: epsd
    DOUBLE PRECISION, DIMENSION(nbdirsmax), INTENT(IN) :: epsdd
    DOUBLE PRECISION, INTENT(OUT) :: sod, sov, beta, mu
    DOUBLE PRECISION, DIMENSION(nbdirsmax), INTENT(OUT) :: mud
    DOUBLE PRECISION, INTENT(OUT) :: sodd, sovd
    DOUBLE PRECISION, DIMENSION(nbdirsmax), INTENT(OUT) :: soddd, sovdd
    DOUBLE PRECISION, PARAMETER :: eps_min=1d-6
    INTRINSIC ABS
    INTRINSIC SQRT
    INTRINSIC ACOS
    INTRINSIC ATAN2
    INTRINSIC SIN
    INTRINSIC HUGE
    DOUBLE PRECISION :: arg1
    DOUBLE PRECISION :: arg1d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg1dd
  DOUBLE PRECISION :: result1
    DOUBLE PRECISION :: result1d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: result1dd
    DOUBLE PRECISION :: pwx1
    DOUBLE PRECISION :: pwx1d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: pwx1dd
    DOUBLE PRECISION :: pwr1
    DOUBLE PRECISION :: pwr1d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: pwr1dd
    DOUBLE PRECISION :: result2
    DOUBLE PRECISION :: result2d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: result2dd
    DOUBLE PRECISION :: abs1d
    DOUBLE PRECISION :: abs4d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: abs4dd
    DOUBLE PRECISION :: abs7d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: abs7dd
    DOUBLE PRECISION :: abs0d
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: abs0dd
    DOUBLE PRECISION :: abs6d
    DOUBLE PRECISION :: abs8
    DOUBLE PRECISION :: abs7
    DOUBLE PRECISION :: abs6
    DOUBLE PRECISION :: abs5
    DOUBLE PRECISION :: abs4
    DOUBLE PRECISION :: abs3
    DOUBLE PRECISION, DIMENSION(nbdirsmax) :: abs3d
    DOUBLE PRECISION :: abs2
    DOUBLE PRECISION :: abs1
    DOUBLE PRECISION :: abs0
    DOUBLE PRECISION :: abs5d
    DOUBLE PRECISION :: abs8d
    INTRINSIC INT
    DOUBLE PRECISION :: result10
    DOUBLE PRECISION :: pwr10
    DOUBLE PRECISION :: arg10
    DOUBLE PRECISION :: arg2
  INTEGER :: nd
  INTEGER :: nbdirs
    IF (eps .GE. 0.) THEN
  DO nd=1,nbdirs
        abs0dd(nd) = epsdd(nd)
  END DO
      abs0d = epsd
      abs0 = eps
    ELSE
      DO nd=1,nbdirs
        abs0dd(nd) = -epsdd(nd)
  END DO
      abs0d = -epsd
      abs0 = -eps
    END IF
    IF (eps .GE. 0.) THEN
  DO nd=1,nbdirs
        abs4dd(nd) = epsdd(nd)
  END DO
      abs4d = epsd
      abs4 = eps
  ELSE
      DO nd=1,nbdirs
        abs4dd(nd) = -epsdd(nd)
      END DO
      abs4d = -epsd
      abs4 = -eps
  END IF
    IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
        abs7dd(nd) = epsdd(nd)
    END DO
      abs7d = epsd
      abs7 = eps
  ELSE
    DO nd=1,nbdirs
        abs7dd(nd) = -epsdd(nd)
    END DO
      abs7d = -epsd
      abs7 = -eps
  END IF
  DO nd=1,nbdirs
! Sommerfeld number for rotation according to Butenschoen 1976
      arg1dd(nd) = 16d0*2*eps*epsdd(nd) - pi**2*2*eps*epsdd(nd)
  END DO
    arg1d = 16d0*2*eps*epsd - pi**2*2*eps*epsd
    arg1 = pi**2*(1d0-eps**2) + 16d0*eps**2
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.d0
    DO nd=1,nbdirs
        result1dd(nd) = 0.D0
    END DO
  ELSE
      result10 = SQRT(arg1)
    DO nd=1,nbdirs
        result1dd(nd) = arg1dd(nd)/(2.0*result10)
    END DO
      result1d = arg1d/(2.0*result10)
  END IF
    result1 = SQRT(arg1)
  DO nd=1,nbdirs
      soddd(nd) = (bdat%a(1)*(bdat%a(2)+abs7)*(result1*(abs4-1d0)*(bdat%&
&       b**2*2d0*(1d0-eps**2)**2*abs0dd(nd)/bdat%d**2+bdat%b**2*abs0*2d0&
&       *2**2*(1d0-eps**2)*eps*epsdd(nd)/bdat%d**2)/(2d0**2*(1d0-eps**2)&
&       **4)+bdat%b**2*abs0*((abs4-1d0)*result1dd(nd)+result1*abs4dd(nd)&
&       )/(bdat%d**2*2d0*(1d0-eps**2)**2))-bdat%b**2*abs0*result1*bdat%a&
&       (1)*(abs4-1d0)*abs7dd(nd)/(bdat%d**2*2d0*(1d0-eps**2)**2))/(bdat&
&       %a(2)+abs7)**2
! Sommerfeld number for displacement according to Butenschoen 1976
      pwx1dd(nd) = -(2*eps*epsdd(nd))
    END DO
    sodd = (bdat%a(1)*((bdat%b**2*abs0d*2d0*(1d0-eps**2)**2/bdat%d**2+&
&     bdat%b**2*abs0*2d0*2**2*(1d0-eps**2)*eps*epsd/bdat%d**2)*result1*(&
&     abs4-1d0)/(2d0**2*(1d0-eps**2)**4)+bdat%b**2*abs0*(result1d*(abs4-&
&     1d0)+result1*abs4d)/(bdat%d**2*2d0*(1d0-eps**2)**2))*(bdat%a(2)+&
&     abs7)-bdat%b**2*abs0*result1*bdat%a(1)*(abs4-1d0)*abs7d/(bdat%d**2&
&     *2d0*(1d0-eps**2)**2))/(bdat%a(2)+abs7)**2
    sod = (bdat%b/bdat%d)**2*abs0/(2d0*(1d0-eps**2)**2)*result1*bdat%a(1&
&     )*(abs4-1d0)/(bdat%a(2)+abs7)
    pwx1d = -(2*eps*epsd)
    pwx1 = 1d0 - eps**2
    IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. -(5d0/2d0) .EQ. INT(-(&
&       5d0/2d0)))) THEN
      pwr10 = pwx1**(-1-5d0/2d0)
    DO nd=1,nbdirs
        pwr1dd(nd) = -(5d0*pwr10*pwx1dd(nd)/2d0)
    END DO
      pwr1d = -(5d0*pwr10*pwx1d/2d0)
    ELSE IF (pwx1 .EQ. 0.0 .AND. -(5d0/2d0) .EQ. 1.0) THEN
    DO nd=1,nbdirs
        pwr1dd(nd) = pwx1dd(nd)
    END DO
      pwr1d = pwx1d
  ELSE
      pwr1d = 0.0
      DO nd=1,nbdirs
        pwr1dd(nd) = 0.D0
      END DO
  END IF
    pwr1 = pwx1**(-(5d0/2d0))
    IF (eps .EQ. 1.0 .OR. eps .EQ. -1.0) THEN
      result1d = 0.d0
      DO nd=1,nbdirs
        result1dd(nd) = 0.D0
      END DO
  ELSE
      arg10 = 1.0 - eps**2
      result10 = SQRT(arg10)
      DO nd=1,nbdirs
        result1dd(nd) = -(epsdd(nd)/result10)
      END DO
      result1d = -(epsd/result10)
  END IF
    result1 = ACOS(eps)
    DO nd=1,nbdirs
      arg1dd(nd) = -(2*eps*epsdd(nd))
    END DO
    arg1d = -(2*eps*epsd)
    arg1 = 1d0 - eps**2
    IF (arg1 .EQ. 0.0) THEN
      result2d = 0.d0
      DO nd=1,nbdirs
        result2dd(nd) = 0.D0
      END DO
  ELSE
      result10 = SQRT(arg1)
      DO nd=1,nbdirs
        result2dd(nd) = arg1dd(nd)/(2.0*result10)
      END DO
      result2d = arg1d/(2.0*result10)
  END IF
    result2 = SQRT(arg1)
    DO nd=1,nbdirs
      sovdd(nd) = (4d0*bdat%b**2*bdat%a(8)*(-bdat%a(9)-eps)*(((pi/2d0-&
&       1d0/2d0*result1)*(1d0+2d0*eps**2)+3d0/2d0*eps*result2)*((1d0-eps&
&       )*pwr1dd(nd)-pwr1*epsdd(nd))+pwr1*(1d0-eps)*((pi/2d0-1d0/2d0*&
&       result1)*2d0*2*eps*epsdd(nd)-(1d0+2d0*eps**2)*result1dd(nd)/2d0+&
&       3d0*(result2*epsdd(nd)+eps*result2dd(nd))/2d0))/bdat%d**2+4d0*&
&       bdat%b**2*pwr1*((pi/2d0-1d0/2d0*result1)*(1d0+2d0*eps**2)+3d0/&
&       2d0*eps*result2)*bdat%a(8)*(1d0-eps)*epsdd(nd)/bdat%d**2)/(-bdat&
&       %a(9)-eps)**2
    END DO
    sovd = (4d0*bdat%b**2*bdat%a(8)*((pwr1d*(1d0-eps)-pwr1*epsd)*((pi/&
&     2d0-1d0/2d0*result1)*(1d0+2d0*eps**2)+3d0/2d0*eps*result2)+pwr1*(&
&     1d0-eps)*((pi/2d0-1d0/2d0*result1)*2d0*2*eps*epsd-result1d*(1d0+&
&     2d0*eps**2)/2d0+3d0*(epsd*result2+eps*result2d)/2d0))*(-bdat%a(9)-&
&     eps)/bdat%d**2+4d0*bdat%b**2*pwr1*((pi/2d0-1d0/2d0*result1)*(1d0+&
&     2d0*eps**2)+3d0/2d0*eps*result2)*bdat%a(8)*(1d0-eps)*epsd/bdat%d**&
&     2)/(-bdat%a(9)-eps)**2
    sov = 4d0*(bdat%b/bdat%d)**2*pwr1*((pi/2d0-1d0/2d0*result1)*(1d0+2d0&
&     *eps**2)+3d0/2d0*eps*result2)*bdat%a(8)*(1d0-eps)/(-bdat%a(9)-eps)
  IF (eps .GE. 0.) THEN
      abs1d = epsd
      abs1 = eps
  ELSE
      abs1d = -epsd
      abs1 = -eps
  END IF
  IF (eps .GE. 0.) THEN
      abs5d = epsd
      abs5 = eps
  ELSE
      abs5d = -epsd
      abs5 = -eps
  END IF
  IF (eps .GE. 0.) THEN
      abs8d = epsd
      abs8 = eps
  ELSE
      abs8d = -epsd
      abs8 = -eps
  END IF
! angle between force for rotation and minimum clearance according to Butenschoen 1976
    arg1d = -(2*eps*epsd)
    arg1 = 1d0 - eps**2
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.d0
  ELSE
      result10 = SQRT(arg1)
      result1d = arg1d/(2.0*result10)
  END IF
    result1 = SQRT(arg1)
    arg10 = pi*result1
    arg2 = 2d0*abs1
    beta = ATAN2(arg10, arg2)*(bdat%a(3)+bdat%a(4)*abs5+bdat%a(5)*eps**2&
&     +bdat%a(6)*abs8**3+bdat%a(7)*eps**4)
  IF (eps .GE. 0.) THEN
      abs2 = eps
  ELSE
      abs2 = -eps
  END IF
    IF (abs2 .LT. eps_min) THEN
! avoid division infinite by infinite in case of zero relative eccentricity
! use analytical limit of abs_MR for eps going to zero
      mu = HUGE(1d0)
      DO nd=1,nbdirs
        mud(nd) = 0.D0
      END DO
  ELSE
      IF ((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot) .GE. 0.&
&     ) THEN
        DO nd=1,nbdirs
          abs3d(nd) = ((omegad(nd, 1)-omegad(nd, 2))*(omega(2)+omega(1)-&
&           2d0*delta_dot)-(omega(1)-omega(2))*(omegad(nd, 2)+omegad(nd&
&           , 1)-2d0*delta_dotd(nd)))/(omega(2)+omega(1)-2d0*delta_dot)&
&           **2
        END DO
        abs3 = (omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot)
    ELSE
        DO nd=1,nbdirs
          abs3d(nd) = -(((omegad(nd, 1)-omegad(nd, 2))*(omega(2)+omega(1&
&           )-2d0*delta_dot)-(omega(1)-omega(2))*(omegad(nd, 2)+omegad(&
&           nd, 1)-2d0*delta_dotd(nd)))/(omega(2)+omega(1)-2d0*delta_dot&
&           )**2)
        END DO
        abs3 = -((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot))
    END IF
    IF (eps .GE. 0.) THEN
        abs6d = epsd
        abs6 = eps
    ELSE
        abs6d = -epsd
        abs6 = -eps
    END IF
! friction coefficient according to Butenschoen
      arg1d = -(2*eps*epsd)
      arg1 = 1d0 - eps**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.d0
  ELSE
        result10 = SQRT(arg1)
        result1d = arg1d/(2.0*result10)
  END IF
      result1 = SQRT(arg1)
      DO nd=1,nbdirs
        mud(nd) = bdat%psi*pi*abs3d(nd)/(result1*sod)
      END DO
      mu = bdat%psi*(abs3*pi/(result1*sod)+SIN(beta)*abs6/2d0)
    END IF
  END SUBROUTINE SOMMERFELD_NUMBERS_D_DV
END MODULE HYDRODYNAMIC_PLAIN_BEARING_DV