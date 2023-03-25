/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2023) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef ___SOLID_SHAPE_H__INCLUDED___
#define ___SOLID_SHAPE_H__INCLUDED___

#include "sp_matrix_base.h"
#include "solid.h"

class Quadrangle4 {
public:
     static const char* ElementName() {
          return "quadrangle4";
     }

     static constexpr sp_grad::index_type iNumNodes = 4;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r);
private:
     static constexpr doublereal ri[] = {1, -1, -1,  1};
     static constexpr doublereal si[] = {1,  1, -1, -1};
};

class Quadrangle8 {
public:
     static const char* ElementName() {
          return "quadrangle8";
     }

     static constexpr sp_grad::index_type iNumNodes = 8;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r);
private:
     static constexpr doublereal ri[] = {1, -1, -1,  1, 0, -1,  0, 1};
     static constexpr doublereal si[] = {1,  1, -1, -1, 1,  0, -1, 0};
};

class Quadrangle8r {
public:
     static const char* ElementName() {
          return "quadrangle8r";
     }

     static constexpr sp_grad::index_type iNumNodes = 8;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r);
private:
     static constexpr doublereal ri[] = {-1,  1, 1, -1,  0, 1, 0, -1};
     static constexpr doublereal si[] = {-1, -1, 1,  1, -1, 0, 1,  0};
};

class Triangle6h {
public:
     static const char* ElementName() {
          return "triangle6";
     }

     static constexpr sp_grad::index_type iNumNodes = 6;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r);
private:
     static constexpr doublereal ri[] = {0, 1, 0, 1./2., 1./2.,     0};
     static constexpr doublereal si[] = {0, 0, 1,    0,  1./2., 1./2.};
};

class Hexahedron8 {
public:
     static const char* ElementName() {
          return "hexahedron8";
     }

     static constexpr sp_grad::index_type iNumNodes = 8;
     static constexpr sp_grad::index_type iNumNodesExtrap = iNumNodes;
     static constexpr bool bHaveDiagMass = true;

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumNodesExtrap, iNumComp>& taune);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r);

private:
     static constexpr doublereal ri[] = {1, -1, -1,  1,  1, -1, -1,  1};
     static constexpr doublereal si[] = {1,  1, -1, -1,  1,  1, -1, -1};
     static constexpr doublereal ti[] = {1,  1,  1,  1, -1, -1, -1, -1};
};

class Hexahedron20 {
public:
     static const char* ElementName() {
          return "hexahedron20";
     }

     static constexpr sp_grad::index_type iNumNodes = 20;
     static constexpr sp_grad::index_type iNumNodesExtrap = iNumNodes;
     static constexpr bool bHaveDiagMass = false; // not applicable to incomplete high order elements

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r);

private:
     static constexpr doublereal ri[] = {1, -1, -1,  1,  1, -1, -1,  1, 0, -1,  0, 1,  0, -1,  0,  1, 1, -1, -1,  1};
     static constexpr doublereal si[] = {1,  1, -1, -1,  1,  1, -1, -1, 1,  0, -1, 0,  1,  0, -1,  0, 1,  1, -1, -1};
     static constexpr doublereal ti[] = {1,  1,  1,  1, -1, -1, -1, -1, 1,  1,  1, 1, -1, -1, -1, -1, 0,  0,  0,  0};
};

class Hexahedron20r {
public:
     static const char* ElementName() {
          return "hexahedron20r";
     }

     static constexpr sp_grad::index_type iNumNodes = 20;
     static constexpr sp_grad::index_type iNumNodesExtrap = 8;
     static constexpr bool bHaveDiagMass = false; // not applicable to incomplete high order elements

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r);

private:
     static constexpr doublereal ri[] = {-1,  1,  1, -1, -1,  1, 1, -1,  0,  1,  0, -1,  0, 1, 0, -1, -1,  1, 1, -1};
     static constexpr doublereal si[] = {-1, -1,  1,  1, -1, -1, 1,  1, -1,  0,  1,  0, -1, 0, 1,  0, -1, -1, 1,  1};
     static constexpr doublereal ti[] = {-1, -1, -1, -1,  1,  1, 1,  1, -1, -1, -1, -1,  1, 1, 1,  1,  0,  0, 0,  0};
};

class Pentahedron15 {
public:
     static const char* ElementName() {
          return "pentahedron15";
     }

     static constexpr sp_grad::index_type iNumNodes = 15;
     static constexpr sp_grad::index_type iNumNodesExtrap = iNumNodes;
     static constexpr bool bHaveDiagMass = false; // to be checked

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r);

private:
     static constexpr doublereal ri[] = { 0,  1,  0, 0, 1, 0, 1./2., 1./2.,     0, 1./2., 1./2.,     0,  0, 1, 0};
     static constexpr doublereal si[] = { 0,  0,  1, 0, 0, 1,     0, 1./2., 1./2.,     0, 1./2., 1./2.,  0, 0, 1};
     static constexpr doublereal ti[] = {-1, -1, -1, 1, 1, 1,    -1,    -1,    -1,     1,     1,     1,  0, 0, 0};
};

class Tetrahedron10h {
public:
     static const char* ElementName() {
          return "tetrahedron10";
     }

     static constexpr sp_grad::index_type iNumNodes = 10;
     static constexpr sp_grad::index_type iNumNodesExtrap = 4;
     static constexpr bool bHaveDiagMass = false; // not supported by collocation rules

     static inline void
     ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                        sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1);

     static inline void
     ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                   sp_grad::SpColVector<doublereal, iNumNodes>& h);

     static inline void
     ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                         sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h);

     template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
     static inline void
     GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                        const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune);

     static inline void
     NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r);

private:
     static constexpr doublereal ri[] = {0, 0, 0,   1,   0,   0,   0, 0.5, 0.5, 0.5};
     static constexpr doublereal si[] = {1, 0, 0,   0, 0.5,   0, 0.5, 0.5,   0,   0};
     static constexpr doublereal ti[] = {0, 1, 0,   0, 0.5, 0.5,   0,   0, 0.5,   0};
};

void
Quadrangle4::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                                sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);

     static_assert(iNumNodes == 4);

     hd(1,1) = (r2+1)/4.0E+0;
     hd(1,2) = (r1+1)/4.0E+0;
     hd(2,1) = -(r2+1)/4.0E+0;
     hd(2,2) = (1-r1)/4.0E+0;
     hd(3,1) = -(1-r2)/4.0E+0;
     hd(3,2) = -(1-r1)/4.0E+0;
     hd(4,1) = (1-r2)/4.0E+0;
     hd(4,2) = -(r1+1)/4.0E+0;
}

void
Quadrangle4::ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                           sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);

     static_assert(iNumNodes == 4);

     h(1) = ((r1+1)*(r2+1))/4.0E+0;
     h(2) = ((1-r1)*(r2+1))/4.0E+0;
     h(3) = ((1-r1)*(1-r2))/4.0E+0;
     h(4) = ((r1+1)*(1-r2))/4.0E+0;
}

void
Quadrangle4::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
}

void
Quadrangle8::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                                sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r1_2 = r1 * r1;
     const doublereal r2_2 = r2 * r2;

     static_assert(iNumNodes == 8);

     hd(1,1) = (r2+1)/4.0E+0-((1-r2_2)/2.0E+0-r1*(r2+1))/2.0E+0;
     hd(1,2) = (r1+1)/4.0E+0-((1-r1_2)/2.0E+0-(r1+1)*r2)/2.0E+0;
     hd(2,1) = (-((-(1-r2_2)/2.0E+0)-r1*(r2+1))/2.0E+0)-(r2+1)/4.0E+0;
     hd(2,2) = (1-r1)/4.0E+0-((1-r1_2)/2.0E+0-(1-r1)*r2)/2.0E+0;
     hd(3,1) = (-((-(1-r2_2)/2.0E+0)-r1*(1-r2))/2.0E+0)-(1-r2)/4.0E+0;
     hd(3,2) = (-((-(1-r1)*r2)-(1-r1_2)/2.0E+0)/2.0E+0)-(1-r1)/4.0E+0;
     hd(4,1) = (1-r2)/4.0E+0-((1-r2_2)/2.0E+0-r1*(1-r2))/2.0E+0;
     hd(4,2) = (-((-(r1+1)*r2)-(1-r1_2)/2.0E+0)/2.0E+0)-(r1+1)/4.0E+0;
     hd(5,1) = -r1*(r2+1);
     hd(5,2) = (1-r1_2)/2.0E+0;
     hd(6,1) = -(1-r2_2)/2.0E+0;
     hd(6,2) = -(1-r1)*r2;
     hd(7,1) = -r1*(1-r2);
     hd(7,2) = -(1-r1_2)/2.0E+0;
     hd(8,1) = (1-r2_2)/2.0E+0;
     hd(8,2) = -(r1+1)*r2;
}

void
Quadrangle8::ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                           sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r1_2 = r1 * r1;
     const doublereal r2_2 = r2 * r2;

     static_assert(iNumNodes == 8);

     h(1) = ((r1+1)*(r2+1))/4.0E+0-(((r1+1)*(1-r2_2))/2.0E+0+((1-r1_2)*(r2+1))/2.0E+0)/2.0E+0;
     h(2) = ((1-r1)*(r2+1))/4.0E+0-(((1-r1)*(1-r2_2))/2.0E+0+((1-r1_2)*(r2+1))/2.0E+0)/2.0E+0;
     h(3) = ((1-r1)*(1-r2))/4.0E+0-(((1-r1)*(1-r2_2))/2.0E+0+((1-r1_2)*(1-r2))/2.0E+0)/2.0E+0;
     h(4) = ((r1+1)*(1-r2))/4.0E+0-(((r1+1)*(1-r2_2))/2.0E+0+((1-r1_2)*(1-r2))/2.0E+0)/2.0E+0;
     h(5) = ((1-r1_2)*(r2+1))/2.0E+0;
     h(6) = ((1-r1)*(1-r2_2))/2.0E+0;
     h(7) = ((1-r1_2)*(1-r2))/2.0E+0;
     h(8) = ((r1+1)*(1-r2_2))/2.0E+0;
}

void
Quadrangle8::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
}

void
Quadrangle8r::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                                 sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);

     static_assert(iNumNodes == 8);

     hd(1,1) = 2.5E-1*(1.0E+0-r2)*(r2+r1+1.0E+0)+2.5E-1*(r1-1.0E+0)*(1.0E+0-r2);
     hd(1,2) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)-2.5E-1*(r1-1.0E+0)*(r2+r1+1.0E+0);
     hd(2,1) = (-2.5E-1*(1.0E+0-r2)*(r2-r1+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2);
     hd(2,2) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)-2.5E-1*((-r1)-1.0E+0)*(r2-r1+1.0E+0);
     hd(3,1) = (-2.5E-1*((-r2)-r1+1.0E+0)*(r2+1.0E+0))-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
     hd(3,2) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)-2.5E-1*((-r1)-1.0E+0)*(r2+1.0E+0);
     hd(4,1) = 2.5E-1*((-r2)+r1+1.0E+0)*(r2+1.0E+0)+2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
     hd(4,2) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)-2.5E-1*(r1-1.0E+0)*(r2+1.0E+0);
     hd(5,1) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(1.0E+0-r2);
     hd(5,2) = -5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
     hd(6,1) = 5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
     hd(6,2) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
     hd(7,1) = 5.0E-1*(1.0E+0-r1)*(r2+1.0E+0)-5.0E-1*(r1+1.0E+0)*(r2+1.0E+0);
     hd(7,2) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0);
     hd(8,1) = -5.0E-1*(1.0E+0-r2)*(r2+1.0E+0);
     hd(8,2) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)-5.0E-1*(1.0E+0-r1)*(r2+1.0E+0);
}

void
Quadrangle8r::ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                            sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);

     static_assert(iNumNodes == 8);

     h(1) = 2.5E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r2+r1+1.0E+0);
     h(2) = 2.5E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r2-r1+1.0E+0);
     h(3) = 2.5E-1*((-r1)-1.0E+0)*((-r2)-r1+1.0E+0)*(r2+1.0E+0);
     h(4) = 2.5E-1*(r1-1.0E+0)*((-r2)+r1+1.0E+0)*(r2+1.0E+0);
     h(5) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
     h(6) = 5.0E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
     h(7) = 5.0E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
     h(8) = 5.0E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
}

void
Quadrangle8r::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
}

void
Triangle6h::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 2>& r,
                               sp_grad::SpMatrix<doublereal, iNumNodes, 2>& hd)
{
     const doublereal zeta = r(1);
     const doublereal eta = r(2);

     static_assert(iNumNodes == 6);

     hd(1,1) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
     hd(1,2) = 2*(zeta+eta-1)-2*((-zeta)-eta+1)+1;
     hd(2,1) = 4*zeta-1;
     hd(2,2) = 0;
     hd(3,1) = 0;
     hd(3,2) = 4*eta-1;
     hd(4,1) = 4*((-zeta)-eta+1)-4*zeta;
     hd(4,2) = -4*zeta;
     hd(5,1) = 4*eta;
     hd(5,2) = 4*zeta;
     hd(6,1) = -4*eta;
     hd(6,2) = 4*((-zeta)-eta+1)-4*eta;
}

void
Triangle6h::ShapeFunction(const sp_grad::SpColVector<doublereal, 2>& r,
                          sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal zeta = r(1);
     const doublereal eta = r(2);

     static_assert(iNumNodes == 6);

     h(1) = (1-2*((-zeta)-eta+1))*(zeta+eta-1);
     h(2) = -(1-2*zeta)*zeta;
     h(3) = -(1-2*eta)*eta;
     h(4) = 4*((-zeta)-eta+1)*zeta;
     h(5) = 4*eta*zeta;
     h(6) = 4*eta*((-zeta)-eta+1);
}

void
Triangle6h::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 2>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
}

void
Hexahedron8::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     static_assert(iNumNodes == 8);

     h0d1(1,1) = ((r(2)+1)*(r(3)+1))/8.0E+0;
     h0d1(1,2) = ((r(1)+1)*(r(3)+1))/8.0E+0;
     h0d1(1,3) = ((r(1)+1)*(r(2)+1))/8.0E+0;
     h0d1(2,1) = -((r(2)+1)*(r(3)+1))/8.0E+0;
     h0d1(2,2) = ((1-r(1))*(r(3)+1))/8.0E+0;
     h0d1(2,3) = ((1-r(1))*(r(2)+1))/8.0E+0;
     h0d1(3,1) = -((1-r(2))*(r(3)+1))/8.0E+0;
     h0d1(3,2) = -((1-r(1))*(r(3)+1))/8.0E+0;
     h0d1(3,3) = ((1-r(1))*(1-r(2)))/8.0E+0;
     h0d1(4,1) = ((1-r(2))*(r(3)+1))/8.0E+0;
     h0d1(4,2) = -((r(1)+1)*(r(3)+1))/8.0E+0;
     h0d1(4,3) = ((r(1)+1)*(1-r(2)))/8.0E+0;
     h0d1(5,1) = ((r(2)+1)*(1-r(3)))/8.0E+0;
     h0d1(5,2) = ((r(1)+1)*(1-r(3)))/8.0E+0;
     h0d1(5,3) = -((r(1)+1)*(r(2)+1))/8.0E+0;
     h0d1(6,1) = -((r(2)+1)*(1-r(3)))/8.0E+0;
     h0d1(6,2) = ((1-r(1))*(1-r(3)))/8.0E+0;
     h0d1(6,3) = -((1-r(1))*(r(2)+1))/8.0E+0;
     h0d1(7,1) = -((1-r(2))*(1-r(3)))/8.0E+0;
     h0d1(7,2) = -((1-r(1))*(1-r(3)))/8.0E+0;
     h0d1(7,3) = -((1-r(1))*(1-r(2)))/8.0E+0;
     h0d1(8,1) = ((1-r(2))*(1-r(3)))/8.0E+0;
     h0d1(8,2) = -((r(1)+1)*(1-r(3)))/8.0E+0;
     h0d1(8,3) = -((r(1)+1)*(1-r(2)))/8.0E+0;
}

void
Hexahedron8::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                           sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     static_assert(iNumNodes == 8);

     h(1) = ((r(1)+1)*(r(2)+1)*(r(3)+1))/8.0E+0;
     h(2) = ((1-r(1))*(r(2)+1)*(r(3)+1))/8.0E+0;
     h(3) = ((1-r(1))*(1-r(2))*(r(3)+1))/8.0E+0;
     h(4) = ((r(1)+1)*(1-r(2))*(r(3)+1))/8.0E+0;
     h(5) = ((r(1)+1)*(r(2)+1)*(1-r(3)))/8.0E+0;
     h(6) = ((1-r(1))*(r(2)+1)*(1-r(3)))/8.0E+0;
     h(7) = ((1-r(1))*(1-r(2))*(1-r(3)))/8.0E+0;
     h(8) = ((r(1)+1)*(1-r(2))*(1-r(3)))/8.0E+0;
}

void
Hexahedron8::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                 sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     ShapeFunction(r, h);
}

template <sp_grad::index_type iNumComp>
void
Hexahedron8::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                const sp_grad::SpMatrix<doublereal, iNumNodesExtrap, iNumComp>& taune)
{
     tauni = taune;
}

void
Hexahedron8::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
     r(3) = ti[iNode - 1];
}

void
Hexahedron20::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                 sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r1_2 = r1 * r1;
     const doublereal r2_2 = r2 * r2;
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 20);

     h0d1(1,1) = ((r2+1)*(r3+1))/8.0E+0-(((r2+1)*(1-r3_2))/4.0E+0+((1-r2_2)*(r3+1))/4.0E+0-(r1*(r2+1)*(r3+1))/2.0E+0)/2.0E+0;
     h0d1(1,2) = ((r1+1)*(r3+1))/8.0E+0-(((r1+1)*(1-r3_2))/4.0E+0-((r1+1)*r2*(r3+1))/2.0E+0+((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0;
     h0d1(1,3) = ((r1+1)*(r2+1))/8.0E+0-((-((r1+1)*(r2+1)*r3)/2.0E+0)+((r1+1)*(1-r2_2))/4.0E+0+((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0;
     h0d1(2,1) = (-((-((r2+1)*(1-r3_2))/4.0E+0)-((1-r2_2)*(r3+1))/4.0E+0-(r1*(r2+1)*(r3+1))/2.0E+0)/2.0E+0)-((r2+1)*(r3+1))/8.0E+0;
     h0d1(2,2) = ((1-r1)*(r3+1))/8.0E+0-(((1-r1)*(1-r3_2))/4.0E+0-((1-r1)*r2*(r3+1))/2.0E+0+((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0;
     h0d1(2,3) = ((1-r1)*(r2+1))/8.0E+0-((-((1-r1)*(r2+1)*r3)/2.0E+0)+((1-r1)*(1-r2_2))/4.0E+0+((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0;
     h0d1(3,1) = (-((-((1-r2)*(1-r3_2))/4.0E+0)-((1-r2_2)*(r3+1))/4.0E+0-(r1*(1-r2)*(r3+1))/2.0E+0)/2.0E+0)-((1-r2)*(r3+1))/8.0E+0;
     h0d1(3,2) = (-((-((1-r1)*(1-r3_2))/4.0E+0)-((1-r1)*r2*(r3+1))/2.0E+0-((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0)-((1-r1)*(r3+1))/8.0E+0;
     h0d1(3,3) = ((1-r1)*(1-r2))/8.0E+0-((-((1-r1)*(1-r2)*r3)/2.0E+0)+((1-r1)*(1-r2_2))/4.0E+0+((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0;
     h0d1(4,1) = ((1-r2)*(r3+1))/8.0E+0-(((1-r2)*(1-r3_2))/4.0E+0+((1-r2_2)*(r3+1))/4.0E+0-(r1*(1-r2)*(r3+1))/2.0E+0)/2.0E+0;
     h0d1(4,2) = (-((-((r1+1)*(1-r3_2))/4.0E+0)-((r1+1)*r2*(r3+1))/2.0E+0-((1-r1_2)*(r3+1))/4.0E+0)/2.0E+0)-((r1+1)*(r3+1))/8.0E+0;
     h0d1(4,3) = ((r1+1)*(1-r2))/8.0E+0-((-((r1+1)*(1-r2)*r3)/2.0E+0)+((r1+1)*(1-r2_2))/4.0E+0+((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0;
     h0d1(5,1) = ((r2+1)*(1-r3))/8.0E+0-(((r2+1)*(1-r3_2))/4.0E+0+((1-r2_2)*(1-r3))/4.0E+0-(r1*(r2+1)*(1-r3))/2.0E+0)/2.0E+0;
     h0d1(5,2) = ((r1+1)*(1-r3))/8.0E+0-(((r1+1)*(1-r3_2))/4.0E+0-((r1+1)*r2*(1-r3))/2.0E+0+((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0;
     h0d1(5,3) = (-((-((r1+1)*(r2+1)*r3)/2.0E+0)-((r1+1)*(1-r2_2))/4.0E+0-((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0)-((r1+1)*(r2+1))/8.0E+0;
     h0d1(6,1) = (-((-((r2+1)*(1-r3_2))/4.0E+0)-((1-r2_2)*(1-r3))/4.0E+0-(r1*(r2+1)*(1-r3))/2.0E+0)/2.0E+0)-((r2+1)*(1-r3))/8.0E+0;
     h0d1(6,2) = ((1-r1)*(1-r3))/8.0E+0-(((1-r1)*(1-r3_2))/4.0E+0-((1-r1)*r2*(1-r3))/2.0E+0+((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0;
     h0d1(6,3) = (-((-((1-r1)*(r2+1)*r3)/2.0E+0)-((1-r1)*(1-r2_2))/4.0E+0-((1-r1_2)*(r2+1))/4.0E+0)/2.0E+0)-((1-r1)*(r2+1))/8.0E+0;
     h0d1(7,1) = (-((-((1-r2)*(1-r3_2))/4.0E+0)-((1-r2_2)*(1-r3))/4.0E+0-(r1*(1-r2)*(1-r3))/2.0E+0)/2.0E+0)-((1-r2)*(1-r3))/8.0E+0;
     h0d1(7,2) = (-((-((1-r1)*(1-r3_2))/4.0E+0)-((1-r1)*r2*(1-r3))/2.0E+0-((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0)-((1-r1)*(1-r3))/8.0E+0;
     h0d1(7,3) = (-((-((1-r1)*(1-r2)*r3)/2.0E+0)-((1-r1)*(1-r2_2))/4.0E+0-((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0)-((1-r1)*(1-r2))/8.0E+0;
     h0d1(8,1) = ((1-r2)*(1-r3))/8.0E+0-(((1-r2)*(1-r3_2))/4.0E+0+((1-r2_2)*(1-r3))/4.0E+0-(r1*(1-r2)*(1-r3))/2.0E+0)/2.0E+0;
     h0d1(8,2) = (-((-((r1+1)*(1-r3_2))/4.0E+0)-((r1+1)*r2*(1-r3))/2.0E+0-((1-r1_2)*(1-r3))/4.0E+0)/2.0E+0)-((r1+1)*(1-r3))/8.0E+0;
     h0d1(8,3) = (-((-((r1+1)*(1-r2)*r3)/2.0E+0)-((r1+1)*(1-r2_2))/4.0E+0-((1-r1_2)*(1-r2))/4.0E+0)/2.0E+0)-((r1+1)*(1-r2))/8.0E+0;
     h0d1(9,1) = -(r1*(r2+1)*(r3+1))/2.0E+0;
     h0d1(9,2) = ((1-r1_2)*(r3+1))/4.0E+0;
     h0d1(9,3) = ((1-r1_2)*(r2+1))/4.0E+0;
     h0d1(10,1) = -((1-r2_2)*(r3+1))/4.0E+0;
     h0d1(10,2) = -((1-r1)*r2*(r3+1))/2.0E+0;
     h0d1(10,3) = ((1-r1)*(1-r2_2))/4.0E+0;
     h0d1(11,1) = -(r1*(1-r2)*(r3+1))/2.0E+0;
     h0d1(11,2) = -((1-r1_2)*(r3+1))/4.0E+0;
     h0d1(11,3) = ((1-r1_2)*(1-r2))/4.0E+0;
     h0d1(12,1) = ((1-r2_2)*(r3+1))/4.0E+0;
     h0d1(12,2) = -((r1+1)*r2*(r3+1))/2.0E+0;
     h0d1(12,3) = ((r1+1)*(1-r2_2))/4.0E+0;
     h0d1(13,1) = -(r1*(r2+1)*(1-r3))/2.0E+0;
     h0d1(13,2) = ((1-r1_2)*(1-r3))/4.0E+0;
     h0d1(13,3) = -((1-r1_2)*(r2+1))/4.0E+0;
     h0d1(14,1) = -((1-r2_2)*(1-r3))/4.0E+0;
     h0d1(14,2) = -((1-r1)*r2*(1-r3))/2.0E+0;
     h0d1(14,3) = -((1-r1)*(1-r2_2))/4.0E+0;
     h0d1(15,1) = -(r1*(1-r2)*(1-r3))/2.0E+0;
     h0d1(15,2) = -((1-r1_2)*(1-r3))/4.0E+0;
     h0d1(15,3) = -((1-r1_2)*(1-r2))/4.0E+0;
     h0d1(16,1) = ((1-r2_2)*(1-r3))/4.0E+0;
     h0d1(16,2) = -((r1+1)*r2*(1-r3))/2.0E+0;
     h0d1(16,3) = -((r1+1)*(1-r2_2))/4.0E+0;
     h0d1(17,1) = ((r2+1)*(1-r3_2))/4.0E+0;
     h0d1(17,2) = ((r1+1)*(1-r3_2))/4.0E+0;
     h0d1(17,3) = -((r1+1)*(r2+1)*r3)/2.0E+0;
     h0d1(18,1) = -((r2+1)*(1-r3_2))/4.0E+0;
     h0d1(18,2) = ((1-r1)*(1-r3_2))/4.0E+0;
     h0d1(18,3) = -((1-r1)*(r2+1)*r3)/2.0E+0;
     h0d1(19,1) = -((1-r2)*(1-r3_2))/4.0E+0;
     h0d1(19,2) = -((1-r1)*(1-r3_2))/4.0E+0;
     h0d1(19,3) = -((1-r1)*(1-r2)*r3)/2.0E+0;
     h0d1(20,1) = ((1-r2)*(1-r3_2))/4.0E+0;
     h0d1(20,2) = -((r1+1)*(1-r3_2))/4.0E+0;
     h0d1(20,3) = -((r1+1)*(1-r2)*r3)/2.0E+0;
}

void
Hexahedron20::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                            sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r1_2 = r1 * r1;
     const doublereal r2_2 = r2 * r2;
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 20);

     h(1) = ((r1+1)*(r2+1)*(r3+1))/8.0E+0-(((r1+1)*(r2+1)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(r2+1)*(r3+1))/4.0E+0)/2.0E+0;
     h(2) = ((1-r1)*(r2+1)*(r3+1))/8.0E+0-(((1-r1)*(r2+1)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(r2+1)*(r3+1))/4.0E+0)/2.0E+0;
     h(3) = ((1-r1)*(1-r2)*(r3+1))/8.0E+0-(((1-r1)*(1-r2)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(1-r2)*(r3+1))/4.0E+0)/2.0E+0;
     h(4) = ((r1+1)*(1-r2)*(r3+1))/8.0E+0-(((r1+1)*(1-r2)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(r3+1))/4.0E+0+((1-r1_2)*(1-r2)*(r3+1))/4.0E+0)/2.0E+0;
     h(5) = ((r1+1)*(r2+1)*(1-r3))/8.0E+0-(((r1+1)*(r2+1)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(r2+1)*(1-r3))/4.0E+0)/2.0E+0;
     h(6) = ((1-r1)*(r2+1)*(1-r3))/8.0E+0-(((1-r1)*(r2+1)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(r2+1)*(1-r3))/4.0E+0)/2.0E+0;
     h(7) = ((1-r1)*(1-r2)*(1-r3))/8.0E+0-(((1-r1)*(1-r2)*(1-r3_2))/4.0E+0+((1-r1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(1-r2)*(1-r3))/4.0E+0)/2.0E+0;
     h(8) = ((r1+1)*(1-r2)*(1-r3))/8.0E+0-(((r1+1)*(1-r2)*(1-r3_2))/4.0E+0+((r1+1)*(1-r2_2)*(1-r3))/4.0E+0+((1-r1_2)*(1-r2)*(1-r3))/4.0E+0)/2.0E+0;
     h(9) = ((1-r1_2)*(r2+1)*(r3+1))/4.0E+0;
     h(10) = ((1-r1)*(1-r2_2)*(r3+1))/4.0E+0;
     h(11) = ((1-r1_2)*(1-r2)*(r3+1))/4.0E+0;
     h(12) = ((r1+1)*(1-r2_2)*(r3+1))/4.0E+0;
     h(13) = ((1-r1_2)*(r2+1)*(1-r3))/4.0E+0;
     h(14) = ((1-r1)*(1-r2_2)*(1-r3))/4.0E+0;
     h(15) = ((1-r1_2)*(1-r2)*(1-r3))/4.0E+0;
     h(16) = ((r1+1)*(1-r2_2)*(1-r3))/4.0E+0;
     h(17) = ((r1+1)*(r2+1)*(1-r3_2))/4.0E+0;
     h(18) = ((1-r1)*(r2+1)*(1-r3_2))/4.0E+0;
     h(19) = ((1-r1)*(1-r2)*(1-r3_2))/4.0E+0;
     h(20) = ((r1+1)*(1-r2)*(1-r3_2))/4.0E+0;
}

void
Hexahedron20::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                  sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     ShapeFunction(r, h);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Hexahedron20::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                 const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune)
{
     static_assert(iNumRhs >= iNumNodes);
     using namespace sp_grad;

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodes; ++i) {
               tauni(i, j) = taune(i, j);
          }
     }
}

void
Hexahedron20::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
     r(3) = ti[iNode - 1];
}

void
Hexahedron20r::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                  sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     static_assert(iNumNodes == 20);

     h0d1(1,1) = 1.25E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0)+1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h0d1(1,2) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0);
     h0d1(1,3) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+r2+r1+2.0E+0);
     h0d1(2,1) = (-1.25E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0))-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h0d1(2,2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0);
     h0d1(2,3) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+r2-r1+2.0E+0);
     h0d1(3,1) = (-1.25E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0))-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(3,2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(3,3) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3-r2-r1+2.0E+0);
     h0d1(4,1) = 1.25E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0)+1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(4,2) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(4,3) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3-r2+r1+2.0E+0);
     h0d1(5,1) = 1.25E-1*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0)+1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(5,2) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-1.25E-1*(r1-1.0E+0)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0);
     h0d1(5,3) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(6,1) = (-1.25E-1*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0))-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(6,2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-1.25E-1*((-r1)-1.0E+0)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0);
     h0d1(6,3) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(7,1) = (-1.25E-1*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0))-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(7,2) = 1.25E-1*((-r1)-1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(7,3) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)-1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(8,1) = 1.25E-1*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0)+1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(8,2) = 1.25E-1*(r1-1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(8,3) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)-1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(9,1) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h0d1(9,2) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r3);
     h0d1(9,3) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
     h0d1(10,1) = 2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(10,2) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(10,3) = -2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(11,1) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(11,2) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r3);
     h0d1(11,3) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
     h0d1(12,1) = -2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(12,2) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3);
     h0d1(12,3) = -2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(13,1) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(13,2) = -2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r3+1.0E+0);
     h0d1(13,3) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2);
     h0d1(14,1) = 2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(14,2) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(14,3) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(15,1) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(15,2) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r3+1.0E+0);
     h0d1(15,3) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0);
     h0d1(16,1) = -2.5E-1*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(16,2) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(16,3) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0);
     h0d1(17,1) = -2.5E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(17,2) = -2.5E-1*(1.0E+0-r1)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(17,3) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(18,1) = 2.5E-1*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(18,2) = -2.5E-1*(r1+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(18,3) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h0d1(19,1) = 2.5E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(19,2) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(19,3) = 2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h0d1(20,1) = -2.5E-1*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(20,2) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r3)*(r3+1.0E+0);
     h0d1(20,3) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)-2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
}

void
Hexahedron20r::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                             sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     static_assert(iNumNodes == 20);

     h(1) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2+r1+2.0E+0);
     h(2) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+r2-r1+2.0E+0);
     h(3) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2-r1+2.0E+0);
     h(4) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3-r2+r1+2.0E+0);
     h(5) = 1.25E-1*(r1-1.0E+0)*(1.0E+0-r2)*((-r3)+r2+r1+2.0E+0)*(r3+1.0E+0);
     h(6) = 1.25E-1*((-r1)-1.0E+0)*(1.0E+0-r2)*((-r3)+r2-r1+2.0E+0)*(r3+1.0E+0);
     h(7) = 1.25E-1*((-r1)-1.0E+0)*(r2+1.0E+0)*((-r3)-r2-r1+2.0E+0)*(r3+1.0E+0);
     h(8) = 1.25E-1*(r1-1.0E+0)*(r2+1.0E+0)*((-r3)-r2+r1+2.0E+0)*(r3+1.0E+0);
     h(9) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h(10) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h(11) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h(12) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0)*(1.0E+0-r3);
     h(13) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h(14) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h(15) = 2.5E-1*(1.0E+0-r1)*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h(16) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r2+1.0E+0)*(r3+1.0E+0);
     h(17) = 2.5E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h(18) = 2.5E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3)*(r3+1.0E+0);
     h(19) = 2.5E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
     h(20) = 2.5E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3)*(r3+1.0E+0);
}

void
Hexahedron20r::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                   sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     static_assert(iNumNodesExtrap == 8);

     h(1) = 1.25E-1*(1.0E+0-r1)*(1.0E+0-r2)*(1.0E+0-r3);
     h(2) = 1.25E-1*(r1+1.0E+0)*(1.0E+0-r2)*(1.0E+0-r3);
     h(3) = 1.25E-1*(r1+1.0E+0)*(r2+1.0E+0)*(1.0E+0-r3);
     h(4) = 1.25E-1*(1.0E+0-r1)*(r2+1.0E+0)*(1.0E+0-r3);
     h(5) = 1.25E-1*(1.0E+0-r1)*(1.0E+0-r2)*(r3+1.0E+0);
     h(6) = 1.25E-1*(r1+1.0E+0)*(1.0E+0-r2)*(r3+1.0E+0);
     h(7) = 1.25E-1*(r1+1.0E+0)*(r2+1.0E+0)*(r3+1.0E+0);
     h(8) = 1.25E-1*(1.0E+0-r1)*(r2+1.0E+0)*(r3+1.0E+0);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Hexahedron20r::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                  const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune)
{
     using namespace sp_grad;

     static_assert(iNumRhs >= iNumNodesExtrap);

     static constexpr struct {
          index_type ico1, ico2, imid;
     } idxint[] = {
          {1, 2,  9},
          {2, 3, 10},
          {3, 4, 11},
          {4, 1, 12},
          {5, 6, 13},
          {6, 7, 14},
          {7, 8, 15},
          {8, 5, 16},
          {1, 5, 17},
          {2, 6, 18},
          {3, 7, 19},
          {4, 8, 20}
     };

     constexpr index_type iNumNodesInterp = sizeof(idxint) / sizeof(idxint[0]);

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodesExtrap; ++i) {
               tauni(i, j) = taune(i, j);
          }
          for (index_type i = 0; i < iNumNodesInterp; ++i) {
               tauni(idxint[i].imid, j) = 0.5 * (taune(idxint[i].ico1, j)
                                                 +  taune(idxint[i].ico2, j));
          }
     }
}

void
Hexahedron20r::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
     r(3) = ti[iNode - 1];
}

void
Pentahedron15::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                  sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 15);

     h0d1(1,1) = (-((1-r3)*((-r3)-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(1-r3);
     h0d1(1,2) = (-((1-r3)*((-r3)-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(1-r3);
     h0d1(1,3) = (-(((-r2)-r1+1)*((-r3)-2*r2-2*r1))/2.0E+0)-(((-r2)-r1+1)*(1-r3))/2.0E+0;
     h0d1(2,1) = ((1-r3)*((-r3)+2*r1-2))/2.0E+0+r1*(1-r3);
     h0d1(2,2) = 0;
     h0d1(2,3) = (-(r1*((-r3)+2*r1-2))/2.0E+0)-(r1*(1-r3))/2.0E+0;
     h0d1(3,1) = 0;
     h0d1(3,2) = ((1-r3)*((-r3)+2*r2-2))/2.0E+0+r2*(1-r3);
     h0d1(3,3) = (-(r2*((-r3)+2*r2-2))/2.0E+0)-(r2*(1-r3))/2.0E+0;
     h0d1(4,1) = (-((r3+1)*(r3-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(r3+1);
     h0d1(4,2) = (-((r3+1)*(r3-2*r2-2*r1))/2.0E+0)-((-r2)-r1+1)*(r3+1);
     h0d1(4,3) = (((-r2)-r1+1)*(r3-2*r2-2*r1))/2.0E+0+(((-r2)-r1+1)*(r3+1))/2.0E+0;
     h0d1(5,1) = ((r3+1)*(r3+2*r1-2))/2.0E+0+r1*(r3+1);
     h0d1(5,2) = 0;
     h0d1(5,3) = (r1*(r3+2*r1-2))/2.0E+0+(r1*(r3+1))/2.0E+0;
     h0d1(6,1) = 0;
     h0d1(6,2) = ((r3+1)*(r3+2*r2-2))/2.0E+0+r2*(r3+1);
     h0d1(6,3) = (r2*(r3+2*r2-2))/2.0E+0+(r2*(r3+1))/2.0E+0;
     h0d1(7,1) = 2*((-r2)-r1+1)*(1-r3)-2*r1*(1-r3);
     h0d1(7,2) = -2*r1*(1-r3);
     h0d1(7,3) = -2*r1*((-r2)-r1+1);
     h0d1(8,1) = 2*r2*(1-r3);
     h0d1(8,2) = 2*r1*(1-r3);
     h0d1(8,3) = -2*r1*r2;
     h0d1(9,1) = -2*r2*(1-r3);
     h0d1(9,2) = 2*((-r2)-r1+1)*(1-r3)-2*r2*(1-r3);
     h0d1(9,3) = -2*((-r2)-r1+1)*r2;
     h0d1(10,1) = 2*((-r2)-r1+1)*(r3+1)-2*r1*(r3+1);
     h0d1(10,2) = -2*r1*(r3+1);
     h0d1(10,3) = 2*r1*((-r2)-r1+1);
     h0d1(11,1) = 2*r2*(r3+1);
     h0d1(11,2) = 2*r1*(r3+1);
     h0d1(11,3) = 2*r1*r2;
     h0d1(12,1) = -2*r2*(r3+1);
     h0d1(12,2) = 2*((-r2)-r1+1)*(r3+1)-2*r2*(r3+1);
     h0d1(12,3) = 2*((-r2)-r1+1)*r2;
     h0d1(13,1) = r3_2-1;
     h0d1(13,2) = r3_2-1;
     h0d1(13,3) = -2*((-r2)-r1+1)*r3;
     h0d1(14,1) = 1-r3_2;
     h0d1(14,2) = 0;
     h0d1(14,3) = -2*r1*r3;
     h0d1(15,1) = 0;
     h0d1(15,2) = 1-r3_2;
     h0d1(15,3) = -2*r2*r3;
}

void
Pentahedron15::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                   sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     ShapeFunction(r, h);
}

void
Pentahedron15::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                             sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);
     const doublereal r3_2 = r3 * r3;

     static_assert(iNumNodes == 15);

     h(1) = (((-r2)-r1+1)*(1-r3)*((-r3)-2*r2-2*r1))/2.0E+0;
     h(2) = (r1*(1-r3)*((-r3)+2*r1-2))/2.0E+0;
     h(3) = (r2*(1-r3)*((-r3)+2*r2-2))/2.0E+0;
     h(4) = (((-r2)-r1+1)*(r3+1)*(r3-2*r2-2*r1))/2.0E+0;
     h(5) = (r1*(r3+1)*(r3+2*r1-2))/2.0E+0;
     h(6) = (r2*(r3+1)*(r3+2*r2-2))/2.0E+0;
     h(7) = 2*r1*((-r2)-r1+1)*(1-r3);
     h(8) = 2*r1*r2*(1-r3);
     h(9) = 2*((-r2)-r1+1)*r2*(1-r3);
     h(10) = 2*r1*((-r2)-r1+1)*(r3+1);
     h(11) = 2*r1*r2*(r3+1);
     h(12) = 2*((-r2)-r1+1)*r2*(r3+1);
     h(13) = ((-r2)-r1+1)*(1-r3_2);
     h(14) = r1*(1-r3_2);
     h(15) = r2*(1-r3_2);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Pentahedron15::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                  const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune) {
     using namespace sp_grad;

     static_assert(iNumRhs >= iNumNodes);

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodes; ++i) {
               tauni(i, j) = taune(i, j);
          }
     }
}

void
Pentahedron15::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
     r(3) = ti[iNode - 1];
}

void
Tetrahedron10h::ShapeFunctionDeriv(const sp_grad::SpColVector<doublereal, 3>& r,
                                   sp_grad::SpMatrix<doublereal, iNumNodes, 3>& h0d1)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     h0d1(1,1) = 0;
     h0d1(1,2) = 4*r2-1;
     h0d1(1,3) = 0;
     h0d1(2,1) = 0;
     h0d1(2,2) = 0;
     h0d1(2,3) = 4*r3-1;
     h0d1(3,1) = 2*r3-2*((-r3)-r2-r1+1)+2*r2+2*r1-1;
     h0d1(3,2) = 2*r3-2*((-r3)-r2-r1+1)+2*r2+2*r1-1;
     h0d1(3,3) = 2*r3-2*((-r3)-r2-r1+1)+2*r2+2*r1-1;
     h0d1(4,1) = 4*r1-1;
     h0d1(4,2) = 0;
     h0d1(4,3) = 0;
     h0d1(5,1) = 0;
     h0d1(5,2) = 4*r3;
     h0d1(5,3) = 4*r2;
     h0d1(6,1) = -4*r3;
     h0d1(6,2) = -4*r3;
     h0d1(6,3) = 4*((-r3)-r2-r1+1)-4*r3;
     h0d1(7,1) = -4*r2;
     h0d1(7,2) = 4*((-r3)-r2-r1+1)-4*r2;
     h0d1(7,3) = -4*r2;
     h0d1(8,1) = 4*r2;
     h0d1(8,2) = 4*r1;
     h0d1(8,3) = 0;
     h0d1(9,1) = 4*r3;
     h0d1(9,2) = 0;
     h0d1(9,3) = 4*r1;
     h0d1(10,1) = 4*((-r3)-r2-r1+1)-4*r1;
     h0d1(10,2) = -4*r1;
     h0d1(10,3) = -4*r1;
}

void
Tetrahedron10h::ShapeFunction(const sp_grad::SpColVector<doublereal, 3>& r,
                              sp_grad::SpColVector<doublereal, iNumNodes>& h)
{
     const doublereal r1 = r(1);
     const doublereal r2 = r(2);
     const doublereal r3 = r(3);

     h(1) = r2*(2*r2-1);
     h(2) = r3*(2*r3-1);
     h(3) = ((-2*r3)-2*r2-2*r1+1)*((-r3)-r2-r1+1);
     h(4) = r1*(2*r1-1);
     h(5) = 4*r2*r3;
     h(6) = 4*((-r3)-r2-r1+1)*r3;
     h(7) = 4*r2*((-r3)-r2-r1+1);
     h(8) = 4*r1*r2;
     h(9) = 4*r1*r3;
     h(10) = 4*r1*((-r3)-r2-r1+1);
}

void
Tetrahedron10h::ShapeFunctionExtrap(const sp_grad::SpColVector<doublereal, 3>& r,
                                    sp_grad::SpColVector<doublereal, iNumNodesExtrap>& h)
{
     h(1) = r(2);
     h(2) = r(3);
     h(3) = 1. - r(1) - r(2) - r(3);
     h(4) = r(1);
}

template <sp_grad::index_type iNumComp, sp_grad::index_type iNumRhs>
void
Tetrahedron10h::GaussToNodalInterp(sp_grad::SpMatrix<doublereal, iNumNodes, iNumComp>& tauni,
                                   const sp_grad::SpMatrix<doublereal, iNumRhs, iNumComp>& taune) {
     using namespace sp_grad;

     static_assert(iNumRhs >= iNumNodesExtrap);

     static constexpr struct {
          index_type ico1, ico2, imid;
     } idxint[] = {
          {0, 1, 4},
          {1, 2, 5},
          {2, 0, 6},
          {0, 3, 7},
          {1, 3, 8},
          {2, 3, 9}
     };

     constexpr index_type iNumNodesInterp = sizeof(idxint) / sizeof(idxint[0]);

     for (index_type j = 1; j <= iNumComp; ++j) {
          for (index_type i = 1; i <= iNumNodesExtrap; ++i) {
               tauni(i, j) = taune(i, j);
          }
          for (index_type i = 0; i < iNumNodesInterp; ++i) {
               tauni(idxint[i].imid + 1, j) = 0.5 * (taune(idxint[i].ico1 + 1, j)
                                                     +  taune(idxint[i].ico2 + 1, j));
          }
     }
}

void
Tetrahedron10h::NodalPosition(sp_grad::index_type iNode, sp_grad::SpColVector<doublereal, 3>& r)
{
     ASSERT(iNode >= 1);
     ASSERT(iNode <= iNumNodes);

     r(1) = ri[iNode - 1];
     r(2) = si[iNode - 1];
     r(3) = ti[iNode - 1];
}

#endif
