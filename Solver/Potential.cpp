/***************************************************************************
*                                                                          *
*   Copyright (c) 2017                                                     *
*   FastFieldSolvers S.R.L.  http://www.fastfieldsolvers.com               *
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU Lesser General Public License (LGPL)     *
*   as published by the Free Software Foundation; either version 2 of      *
*   the License, or (at your option) any later version.                    *
*   for detail see the LICENCE text file.                                  *
*                                                                          *
*   This program is distributed in the hope that it will be useful,        *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
*   GNU Library General Public License for more details.                   *
*                                                                          *
*   You should have received a copy of the GNU Library General Public      *
*   License along with this program; if not, write to the Free Software    *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307   *
*   USA                                                                    *
*                                                                          *
***************************************************************************/


// Potential.cpp : potential class
//
// Potential computations
//
// Enrico Di Lorenzo, 2002/08/21

//#include "stdafx.h"

#include <math.h>

#include "Potential.h"

#include "FasterCapGlobal.h"

#define LOG_FOUR_PLUS_TWO   3.3862943611198906188344642429164

// these are the number of points for the n-th order
// quadrature rule and the number of weights * weigths coefficients
#define POTENTIAL_PTS_2TH_ORDER		3
#define POTENTIAL_W2_2TH_ORDER		9
#define POTENTIAL_PTS_3TH_ORDER		4
#define POTENTIAL_W2_3TH_ORDER		16
#define POTENTIAL_PTS_4TH_ORDER		6
#define POTENTIAL_W2_4TH_ORDER		36

CPotential::CPotential()
{
	// init points and weights for numerical quadrature
	InitNumerical();
}

// Initialization of integration formula points and weights:
//
//    1, NORDER =  1, precision 1, Zienkiewicz #1.
//    2, NORDER =  3, precision 2, Strang and Fix formula #1.
//    3, NORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
//    4, NORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
//    5, NORDER =  6, precision 3, Strang and Fix formula #4.
//    6, NORDER =  6, precision 3, Stroud formula T2:3-1.
//    7, NORDER =  6, precision 4, Strang and Fix formula #5.
//    8, NORDER =  7, precision 4, Strang and Fix formula #6.
//    9, NORDER =  7, precision 5, Strang and Fix formula #7,
//       Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
//   10, NORDER =  9, precision 6, Strang and Fix formula #8.
//   11, NORDER = 12, precision 6, Strang and Fix formula #9.
//   12, NORDER = 13, precision 7, Strang and Fix formula #10.
//   13, NORDER =  7, precision ?.
//   14, NORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
//   15, NORDER = 64, precision 15, triangular product Gauss rule.
//   16, NORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
//   17, NORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
//   18, NORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
//   19, NORDER = 37, precision 13, from ACM TOMS #706.
//
// Remark: 14, 15, 18, 19 not implemented at the moment
// (Best formulae seem 2, 7, 9, 10)
void CPotential::InitNumerical()
{
	double a, b, c, d, e, f, g, h;
	double p, q, r, s, t, u, v, w;
	double w1, w2, w3, w4, w5, w6, z;
	unsigned int rule, i, j;

	// matrices of integration points and weights
	// (in homogeneous coordinates)
	for(i=0; i<POTENTIAL_RULE_MAX; i++) {
		for(j=0; j<POTENTIAL_PTS_MAX; j++) {
			m_clsPoints[i][j].z = 1.0;
		}
	}

	//
	//  1 point, precision 1.
	//
	rule = 1;

	a = 1.0 / 3.0;
	w = 1.0;

	m_uiNorder[rule] = 1;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = a;
	m_dWeight[rule][0] = w;

	//
	//  3 points, precision 2, Strang and Fix formula #1.
	//
	rule = 2;

	a = 1.0;
	b = 3.0;
	c = 4.0;
	d = 6.0;

	m_uiNorder[rule] = 3;
	m_clsPoints[rule][0].x = c / d;
	m_clsPoints[rule][0].y = a / d;
	m_clsPoints[rule][1].x = a / d;
	m_clsPoints[rule][1].y = c / d;
	m_clsPoints[rule][2].x = a / d;
	m_clsPoints[rule][2].y = a / d;
	m_dWeight[rule][0] = a / b;
	m_dWeight[rule][1] = a / b;
	m_dWeight[rule][2] = a / b;

	//
	//  3 points, precision 2, Strang and Fix formula #2.
	//
	rule = 3;

	a = 0.5;
	b = 1.0;
	c = 1.0 / 3.0;
	z = 0.0;

	m_uiNorder[rule] = 3;
	m_clsPoints[rule][0].x = z;
	m_clsPoints[rule][0].y = a;
	m_clsPoints[rule][1].x = a;
	m_clsPoints[rule][1].y = z;
	m_clsPoints[rule][2].x = a;
	m_clsPoints[rule][2].y = a;
	m_dWeight[rule][0] = c;
	m_dWeight[rule][1] = c;
	m_dWeight[rule][2] = c;

	//
	//  4 points, precision 3, Strang and Fix formula #3.
	//
	rule = 4;

	a =   6.0;
	b =  10.0;
	c =  18.0;
	d =  25.0;
	e = -27.0;
	f =  30.0;
	g =  48.0;

	m_uiNorder[rule] = 4;
	m_clsPoints[rule][0].x = b / f;
	m_clsPoints[rule][0].y = b / f;
	m_clsPoints[rule][1].x = c / f;
	m_clsPoints[rule][1].y = a / f;
	m_clsPoints[rule][2].x = a / f;
	m_clsPoints[rule][2].y = c / f;
	m_clsPoints[rule][3].x = a / f;
	m_clsPoints[rule][3].y = a / f;
	m_dWeight[rule][0] = e / g;
	m_dWeight[rule][1] = d / g;
	m_dWeight[rule][2] = d / g;
	m_dWeight[rule][3] = d / g;

	//
	//  6 points, precision 3, Strang and Fix formula #4.
	//
	rule = 5;

	a = 0.659027622374092;
	b = 0.231933368553031;
	c = 0.109039009072877;
	w = 1.0 / 6.0;

	m_uiNorder[rule] = 6;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = b;
	m_clsPoints[rule][1].x = a;
	m_clsPoints[rule][1].y = c;
	m_clsPoints[rule][2].x = b;
	m_clsPoints[rule][2].y = a;
	m_clsPoints[rule][3].x = b;
	m_clsPoints[rule][3].y = c;
	m_clsPoints[rule][4].x = c;
	m_clsPoints[rule][4].y = a;
	m_clsPoints[rule][5].x = c;
	m_clsPoints[rule][5].y = b;
	m_dWeight[rule][0] = w;
	m_dWeight[rule][1] = w;
	m_dWeight[rule][2] = w;
	m_dWeight[rule][3] = w;
	m_dWeight[rule][4] = w;
	m_dWeight[rule][5] = w;

	//
	//  6 points, precision 3, Stroud T2:3-1.
	//
	rule = 6;

	a = 0.0;
	b = 0.5;
	c = 2.0 /  3.0;
	d = 1.0 /  6.0;
	v = 1.0 / 30.0;
	w = 3.0 / 10.0;

	m_uiNorder[rule] = 6;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = b;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = a;
	m_clsPoints[rule][2].x = b;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = c;
	m_clsPoints[rule][3].y = d;
	m_clsPoints[rule][4].x = d;
	m_clsPoints[rule][4].y = c;
	m_clsPoints[rule][5].x = d;
	m_clsPoints[rule][5].y = d;
	m_dWeight[rule][0] = v;
	m_dWeight[rule][1] = v;
	m_dWeight[rule][2] = v;
	m_dWeight[rule][3] = w;
	m_dWeight[rule][4] = w;
	m_dWeight[rule][5] = w;

	//
	//  6 points, precision 4, Strang and Fix, formula #5.
	//
	rule = 7;

	a = 0.816847572980459;
	b = 0.091576213509771;
	c = 0.108103018168070;
	d = 0.445948490915965;
	v = 0.109951743655322;
	w = 0.223381589678011;

	m_uiNorder[rule] = 6;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = b;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = a;
	m_clsPoints[rule][2].x = b;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = c;
	m_clsPoints[rule][3].y = d;
	m_clsPoints[rule][4].x = d;
	m_clsPoints[rule][4].y = c;
	m_clsPoints[rule][5].x = d;
	m_clsPoints[rule][5].y = d;
	m_dWeight[rule][0] = v;
	m_dWeight[rule][1] = v;
	m_dWeight[rule][2] = v;
	m_dWeight[rule][3] = w;
	m_dWeight[rule][4] = w;
	m_dWeight[rule][5] = w;

	//
	//  7 points, precision 4, Strang and Fix formula #6.
	//
	rule = 8;

	a = 1.0 / 3.0;
	c = 0.736712498968435;
	d = 0.237932366472434;
	e = 0.025355134551932;
	v = 0.375;
	w = 0.104166666666667;

	m_uiNorder[rule] = 7;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = a;
	m_clsPoints[rule][1].x = c;
	m_clsPoints[rule][1].y = d;
	m_clsPoints[rule][2].x = c;
	m_clsPoints[rule][2].y = e;
	m_clsPoints[rule][3].x = d;
	m_clsPoints[rule][3].y = c;
	m_clsPoints[rule][4].x = d;
	m_clsPoints[rule][4].y = e;
	m_clsPoints[rule][5].x = e;
	m_clsPoints[rule][5].y = c;
	m_clsPoints[rule][6].x = e;
	m_clsPoints[rule][6].y  = d;
	m_dWeight[rule][0] = v;
	m_dWeight[rule][1] = w;
	m_dWeight[rule][2] = w;
	m_dWeight[rule][3] = w;
	m_dWeight[rule][4] = w;
	m_dWeight[rule][5] = w;
	m_dWeight[rule][6] = w;

	//
	//  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
	//
	rule = 9;

	a = 1.0 / 3.0;
	b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
	c = ( 6.0 -           sqrt ( 15.0 ) ) / 21.0;
	d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
	e = ( 6.0 +           sqrt ( 15.0 ) ) / 21.0;
	u = 0.225;
	v = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
	w = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;

	m_uiNorder[rule] = 7;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = a;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = c;
	m_clsPoints[rule][2].x = c;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = c;
	m_clsPoints[rule][3].y = c;
	m_clsPoints[rule][4].x = d;
	m_clsPoints[rule][4].y = e;
	m_clsPoints[rule][5].x = e;
	m_clsPoints[rule][5].y = d;
	m_clsPoints[rule][6].x = e;
	m_clsPoints[rule][6].y  = e;
	m_dWeight[rule][0] = u;
	m_dWeight[rule][1] = v;
	m_dWeight[rule][2] = v;
	m_dWeight[rule][3] = v;
	m_dWeight[rule][4] = w;
	m_dWeight[rule][5] = w;
	m_dWeight[rule][6] = w;

	//
	//  9 points, precision 6, Strang and Fix formula #8.
	//
	rule = 10;

	a = 0.124949503233232;
	b = 0.437525248383384;
	c = 0.797112651860071;
	d = 0.165409927389841;
	e = 0.037477420750088;

	u = 0.205950504760887;
	v = 0.063691414286223;

	m_uiNorder[rule] = 9;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = b;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = a;
	m_clsPoints[rule][2].x = b;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = c;
	m_clsPoints[rule][3].y = d;
	m_clsPoints[rule][4].x = c;
	m_clsPoints[rule][4].y = e;
	m_clsPoints[rule][5].x = d;
	m_clsPoints[rule][5].y = c;
	m_clsPoints[rule][6].x = d;
	m_clsPoints[rule][6].y = e;
	m_clsPoints[rule][7].x = e;
	m_clsPoints[rule][7].y = c;
	m_clsPoints[rule][8].x = e;
	m_clsPoints[rule][8].y  = d;
	m_dWeight[rule][0] = u;
	m_dWeight[rule][1] = u;
	m_dWeight[rule][2] = u;
	m_dWeight[rule][3] = v;
	m_dWeight[rule][4] = v;
	m_dWeight[rule][5] = v;
	m_dWeight[rule][6] = v;
	m_dWeight[rule][7] = v;
	m_dWeight[rule][8] = v;


	//
	//  12 points, precision 6, Strang and Fix, formula #9.
	//
	rule = 11;

	a = 0.873821971016996;
	b = 0.063089014491502;
	c = 0.501426509658179;
	d = 0.249286745170910;
	e = 0.636502499121399;
	f = 0.310352451033785;
	g = 0.053145049844816;

	u = 0.050844906370207;
	v = 0.116786275726379;
	w = 0.082851075618374;

	m_uiNorder[rule] = 12;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = b;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = a;
	m_clsPoints[rule][2].x = b;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = d;
	m_clsPoints[rule][3].y = c;
	m_clsPoints[rule][4].x = c;
	m_clsPoints[rule][4].y = d;
	m_clsPoints[rule][5].x = d;
	m_clsPoints[rule][5].y = d;
	m_clsPoints[rule][6].x = e;
	m_clsPoints[rule][6].y = f;
	m_clsPoints[rule][7].x = e;
	m_clsPoints[rule][7].y = g;
	m_clsPoints[rule][8].x = f;
	m_clsPoints[rule][8].y = e;
	m_clsPoints[rule][9].x = f;
	m_clsPoints[rule][9].y = g;
	m_clsPoints[rule][10].x = g;
	m_clsPoints[rule][10].y = e;
	m_clsPoints[rule][11].x = g;
	m_clsPoints[rule][11].y = f;
	m_dWeight[rule][0] = u;
	m_dWeight[rule][1] = u;
	m_dWeight[rule][2] = u;
	m_dWeight[rule][3] = v;
	m_dWeight[rule][4] = v;
	m_dWeight[rule][5] = v;
	m_dWeight[rule][6] = w;
	m_dWeight[rule][7] = w;
	m_dWeight[rule][8] = w;
	m_dWeight[rule][9] = w;
	m_dWeight[rule][10] = w;
	m_dWeight[rule][11] = w;

	//
	//  13 points, precision 7, Strang and Fix, formula #10.
	//
	rule = 12;

	a = 0.479308067841923;
	b = 0.260345966079038;
	c = 0.869739794195568;
	d = 0.065130102902216;
	e = 0.638444188569809;
	f = 0.312865496004875;
	g = 0.048690315425316;
	h = 1.0 / 3.0;
	t = 0.175615257433204;
	u = 0.053347235608839;
	v = 0.077113760890257;
	w = -0.149570044467670;

	m_uiNorder[rule] = 13;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = b;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = a;
	m_clsPoints[rule][2].x = b;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = c;
	m_clsPoints[rule][3].y = d;
	m_clsPoints[rule][4].x = d;
	m_clsPoints[rule][4].y = c;
	m_clsPoints[rule][5].x = d;
	m_clsPoints[rule][5].y = d;
	m_clsPoints[rule][6].x = e;
	m_clsPoints[rule][6].y = f;
	m_clsPoints[rule][7].x = e;
	m_clsPoints[rule][7].y = g;
	m_clsPoints[rule][8].x = f;
	m_clsPoints[rule][8].y = e;
	m_clsPoints[rule][9].x = f;
	m_clsPoints[rule][9].y = g;
	m_clsPoints[rule][10].x = g;
	m_clsPoints[rule][10].y = e;
	m_clsPoints[rule][11].x = g;
	m_clsPoints[rule][11].y = f;
	m_clsPoints[rule][12].x = h;
	m_clsPoints[rule][12].y = h;
	m_dWeight[rule][0] = t;
	m_dWeight[rule][1] = t;
	m_dWeight[rule][2] = t;
	m_dWeight[rule][3] = u;
	m_dWeight[rule][4] = u;
	m_dWeight[rule][5] = u;
	m_dWeight[rule][6] = v;
	m_dWeight[rule][7] = v;
	m_dWeight[rule][8] = v;
	m_dWeight[rule][9] = v;
	m_dWeight[rule][10] = v;
	m_dWeight[rule][11] = v;
	m_dWeight[rule][12] = w;

	//
	//  7 points, precision ?.
	//
	rule = 13;

	a = 1.0 / 3.0;
	b = 1.0;
	c = 0.5;
	z = 0.0;

	u = 27.0 / 60.0;
	v =  3.0 / 60.0;
	w =  8.0 / 60.0;

	m_uiNorder[rule] = 7;
	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = a;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = z;
	m_clsPoints[rule][2].x = z;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = z;
	m_clsPoints[rule][3].y = z;
	m_clsPoints[rule][4].x = z;
	m_clsPoints[rule][4].y = c;
	m_clsPoints[rule][5].x = c;
	m_clsPoints[rule][5].y = z;
	m_clsPoints[rule][6].x = c;
	m_clsPoints[rule][6].y  = c;
	m_dWeight[rule][0] = u;
	m_dWeight[rule][1] = v;
	m_dWeight[rule][2] = v;
	m_dWeight[rule][3] = v;
	m_dWeight[rule][4] = w;
	m_dWeight[rule][5] = w;
	m_dWeight[rule][6] = w;

	//
	//  16 points.
	//
	// not implemented yet

	rule = 14;

	m_uiNorder[rule] = 0;

	/*
	norder = 16

	  norder2 = 4

		call legendre_set ( norder2, points[].x1, weight1 )

		  points[].x1(1:norder2) = 0.5 * ( points[].x1(1:norder2) + 1.0 )

			weight2(1) = 0.1355069134
			weight2(2) = 0.2034645680
			weight2(3) = 0.1298475476
			weight2(4) = 0.0311809709

			  points[].x2(1) = 0.0571041961
			  points[].x2(2) = 0.2768430136
			  points[].x2(3) = 0.5835904324
			  points[].x2(4) = 0.8602401357

				k = 0
				do i = 1, norder2
				do j = 1, norder2
				k = k + 1
				points[].x(k) = points[].x2(j)
				points[].y(k) = points[].x1(i) * ( 1.0 - points[].x2(j) )
				weight(k) = weight1(i) * weight2(j)
				end do
				end do
	*/

	//
	//  64 points, precision 15.
	//
	// not implemented yet

	rule = 15;

	m_uiNorder[rule] = 0;
	/*
	norder = 64

	  weight2(1) = 0.00329519144
	  weight2(2) = 0.01784290266
	  weight2(3) = 0.04543931950
	  weight2(4) = 0.07919959949
	  weight2(5) = 0.10604735944
	  weight2(6) = 0.11250579947
	  weight2(7) = 0.09111902364
	  weight2(8) = 0.04455080436

		points[].x2(1) = 0.04463395529
		points[].x2(2) = 0.14436625704
		points[].x2(3) = 0.28682475714
		points[].x2(4) = 0.45481331520
		points[].x2(5) = 0.62806783542
		points[].x2(6) = 0.78569152060
		points[].x2(7) = 0.90867639210
		points[].x2(8) = 0.98222008485

		  norder2 = 8
		  call legendre_set ( norder2, points[].x1, weight1 )

			k = 0
			do j = 1, norder2
			do i = 1, norder2
			k = k + 1
			points[].x(k) = 1.0 - points[].x2(j)
			points[].y(k) = 0.5 * ( 1.0 + points[].x1(i) ) * points[].x2(j)
			weight(k) = weight1(i) * weight2(j)
			end do
			end do
	*/

	//
	//  19 points, precision 8.
	//
	rule = 16;

	a = 1.0 / 3.0;
	b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
	c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
	d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
	e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
	f = ( 40.0 - 10.0 * sqrt ( 15.0 ) + 10.0 * sqrt ( 7.0 ) + 2.0 * sqrt ( 105.0 ) ) / 90.0;
	g = ( 25.0 +  5.0 * sqrt ( 15.0 ) -  5.0 * sqrt ( 7.0 ) - sqrt ( 105.0 ) ) / 90.0;
	p = ( 40.0 + 10.0 * sqrt ( 15.0 ) + 10.0 * sqrt ( 7.0 ) - 2.0 * sqrt ( 105.0 ) ) / 90.0;
	q = ( 25.0 -  5.0 * sqrt ( 15.0 ) -  5.0 * sqrt ( 7.0 ) + sqrt ( 105.0 ) ) / 90.0;
	r = ( 40.0 + 10.0 * sqrt ( 7.0 ) ) / 90.0;
	s = ( 25.0 +  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) - sqrt ( 105.0 ) ) / 90.0;
	t = ( 25.0 -  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) + sqrt ( 105.0 ) ) / 90.0;

	w1 = ( 7137.0 - 1800.0 * sqrt ( 7.0 ) ) / 62720.0;
	w2 = - 9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 ) / 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0 + 198763.0 * sqrt ( 105.0 ) / 939008.0;
	w2 = w2 / 3.0;
	w3 = -9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 ) / 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0 - 198763.0 * sqrt ( 105.0 ) / 939008.0;
	w3 = w3 / 3.0;
	w4 = ( 102791225.0 - 23876225.0 * sqrt ( 15.0 ) - 34500875.0 * sqrt ( 7.0 ) + 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0;
	w4 = w4 / 3.0;
	w5 = ( 102791225.0 + 23876225.0 * sqrt ( 15.0 ) - 34500875.0 * sqrt ( 7.0 ) - 9914825 * sqrt ( 105.0 ) ) / 59157504.0;
	w5 = w5 / 3.0;
	w6 = ( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0;
	w6 = w6 / 6.0;

	m_uiNorder[rule] = 19;

	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = a;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = c;
	m_clsPoints[rule][2].x = c;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = c;
	m_clsPoints[rule][3].y = c;
	m_clsPoints[rule][4].x = d;
	m_clsPoints[rule][4].y = e;
	m_clsPoints[rule][5].x = e;
	m_clsPoints[rule][5].y = d;
	m_clsPoints[rule][6].x = e;
	m_clsPoints[rule][6].y = e;
	m_clsPoints[rule][7].x = f;
	m_clsPoints[rule][7].y = g;
	m_clsPoints[rule][8].x = g;
	m_clsPoints[rule][8].y = f;
	m_clsPoints[rule][9].x = g;
	m_clsPoints[rule][9].y = g;
	m_clsPoints[rule][10].x = p;
	m_clsPoints[rule][10].y = q;
	m_clsPoints[rule][11].x = q;
	m_clsPoints[rule][11].y = p;
	m_clsPoints[rule][12].x = q;
	m_clsPoints[rule][12].y = q;
	m_clsPoints[rule][13].x = r;
	m_clsPoints[rule][13].y = s;
	m_clsPoints[rule][14].x = r;
	m_clsPoints[rule][14].y = t;
	m_clsPoints[rule][15].x = s;
	m_clsPoints[rule][15].y = r;
	m_clsPoints[rule][16].x = s;
	m_clsPoints[rule][16].y = t;
	m_clsPoints[rule][17].x = t;
	m_clsPoints[rule][17].y = r;
	m_clsPoints[rule][18].x = t;
	m_clsPoints[rule][18].y = s;
	m_dWeight[rule][0] = w1;
	m_dWeight[rule][1] = w2;
	m_dWeight[rule][2] = w2;
	m_dWeight[rule][3] = w2;
	m_dWeight[rule][4] = w3;
	m_dWeight[rule][5] = w3;
	m_dWeight[rule][6] = w3;
	m_dWeight[rule][7] = w4;
	m_dWeight[rule][8] = w4;
	m_dWeight[rule][9] = w4;
	m_dWeight[rule][10] = w5;
	m_dWeight[rule][11] = w5;
	m_dWeight[rule][12] = w5;
	m_dWeight[rule][13] = w6;
	m_dWeight[rule][14] = w6;
	m_dWeight[rule][15] = w6;
	m_dWeight[rule][16] = w6;
	m_dWeight[rule][17] = w6;
	m_dWeight[rule][18] = w6;


	//
	//  19 points, precision 9.
	//
	rule = 17;

	m_uiNorder[rule] = 19;

	a = 1.0 / 3.0;
	b = 0.02063496160252593;
	c = 0.4896825191987370;
	d = 0.1258208170141290;
	e = 0.4370895914929355;
	f = 0.6235929287619356;
	g = 0.1882035356190322;
	r = 0.9105409732110941;
	s = 0.04472951339445297;
	t = 0.7411985987844980;
	u = 0.03683841205473626;
	v = 0.22196288916076574;

	w1 = 0.09713579628279610;
	w2 = 0.03133470022713983;
	w3 = 0.07782754100477543;
	w4 = 0.07964773892720910;
	w5 = 0.02557767565869810;
	w6 = 0.04328353937728940;


	m_clsPoints[rule][0].x = a;
	m_clsPoints[rule][0].y = a;
	m_clsPoints[rule][1].x = b;
	m_clsPoints[rule][1].y = c;
	m_clsPoints[rule][2].x = c;
	m_clsPoints[rule][2].y = b;
	m_clsPoints[rule][3].x = c;
	m_clsPoints[rule][3].y = c;
	m_clsPoints[rule][4].x = d;
	m_clsPoints[rule][4].y = e;
	m_clsPoints[rule][5].x = e;
	m_clsPoints[rule][5].y = d;
	m_clsPoints[rule][6].x = e;
	m_clsPoints[rule][6].y = e;
	m_clsPoints[rule][7].x = f;
	m_clsPoints[rule][7].y = g;
	m_clsPoints[rule][8].x = g;
	m_clsPoints[rule][8].y = f;
	m_clsPoints[rule][9].x = g;
	m_clsPoints[rule][9].y = g;
	m_clsPoints[rule][10].x = r;
	m_clsPoints[rule][10].y = s;
	m_clsPoints[rule][11].x = s;
	m_clsPoints[rule][11].y = r;
	m_clsPoints[rule][12].x = s;
	m_clsPoints[rule][12].y = s;
	m_clsPoints[rule][13].x = t;
	m_clsPoints[rule][13].y = u;
	m_clsPoints[rule][14].x = t;
	m_clsPoints[rule][14].y = v;
	m_clsPoints[rule][15].x = u;
	m_clsPoints[rule][15].y = t;
	m_clsPoints[rule][16].x = u;
	m_clsPoints[rule][16].y = v;
	m_clsPoints[rule][17].x = v;
	m_clsPoints[rule][17].y = t;
	m_clsPoints[rule][18].x = v;
	m_clsPoints[rule][18].y = u;
	m_dWeight[rule][0] = w1;
	m_dWeight[rule][1] = w2;
	m_dWeight[rule][2] = w2;
	m_dWeight[rule][3] = w2;
	m_dWeight[rule][4] = w3;
	m_dWeight[rule][5] = w3;
	m_dWeight[rule][6] = w3;
	m_dWeight[rule][7] = w4;
	m_dWeight[rule][8] = w4;
	m_dWeight[rule][9] = w4;
	m_dWeight[rule][10] = w5;
	m_dWeight[rule][11] = w5;
	m_dWeight[rule][12] = w5;
	m_dWeight[rule][13] = w6;
	m_dWeight[rule][14] = w6;
	m_dWeight[rule][15] = w6;
	m_dWeight[rule][16] = w6;
	m_dWeight[rule][17] = w6;
	m_dWeight[rule][18] = w6;


	//
	//  28 points, precision 11.
	//
	// not implemented yet

	rule = 18;

	m_uiNorder[rule] = 0;

	/*
    a = 1.0 / 3.0
    b = 0.9480217181434233
    c = 0.02598914092828833
    d = 0.8114249947041546
    e = 0.09428750264792270
    f = 0.01072644996557060
    g = 0.4946367750172147
    p = 0.5853132347709715
    q = 0.2073433826145142
    r = 0.1221843885990187
    s = 0.4389078057004907
    t = 0.6779376548825902
    u = 0.04484167758913055
    v = 0.27722066752827925
    w = 0.8588702812826364
    x = 0.0
    y = 0.1411297187173636

    w1 = 0.08797730116222190
    w2 = 0.008744311553736190
    w3 = 0.03808157199393533
    w4 = 0.01885544805613125
    w5 = 0.07215969754474100
    w6 = 0.06932913870553720
    w7 = 0.04105631542928860
    w8 = 0.007362383783300573

    norder = 28
    points[].x(1:28) =   (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  p,  q,  q, &
       r,  s,  s,  t,  t,  u,  u,  v,  v,  w,  w,  x,  x,  y,  y /)
    points[].y(1:28) =   (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  q,  p,  q, &
       s,  r,  s,  u,  v,  t,  v,  t,  u,  x,  y,  w,  y,  w,  x /)
    weight(1:28) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
      w6, w6, w6, w7, w7, w7, w7, w7, w7, w8, w8, w8, w8, w8, w8 /)
	*/

	//
	//  37 points, precision 13.
	//
	// not implemented yet

	rule = 19;

	m_uiNorder[rule] = 0;

	/*
    a = 1.0 / 3.0
    b = 0.950275662924105565450352089520
    c = 0.024862168537947217274823955239
    d = 0.171614914923835347556304795551
    e = 0.414192542538082326221847602214
    f = 0.539412243677190440263092985511
    g = 0.230293878161404779868453507244

    w1 = 0.051739766065744133555179145422
    w2 = 0.008007799555564801597804123460
    w3 = 0.046868898981821644823226732071
    w4 = 0.046590940183976487960361770070
    w5 = 0.031016943313796381407646220131
    w6 = 0.010791612736631273623178240136
    w7 = 0.032195534242431618819414482205
    w8 = 0.015445834210701583817692900053
    w9 = 0.017822989923178661888748319485
    wx = 0.037038683681384627918546472190

    norder = 37
    points[].x(1:10) =   (/ a, b, c, c, d, e, e, f, g, g /)
    points[].y(1:10) =   (/ a, c, b, c, e, d, e, g, f, g /)
    weight(1:37) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, w5, w5, &
                      w6, w6, w6, w7, w7, w7, w8, w8, w8, w8, w8, w8, w9, &
                      w9, w9, w9, w9, w9, wx, wx, wx, wx, wx, wx /)

    a = 0.772160036676532561750285570113
    b = 0.113919981661733719124857214943

    points[].x(11) = a
    points[].y(11) = b

    points[].x(12) = b
    points[].y(12) = a

    points[].x(13) = b
    points[].y(13) = b

    a = 0.009085399949835353883572964740
    b = 0.495457300025082323058213517632

    points[].x(14) = a
    points[].y(14) = b

    points[].x(15) = b
    points[].y(15) = a

    points[].x(16) = b
    points[].y(16) = b

    a = 0.062277290305886993497083640527
    b = 0.468861354847056503251458179727

    points[].x(17) = a
    points[].y(17) = b

    points[].x(18) = b
    points[].y(18) = a

    points[].x(19) = b
    points[].y(19) = b

    a = 0.022076289653624405142446876931
    b = 0.851306504174348550389457672223
    c = 1.0 - a - b

    points[].x(20) = a
    points[].y(20) = b

    points[].x(21) = a
    points[].y(21) = c

    points[].x(22) = b
    points[].y(22) = a

    points[].x(23) = b
    points[].y(23) = c

    points[].x(24) = c
    points[].y(24) = a

    points[].x(25) = c
    points[].y(25) = b

    a = 0.018620522802520968955913511549
    b = 0.689441970728591295496647976487
    c = 1.0 - a - b

    points[].x(26) = a
    points[].y(26) = b

    points[].x(27) = a
    points[].y(27) = c

    points[].x(28) = b
    points[].y(28) = a

    points[].x(29) = b
    points[].y(29) = c

    points[].x(30) = c
    points[].y(30) = a

    points[].x(31) = c
    points[].y(31) = b

    a = 0.096506481292159228736516560903
    b = 0.635867859433872768286976979827
    c = 1.0 - a - b

    points[].x(32) = a
    points[].y(32) = b

    points[].x(33) = a
    points[].y(33) = c

    points[].x(34) = b
    points[].y(34) = a

    points[].x(35) = b
    points[].y(35) = c

    points[].x(36) = c
    points[].y(36) = a

    points[].x(37) = c
    points[].y(37) = b

*/
}

// Auto cofficient of potential of a triangular patch
// with uniform charge
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
//
// computation is fully analytic
double CPotential::Auto(C3DVector_float vertexes[3])
{
	C3DVector side1, side2, side3;
	double a, b, c, phalf, ca;

	// compute side vectors of triangle
	side1 = vertexes[1] - vertexes[0];
	side2 = vertexes[2] - vertexes[1];
	side3 = vertexes[0] - vertexes[2];

	// side lenghts
	a = Mod(side1);
	b = Mod(side2);
	c = Mod(side3);

	// semiperimeter
	phalf = (a+b+c)/2;

	// coefficient of potential
	ca = -((1/a)*log(1-a/phalf) + (1/b)*log(1-b/phalf) + (1/c)*log(1-c/phalf))/(3*PI_TIMES_E0);

	return ca;
}

// Auto cofficient of potential of a triangular patch
// with uniform charge
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
//
// computation is fully analytic
//
// Overloaded version using CAutoPanel directly
double CPotential::Auto(CAutoPanel &panel)
{
	C3DVector side1, side2, side3;
	double a, b, c, phalf, ca;

	// compute side vectors of triangle
	side1 = panel.m_clsVertex[1] - panel.m_clsVertex[0];
	side2 = panel.m_clsVertex[2] - panel.m_clsVertex[1];
	side3 = panel.m_clsVertex[0] - panel.m_clsVertex[2];

	// side lenghts
	a = Mod(side1);
	b = Mod(side2);
	c = Mod(side3);

	// semiperimeter
	phalf = (a+b+c)/2;

	// coefficient of potential
	ca = -((1/a)*log(1-a/phalf) + (1/b)*log(1-b/phalf) + (1/c)*log(1-c/phalf))/(3*PI_TIMES_E0);

	return ca;
}

double CPotential::QAuto(C3DVector vertexes[4])
{
	C3DVector_float tri1[3], tri2[3];
	double potestim1, potestim2, potestim3;
	double area;

	area = Mod(vertexes[1] - vertexes[0]) * Mod(vertexes[2] - vertexes[1]);

	tri1[0] = vertexes[0];
	tri1[1] = vertexes[1];
	tri1[2] = vertexes[2];

	tri2[0] = vertexes[2];
	tri2[1] = vertexes[0];
	tri2[2] = vertexes[3];

	potestim1 = Auto(tri1)*area*area/4.0;
	potestim2 = Auto(tri2)*area*area/4.0;
	potestim3 = Mutual_9thOrd_HalfNum(tri1, tri2, false);

	return (potestim1+potestim2+2.0*potestim3)/(area*area);
}

// Auto cofficient of potential of a triangular patch
// with uniform charge
//
// This version is half numerical, half analytical
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
// rule is the integration formula to use:
//
//    1, NORDER =  1, precision 1, Zienkiewicz #1.
//    2, NORDER =  3, precision 2, Strang and Fix formula #1.
//    3, NORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
//    4, NORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
//    5, NORDER =  6, precision 3, Strang and Fix formula #4.
//    6, NORDER =  6, precision 3, Stroud formula T2:3-1.
//    7, NORDER =  6, precision 4, Strang and Fix formula #5.
//    8, NORDER =  7, precision 4, Strang and Fix formula #6.
//    9, NORDER =  7, precision 5, Strang and Fix formula #7,
//       Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
//   10, NORDER =  9, precision 6, Strang and Fix formula #8.
//   11, NORDER = 12, precision 6, Strang and Fix formula #9.
//   12, NORDER = 13, precision 7, Strang and Fix formula #10.
//   13, NORDER =  7, precision ?.
//   14, NORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
//   15, NORDER = 64, precision 15, triangular product Gauss rule.
//   16, NORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
//   17, NORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
//   18, NORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
//   19, NORDER = 37, precision 13, from ACM TOMS #706.
//
// (Best formulae seem 2, 7, 9, 10)
//
// formula is the potential formula to use:
//
//   0, InsidePotential() (specialized)
//   1, Potential()   (generic)
//
double CPotential::AutoNumerical(C3DVector_float vertexes[3], int rule, int formula, bool divideByArea)
{
	C3DVector side1, side2, side3, normal, x, y, tz;
	C3DVector p3d;
	C3DVector vert2d[2], transform[2];
	double jacobian, pot, area1, result, point[2];
	unsigned int i;

	//
	// calculate local coordinate frame for triangle
	//

	// compute side vectors of triangle
	side1 = vertexes[1] - vertexes[0];
	side2 = vertexes[2] - vertexes[1];
	side3 = vertexes[0] - vertexes[2];

	// compute unit vector normal to triangle plane (uses vector product)
	normal = CrossProd(side1, side2);
	tz = normal/Mod(normal);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	x = side1/Mod(side1);
	// note that, since x and z are unit vectors, y is too
	y = CrossProd(x, tz);


	// calculates the coordinates of the first triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	vert2d[0].y = Mod(side1);
	vert2d[0].z = -DotProd(side3,x);
	vert2d[1].z = -DotProd(side3,y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = vert2d * unittranmtx;
	//
	// where unittranmtx is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)
	transform[0].x = vert2d[0].y-vert2d[0].z;
	transform[0].y = -vert2d[0].z;
	transform[0].z = vert2d[0].z;
	transform[1].x = -vert2d[1].z;
	transform[1].y = -vert2d[1].z;
	transform[1].z = vert2d[1].z;

	// this is the jacobian of the transform (a determinant)
	jacobian = transform[0][0]*transform[1][1]-
				transform[1][0]*transform[0][1];

	result = 0;
	for(i=0; i<m_uiNorder[rule]; i++) {
		// calculates the transformed numerical integration point
		point[0] = transform[0].x * m_clsPoints[rule][i][0] +
					transform[0].y * m_clsPoints[rule][i][1] +
					transform[0].z * m_clsPoints[rule][i][2];
		point[1] = transform[1].x * m_clsPoints[rule][i][0] +
					transform[1].y * m_clsPoints[rule][i][1] +
					transform[1].z * m_clsPoints[rule][i][2];

		// calculate 3D point location
		p3d = vertexes[0] + x*point[0] + y*point[1];
		// calculate the potential at given point due to second triangle
		if(formula == 0) {
			pot = InsidePotential(p3d, vertexes);
		}
		else {
			pot = Potential(p3d, vertexes);
		}

		// sum up the result
		result = result + m_dWeight[rule][i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result = result * jacobian / 2.0;

	// area of the triangle
	area1 = Mod(normal)/2.0;

	if( divideByArea == true )
		result /= area1;

	// return result
	return( result / FOUR_PI_TIMES_E0 );
}

// Analytical computation of the potential at a point
// inside or on the boundary of a triangle
// Uses CornerPotential()
//
// Remark: no check if the point is outside the triangle
// of a triangle due to the triangle itself
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
// point     is the evaluation point
//
double CPotential::InsidePotential(C3DVector point, C3DVector_float vertexes[3])
{
	C3DVector side1, side2, side3, join1, join2, join3;
	C3DVector_float newVertexes[3];
	double a, b, c, join1mod, join2mod, join3mod, pot1, pot2, pot3, pot, area;

	// compute side vectors of triangle
	side1 = vertexes[1] - vertexes[0];
	side2 = vertexes[2] - vertexes[1];
	side3 = vertexes[0] - vertexes[2];

	// side lenghts
	a = Mod(side1);
	b = Mod(side2);
	c = Mod(side3);

	// compute vectors joining the vertexes of the triangle
	// with the eval point
	join1 = point - vertexes[0];
	join2 = point - vertexes[1];
	join3 = point - vertexes[2];

	// join lenghts
	join1mod = Mod(join1);
	join2mod = Mod(join2);
	join3mod = Mod(join3);

	// if evaluation point is on one of the vertexes
	if(join1mod < POTENTIAL_GEO_EPS) {
		pot = CornerPotential(vertexes);
	}
	else if(join2mod < POTENTIAL_GEO_EPS) {
		newVertexes[0] = vertexes[1];
		newVertexes[1] = vertexes[2];
		newVertexes[2] = vertexes[0];
		pot = CornerPotential(newVertexes);
	}
	else if(join3mod < POTENTIAL_GEO_EPS) {
		newVertexes[0] = vertexes[2];
		newVertexes[1] = vertexes[0];
		newVertexes[2] = vertexes[1];
		pot = CornerPotential(newVertexes);
	}
	// if evaluation point is on one of the sides
	else if(fabs(CrossProd(join1, side1).z)/(a*join1mod) < POTENTIAL_GEO_EPS) {
		newVertexes[0] = point;
		newVertexes[1] = vertexes[2];
		newVertexes[2] = vertexes[0];
		pot1 = CornerPotential(newVertexes, false);
		newVertexes[0] = point;
		newVertexes[1] = vertexes[1];
		newVertexes[2] = vertexes[2];
		pot2 = CornerPotential(newVertexes, false);

		// area of triangle (norm of vector product of two sides divided by 2)
		area = Mod(CrossProd(side1,side2))/2.0;

		pot = (pot1 + pot2) / area;
	}
	else if(fabs(CrossProd(join2, side2).z)/(b*join2mod) < POTENTIAL_GEO_EPS) {
		newVertexes[0] = point;
		newVertexes[1] = vertexes[0];
		newVertexes[2] = vertexes[1];
		pot1 = CornerPotential(newVertexes, false);
		newVertexes[0] = point;
		newVertexes[1] = vertexes[2];
		newVertexes[2] = vertexes[0];
		pot2 = CornerPotential(newVertexes, false);

		// area of triangle (norm of vector product of two sides divided by 2)
		area = Mod(CrossProd(side1,side2))/2.0;

		pot = (pot1 + pot2) / area;
	}
	else if(fabs(CrossProd(join3, side3).z)/(c*join3mod) < POTENTIAL_GEO_EPS) {
		newVertexes[0] = point;
		newVertexes[1] = vertexes[1];
		newVertexes[2] = vertexes[2];
		pot1 = CornerPotential(newVertexes, false);
		newVertexes[0] = point;
		newVertexes[1] = vertexes[0];
		newVertexes[2] = vertexes[1];
		pot2 = CornerPotential(newVertexes, false);

		// area of triangle (norm of vector product of two sides divided by 2)
		area = Mod(CrossProd(side1,side2))/2.0;

		pot = (pot1 + pot2) / area;
	}
	// if inside the triangle
	else {
		newVertexes[0] = point;
		newVertexes[1] = vertexes[0];
		newVertexes[2] = vertexes[1];
		pot1 = CornerPotential(newVertexes, false);
		newVertexes[0] = point;
		newVertexes[1] = vertexes[1];
		newVertexes[2] = vertexes[2];
		pot2 = CornerPotential(newVertexes, false);
		newVertexes[0] = point;
		newVertexes[1] = vertexes[2];
		newVertexes[2] = vertexes[0];
		pot3 = CornerPotential(newVertexes, false);

		// area of triangle (norm of vector product of two sides divided by 2)
		area = Mod(CrossProd(side1,side2))/2.0;

		pot = (pot1 + pot2 + pot3) / area;
	}

	return pot;
}

// Analytical computation of the potential at the corner
// of a triangle due to the triangle itself
//
// Uses formula in 'Analytical Integration of the Fundamental Solution
// 1/R Over Panel Boundary Element', R.Rangogni, ENEL-CRIS,
// 20159 Milano, Italy, in 'Advances in Boundary Elements'
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
// the first vertex is taken as the desired evaluation point;
// computation is fully analytic
double CPotential::CornerPotential(C3DVector_float vertexes[3], bool divideByArea )
{
	C3DVector side1, side2, side3;
	double a, c, m, theta, sqr1plusm2, m_by_tghalftheta, ca;
	double area;

	// compute side vectors of triangle
	side1 = vertexes[1] - vertexes[0];
	side2 = vertexes[2] - vertexes[1];
	side3 = vertexes[0] - vertexes[2];

	// side lenghts
	a = Mod(side1);
	c = Mod(side3);

	// calculate angular coefficient of side2 w.r.t. side1
	m = DotProd(side1, side2);
	// Remark: this is only a patch for when the angle is close
	// to ninety degrees; but in this case the result is not
	// accurate !
	if(fabs(m) < POTENTIAL_GEO_EPS) {
		m = 1E+24;
	}
	else {
		m = CrossProd(side1, side2).z / m;
	}

	// calculate theta
	theta = acos( DotProd(side1, -side3) / (a*c) );

	// now the recurring terms for the formula
	sqr1plusm2 = sqrt(1+m*m);
	m_by_tghalftheta = m * tan(theta/2.0);

	// integral value
	ca = ((m*a)/sqr1plusm2) * log( ((sqr1plusm2+1+m_by_tghalftheta)*(sqr1plusm2-1)) / ((sqr1plusm2-1-m_by_tghalftheta)*(sqr1plusm2+1)) );


	if(divideByArea == true) {
		// area of triangle (norm of vector product of two sides divided by 2)
		area = Mod(CrossProd(side1,side2))/2.0;
		ca /= area;
	}

	return ca;
}

// Potential generated at an observation point
// by a triangular patch with uniform unit charge
//
// r         is the observation point (can also be
//           on boundary or inside the triangle)
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
double CPotential::Potential(C3DVector r, C3DVector_float vertexes[3])
{
	C3DVector side1, side2, side3;
	CPotParam param;
	double dnorm, Rplusnorm, Rminusnorm;
	double result, tan1, tan2, ln, nmod, nr_dot_pr;
	int i, j;

	// copy evaluation point in parameter structure
	param.r = r;

	// compute side vectors of triangle
	side1 = vertexes[1] - vertexes[0];
	side2 = vertexes[2] - vertexes[1];
	side3 = vertexes[0] - vertexes[2];

	// compute unit vector normal to triangle plane (uses vector product)
	param.n = CrossProd(side1, side2);
	// optimized operation
	//	n = n/Mod(n);
	nmod = Mod(param.n);
	param.n.x = param.n.x / nmod;
	param.n.y = param.n.y / nmod;
	param.n.z = param.n.z / nmod;

	// compute projection on triangle plane of vector to observation point (r)
	// (remember that n*r' is the dot product of n and r)
	// optimized operation
	//	param.rho = r - n * DotProd(n, r);
	nr_dot_pr = DotProd(param.n, r);
	param.rho.x = r.x - param.n.x * nr_dot_pr;
	param.rho.y = r.y - param.n.y * nr_dot_pr;
	param.rho.z = r.z - param.n.z * nr_dot_pr;

	result = 0;

	for (i=0; i<3; i++) {

		// circular scan of vertexes
		i == 2 ? j = 0 : j = i+1;
		param.vertex1 = &vertexes[i];
		param.vertex2 = &vertexes[j];

		// compute geometric values
		Sidecontrib(&param);

		// compute contribution of the current edge to potential
		// (note: test in P0 is done in order to vanish the contribute of an edge
		// in case the observation point or its projection on the triangle plane
		// lies on the edge or its extension; EPS has been chosen for the dot product
		// used to form P0 to be != 0)
		if (param.P0 > POTENTIAL_GEO_EPS) {

			dnorm = fabs(param.d);
			Rplusnorm = Mod(param.Rplus);
			Rminusnorm = Mod(param.Rminus);

			tan1 = atan(param.P0*param.lplus/(param.R0*param.R0+dnorm*Rplusnorm));
			tan2 = atan(param.P0*param.lminus/(param.R0*param.R0+dnorm*Rminusnorm));
			ln = log((Rplusnorm+param.lplus)/(Rminusnorm+param.lminus));

			// optimized operation
			//result = result + DotProd(param.P0v, param.u)*(param.P0*ln-dnorm*(tan1-tan2));
			result += (param.P0v.x*param.u.x + param.P0v.y*param.u.y + param.P0v.z*param.u.z)*(param.P0*ln-dnorm*(tan1-tan2));
		}
	}

	// divide by area of triangle (norm of vector product of two sides divided by 2)
	return (result * 2.0) / nmod;
}

// Contribution of a single side of a polygon to the potential
// generated by a uniform charge distributed on the polygon
void CPotential::Sidecontrib(CPotParam *param)
{

	C3DVector edge, rhominus, rhoplus, l, Pplus, Pminus;

	// compute projection on triangle plane of vector to first endpoint (the '-' one)
	rhominus = *(param->vertex1) - param->n * DotProd(param->n, *(param->vertex1));
	// compute projection on triangle plane of vector to second endpoint (the '+' one)
	rhoplus = *(param->vertex2) - param->n * DotProd(param->n, *(param->vertex2));

	// compute unit vector along edge
	edge = rhoplus - rhominus;
	l = edge/Mod(edge);

	// third element of the cartesian tern (cartesian tern: is this term correct?)
	param->u = CrossProd(l,param->n);

	// vectors joining observation point with end points
	param->Rminus = *(param->vertex1) - param->r;
	param->Rplus = *(param->vertex2) - param->r;

	// projections of Rplus, Rminus on polygon plane
	Pplus = rhoplus - param->rho;
	Pminus = rhominus - param->rho;

	// height of observation point above polygon plane
	// (could also be -n*Rminus')
	param->d = -DotProd(param->n, param->Rplus);

	// rectangular coordinates on polygon plane
	// (could also be fabs(Pminus*u') )
	param->P0 = fabs(DotProd(Pplus, param->u));
	param->lplus = DotProd(Pplus, l);
	param->lminus = DotProd(Pminus, l);

	// lenght of vector on plane perp. to polygon edge
	// joining obs. point with extension of edge segment
	param->R0 = sqrt(param->P0*param->P0 + param->d*param->d);

	// unit vector perpendicular to edge
	// (could also be P0v = (Pminus - lminus*l)/P0)
	// (note: test for P0 approx. null with EPS tolerance)
	if (param->P0 < POTENTIAL_GEO_EPS) {
		// dummy, not needed in fact; in this case, the contribute
		// from this edge vanishes (see main routine)
		param->P0v[0] = 1;
		param->P0v[1] = 0;
		param->P0v[2] = 0;
	}
	else
		param->P0v = (Pplus - l*param->lplus) / param->P0;

}


// Potential generated at an observation point
// by a triangular patch with uniform unit charge
//
// r         is the observation point (can also be
//           on boundary or inside the triangle)
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
// Optimized version, not calling subroutine, which is inline
double CPotential::PotentialOpt(C3DVector r, C3DVector_float vertexes[3], bool divideByArea)
{
	C3DVector side[3];
	C3DVector n, rho, rhomp[3], u, Rmp[3], P0v, Pmp[3], l;
	double dnorm, Rmp_norm[3];
	double d, P0, P0signed, P0sign, lplus, lminus, R0, R0square;
	double lplusabs, lminusabs, lmax;
	double result, tan1, tan2, ln, nmod, nr_dot_pr, nv_dot_pr;
	int i, j;

	// compute side vectors of triangle
	side[0] = vertexes[1] - vertexes[0];
	side[1] = vertexes[2] - vertexes[1];
	side[2] = vertexes[0] - vertexes[2];

	// compute unit vector normal to triangle plane (uses vector product)
	n = CrossProd(side[0], side[1]);
	// optimized operation
	//	n = n/Mod(n);
	nmod = Mod(n);
	n.x = n.x / nmod;
	n.y = n.y / nmod;
	n.z = n.z / nmod;

	// compute projection on triangle plane of vector to observation point (r)
	// (remember that n*r' is the dot product of n and r)
	// optimized operation
	//	param.rho = r - n * DotProd(n, r);
	nr_dot_pr = DotProd(n, r);
	rho.x = r.x - n.x * nr_dot_pr;
	rho.y = r.y - n.y * nr_dot_pr;
	rho.z = r.z - n.z * nr_dot_pr;

	for (i=0; i<3; i++) {
		// compute projection on triangle plane of vector to the sides endpoints
		// optimized operation
		//  rhomp[i] = vertexes[i] - n * DotProd(n, vertexes[i]);
		nv_dot_pr = DotProd(n, vertexes[i]);
		rhomp[i].x = vertexes[i].x - n.x * nv_dot_pr;
		rhomp[i].y = vertexes[i].y - n.y * nv_dot_pr;
		rhomp[i].z = vertexes[i].z - n.z * nv_dot_pr;

		// projections of Rplus, Rminus on polygon plane
		// optimized operation
		// Pmp[i] = rhomp[i] - rho;
		Pmp[i].x = rhomp[i].x - rho.x;
		Pmp[i].y = rhomp[i].y - rho.y;
		Pmp[i].z = rhomp[i].z - rho.z;

		// vectors joining observation point with end points
		// optimized operation
		//Rmp[i] = vertexes[i] - r;
		Rmp[i].x = vertexes[i].x - r.x;
		Rmp[i].y = vertexes[i].y - r.y;
		Rmp[i].z = vertexes[i].z - r.z;
		// and their modules
		Rmp_norm[i] = Mod(Rmp[i]);
	}

	// height of observation point above polygon plane
	// (could also be -n*Rmp[1] or [2]')
	d = -DotProd(n, Rmp[0]);
	dnorm = fabs(d);

	result = 0;

	for (i=0; i<3; i++) {

		// circular scan of vertexes
		i == 2 ? j = 0 : j = i+1;

		// compute geometric values (old call to Sidecontrib() )

		// compute unit vector along edge
		l = side[i]/Mod(side[i]);

		// third element of the cartesian tern (cartesian tern: is this term correct?)
		// optimized operation
		//u = CrossProd(l,n);
		C3D_CROSSPROD(u,l,n)

		// rectangular coordinates on polygon plane
		// (could also be fabs(Pminus*u') )
		P0signed = DotProd(Pmp[j], u);
		P0 = fabs(P0signed);
		lplus = DotProd(Pmp[j], l);
		lminus = DotProd(Pmp[i], l);

		// lenght of vector on plane perp. to polygon edge
		// joining obs. point with extension of edge segment
		R0 = sqrt(P0*P0 + d*d);

		// compute contribution of the current edge to potential
		// Note: in case the observation point or its projection on the triangle plane
		// lies on the edge or its extension, P0 is close to zero. However,
		// no need to test P0: if P0 is close to zero, the contribution
        // of the edge will vanish, and no numerical issues arise (no big numbers cancelling,
        // no division by zero). However, if the observation point is on the triangle plane,
        // and lies on the edge or its extension, the logarithmic term will diverge,
        // since both (Rmp_norm[j]+lplus) and Rmp_norm[i]+lminus) are close to zero.
        // However, to avoid penalizing small panels, we cannot test against a fixed threshold
        // (unless close to numerical precision of the machine). Let's weight, using the
        // greatest of 'lplus' and 'lminus' (one of the two must be different from zero)
        lplusabs = fabs(lplus);
        lminusabs = fabs(lminus);
        lplus > lminus ? lmax = lplusabs : lmax = lminusabs;
		if (P0/lmax > POTENTIAL_GEO_EPS && fabs((Rmp_norm[j]+lplus)/lmax) > POTENTIAL_GEO_EPS && fabs((Rmp_norm[i]+lminus)/lmax) > POTENTIAL_GEO_EPS) {

			R0square = R0*R0;
			tan1 = atan(P0*lplus / (R0square + dnorm*Rmp_norm[j]));
			tan2 = atan(P0*lminus / (R0square + dnorm*Rmp_norm[i]));
			ln = log((Rmp_norm[j]+lplus)/(Rmp_norm[i]+lminus));

			_ASSERT(!isnan(tan1));
			_ASSERT(!isnan(tan2));
			_ASSERT(!isnan(ln));
			_ASSERT(isfinite(tan1));
			_ASSERT(isfinite(tan2));
			_ASSERT(isfinite(ln));

			// unit vector perpendicular to edge
			// (could also be P0v = (Pminus - lminus*l)/P0)
			// (note: test for P0 approx. null with EPS tolerance)
			// optimized operation
			//P0v = (Pmp[j] - l*lplus) / P0;
			// no need to explicitly calculate P0, see below
			//P0v.x = (Pmp[j].x - l.x*lplus) / P0;
			//P0v.y = (Pmp[j].y - l.y*lplus) / P0;
			//P0v.z = (Pmp[j].z - l.z*lplus) / P0;

			// optimized operation
			//result = result + DotProd(param.P0v, param.u)*(param.P0*ln-dnorm*(tan1-tan2));
			//result += C3D_DOTPROD(P0v,u)*(P0*ln-dnorm*(tan1-tan2));
			// C3D_DOTPROD(P0v,u) is only used to define the sign of the current term, since
			// both P0v and u are unit vectors. We can understand the sign also without
			// explicitly calculate P0v. So get the sign from P0signed
			// (P0signed > 0 if true returns 1, if false returns 0) and use it
			P0sign = (P0signed > 0) - (P0signed < 0);
			result += P0sign*(P0*ln-dnorm*(tan1-tan2));
		}
	}

	if(divideByArea == true) {
		// divide by area of triangle (norm of vector product of two sides divided by 2)
		return (result * 2.0) / nmod;
	}
	else {
		return result;
	}

}

// Potential generated at an observation point
// by a triangular patch with uniform unit charge
//
// r         is the observation point (can also be
//           on boundary or inside the triangle)
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
// Optimized version, not calling subroutine, which is inline
// This overloaded version uses CAutoPanel to avoid re-calculation
// of panel parameters
double CPotential::PotentialOpt(C3DVector r, CAutoPanel panel, bool divideByArea)
{
	C3DVector side[3];
	C3DVector n, rho, rhomp[3], u, Rmp[3], P0v, Pmp[3], l;
	double dnorm, Rmp_norm[3];
	double d, P0, lplus, lminus, R0, R0square;
	double lplusabs, lminusabs, lmax;
	double result, tan1, tan2, ln, nr_dot_pr, nv_dot_pr;
	int i, j;

	// compute side vectors of triangle
	side[0] = panel.m_clsVertex[1] - panel.m_clsVertex[0];
	side[1] = panel.m_clsVertex[2] - panel.m_clsVertex[1];
	side[2] = panel.m_clsVertex[0] - panel.m_clsVertex[2];

	// compute projection on triangle plane of vector to observation point (r)
	// (remember that n*r' is the dot product of n and r)
	// optimized operation
	//	param.rho = r - n * DotProd(n, r);
	nr_dot_pr = DotProd(panel.m_clsNormal, r);
	rho.x = r.x - panel.m_clsNormal.x * nr_dot_pr;
	rho.y = r.y - panel.m_clsNormal.y * nr_dot_pr;
	rho.z = r.z - panel.m_clsNormal.z * nr_dot_pr;

	for (i=0; i<3; i++) {
		// compute projection on triangle plane of vector to the sides endpoints
		// optimized operation
		//  rhomp[i] = vertexes[i] - n * DotProd(n, vertexes[i]);
		nv_dot_pr = DotProd(panel.m_clsNormal, panel.m_clsVertex[i]);
		rhomp[i].x = panel.m_clsVertex[i].x - panel.m_clsNormal.x * nv_dot_pr;
		rhomp[i].y = panel.m_clsVertex[i].y - panel.m_clsNormal.y * nv_dot_pr;
		rhomp[i].z = panel.m_clsVertex[i].z - panel.m_clsNormal.z * nv_dot_pr;

		// projections of Rplus, Rminus on polygon plane
		// optimized operation
		// Pmp[i] = rhomp[i] - rho;
		Pmp[i].x = rhomp[i].x - rho.x;
		Pmp[i].y = rhomp[i].y - rho.y;
		Pmp[i].z = rhomp[i].z - rho.z;

		// vectors joining observation point with end points
		// optimized operation
		//Rmp[i] = vertexes[i] - r;
		Rmp[i].x = panel.m_clsVertex[i].x - r.x;
		Rmp[i].y = panel.m_clsVertex[i].y - r.y;
		Rmp[i].z = panel.m_clsVertex[i].z - r.z;
		// and their modules
		Rmp_norm[i] = Mod(Rmp[i]);
	}

	// height of observation point above polygon plane
	// (could also be -n*Rmp[1] or [2]')
	d = -DotProd(panel.m_clsNormal, Rmp[0]);
	dnorm = fabs(d);

	result = 0;

	for (i=0; i<3; i++) {

		// circular scan of vertexes
		i == 2 ? j = 0 : j = i+1;

		// compute geometric values (old call to Sidecontrib() )

		// compute unit vector along edge
		l = side[i]/Mod(side[i]);

		// third element of the cartesian tern (cartesian tern: is this term correct?)
		// optimized operation
		//u = CrossProd(l,panel.m_clsNormal);
		C3D_CROSSPROD(u,l,panel.m_clsNormal)

		// rectangular coordinates on polygon plane
		// (could also be fabs(Pminus*u') )
		P0 = fabs(DotProd(Pmp[j], u));
		lplus = DotProd(Pmp[j], l);
		// optimized operation
		lminus = DotProd(Pmp[i], l);

		// lenght of vector on plane perp. to polygon edge
		// joining obs. point with extension of edge segment
		R0 = sqrt(P0*P0 + d*d);

		// compute contribution of the current edge to potential
		// Note: in case the observation point or its projection on the triangle plane
		// lies on the edge or its extension, P0 is close to zero. However,
		// no need to test P0: if P0 is close to zero, the contribution
        // of the edge will vanish, and no numerical issues arise (no big numbers cancelling,
        // no division by zero). However, if the observation point is on the triangle plane,
        // and lies on the edge or its extension, the logarithmic term will diverge,
        // since both (Rmp_norm[j]+lplus) and Rmp_norm[i]+lminus) are close to zero.
        // However, to avoid penalizing small panels, we cannot test against a fixed threshold
        // (unless close to numerical precision of the machine). Let's weight, using the
        // greatest of 'lplus' and 'lminus' (one of the two must be different from zero)
        lplusabs = fabs(lplus);
        lminusabs = fabs(lminus);
        lplus > lminus ? lmax = lplusabs : lmax = lminusabs;
		if (P0/lmax > POTENTIAL_GEO_EPS && fabs((Rmp_norm[j]+lplus)/lmax) > POTENTIAL_GEO_EPS && fabs((Rmp_norm[i]+lminus)/lmax) > POTENTIAL_GEO_EPS) {

			R0square = R0*R0;
			tan1 = atan(P0*lplus / (R0square + dnorm*Rmp_norm[j]));
			tan2 = atan(P0*lminus / (R0square + dnorm*Rmp_norm[i]));
			ln = log((Rmp_norm[j]+lplus)/(Rmp_norm[i]+lminus));

			// unit vector perpendicular to edge
			// (could also be P0v = (Pminus - lminus*l)/P0)
			// (note: test for P0 approx. null with EPS tolerance)
			// optimized operation
			//P0v = (Pmp[j] - l*lplus) / P0;
			P0v.x = (Pmp[j].x - l.x*lplus) / P0;
			P0v.y = (Pmp[j].y - l.y*lplus) / P0;
			P0v.z = (Pmp[j].z - l.z*lplus) / P0;

			// optimized operation
			//result = result + DotProd(param.P0v, param.u)*(param.P0*ln-dnorm*(tan1-tan2));
			result += C3D_DOTPROD(P0v,u)*(P0*ln-dnorm*(tan1-tan2));
		}
	}

	if(divideByArea == true) {
		return result / panel.GetDimension();
	}
	else {
		return result;
	}
}

// Potential generated at an observation point
// by a linear segment with uniform unit charge
//
// r         is the observation point (can also be
//           inside the segment but not on the end points)
//
// vertexes  is an array of 2 2D point coordinates,
//           corresponding to the segment end points
//
// full analytic integration
// formula used is from A. R. Djordjevic, R. F. Harrington, T. K. Sarkar,
// "Evaluation of quasi-static matrix parameters for multiconductor transmission
// lines using Galerkin's method", IEEE Transactions on Microwave Theory and Techniques,
// Vol. 42, No. 7, 1994
double CPotential::PotentialOpt(C2DVector r, C2DVector_float vertexes[2], bool divideByLen)
{
	C2DVector segVector, xUnitVector, yUnitVector, origin, n, rho, ro;
	double segMod, x_plus_a, x_minus_a, add1, add2, add3, add4, result, a, y_square;
	bool done;

    done = false;

    // coordinate reference system is referred to the segment.
    // the origin is in the center of the segment,
    // x is along the segment (from vertex[0] to vertex[1]),
    // y is normal to the segment and the smallest angle
    // is in ccw direction from x

	// compute vector along the segment
	segVector = vertexes[1] - vertexes[0];
	segMod = Mod(segVector);
	ASSERT(segMod != 0.0);
	// compute x unit vector
	// optimized operation
	//	n = n/Mod(n);
	xUnitVector.x = segVector.x / segMod;
	xUnitVector.y = segVector.y / segMod;
	// compute y unit vector
	yUnitVector.x = -xUnitVector.y;
	yUnitVector.y = xUnitVector.x;
	// origin is in the middle of the segment
	origin.x = vertexes[0].x + segVector.x / 2.0;
	origin.y = vertexes[0].y + segVector.y / 2.0;
	// x-coordinate of the second end point
	a = segMod / 2.0;
	// get the field point in the new coordinate system
	rho.x = r.x - origin.x;
	rho.y = r.y - origin.y;
	ro.x = DotProd(rho, xUnitVector);
	ro.y = DotProd(rho, yUnitVector);

    // ready for calculating the integral, where (x,y) is the field point and 'a' is 1/2 of segment length
    //
    //                   2          2                       2          2
    //              log(y  + (x + a) )                 log(y  + (x - a) )
    //   - (x + a) (------------------ - 1) + (x - a) (------------------ - 1) + y atan2(y, x + a)  - y atan2(y, x - a)
    //                      2                                   2

    x_plus_a = ro.x + a;
    x_minus_a = ro.x - a;

    done = false;

   // first check singular condition
    if(fabs(ro.y) < POTENTIAL_TOL_EPS) {
        if(fabs(x_plus_a) < POTENTIAL_TOL_EPS || fabs(x_minus_a) < POTENTIAL_TOL_EPS) {
            result = -a*(log(a*a)+LOG_FOUR_PLUS_TWO);
            done = true;
        }
    }

    if(done == false) {
        y_square = ro.y * ro.y;
        add1 = -x_plus_a * (log(y_square + x_plus_a * x_plus_a) / 2.0 - 1.0);
        add2 = x_minus_a * (log(y_square + x_minus_a * x_minus_a) / 2.0 - 1.0);
        if(fabs(x_plus_a) < POTENTIAL_TOL_EPS) {
            add3 = fabs(ro.y) * PI_HALF;
        }
        else {
            add3 = ro.y * atan2(ro.y, x_plus_a);
        }
        if(fabs(x_minus_a) < POTENTIAL_TOL_EPS) {
            add4 = -fabs(ro.y) * PI_HALF;
        }
        else {
            add4 = -ro.y * atan2(ro.y, x_minus_a);
        }
        result = add1 + add2 + add3 + add4;
    }

	if(divideByLen == true) {
		// divide by length of segment
		return (result / segMod);
	}
	else {
		return result;
	}
}

// Electric field component at point 'r' in direction 'pnormal'
// due to a uniformly charged segment
//
// vertexes is an array of 2 2D point coordinates,
//           corresponding to the vertexes of the segment
//
// r is the point at which we'll calculate the electric field
//
// pnormal is the direction along which the magnitude of the electric field
//          must be calculated
//
// full analytic integration
// formula used is from A. R. Djordjevic, R. F. Harrington, T. K. Sarkar,
// "Evaluation of quasi-static matrix parameters for multiconductor transmission
// lines using Galerkin's method", IEEE Transactions on Microwave Theory and Techniques,
// Vol. 42, No. 7, 1994
double CPotential::EnField(C2DVector r, C2DVector_float vertexes[2], C2DVector pnormal, bool divideByLen)
{
	C2DVector segVector, xUnitVector, yUnitVector, origin, n, rho, ro, efield, normal;
	double segMod, x_plus_a, x_minus_a, result, a, y_square;
	bool done;

//debug
//LogMsg("r %g,%g; v1 %g,%g, v2 %g, %g; n %g %g\n", r.x, r.y, vertexes[0].x, vertexes[0].y, vertexes[1].x, vertexes[1].y, pnormal.x, pnormal.y);

    // coordinate reference system is referred to the segment.
    // the origin is in the center of the segment,
    // x is along the segment (from vertex[0] to vertex[1]),
    // y is normal to the segment and the smallest angle
    // is in ccw direction from x

	// compute vector along the segment
	segVector = vertexes[1] - vertexes[0];
	segMod = Mod(segVector);
	ASSERT(segMod != 0.0);
	// compute x unit vector
	// optimized operation
	//	n = n/Mod(n);
	xUnitVector.x = segVector.x / segMod;
	xUnitVector.y = segVector.y / segMod;
	// compute y unit vector
	yUnitVector.x = -xUnitVector.y;
	yUnitVector.y = xUnitVector.x;
	// origin is in the middle of the segment
	origin.x = vertexes[0].x + segVector.x / 2.0;
	origin.y = vertexes[0].y + segVector.y / 2.0;
	// x-coordinate of the second end point
	a = segMod / 2.0;
	// get the field point in the new coordinate system
	rho.x = r.x - origin.x;
	rho.y = r.y - origin.y;
	ro.x = DotProd(rho, xUnitVector);
	ro.y = DotProd(rho, yUnitVector);
	// get the normal in the new coordinate system
	normal.x = DotProd(pnormal, xUnitVector);
	normal.y = DotProd(pnormal, yUnitVector);
	normal.Normalize();

// debug
//LogMsg("normal %g,%g; rho %g, %g; ro %g, %g\n", normal.x, normal.y, rho.x, rho.y, ro.x, ro.y);

    // ready for calculating the integral, where (x,y) is the field point and 'a' is 1/2 of segment length.
    // The electric field has two components, corresponding to the real and imaginary part of the result
    //
    // real part:
    //                        2          2         2          2
    //                   log(y  + (x + a) )   log(y  + (x - a) )
    //                   ------------------ - ------------------
    //                            2                    2
    //
    // imaginary part:
    //
    //                    atan2(y, x - a) - atan2(y, x + a)
    //

    x_plus_a = ro.x + a;
    x_minus_a = ro.x - a;

// debug
//LogMsg("x_plus_a %g, x_minus_a %g\n", x_plus_a, x_minus_a);

    done = false;

   // first check singular condition
    if(fabs(ro.y) < POTENTIAL_TOL_EPS) {
        if(fabs(x_plus_a) < POTENTIAL_TOL_EPS || fabs(x_minus_a) < POTENTIAL_TOL_EPS) {
            result = -a*(log(a*a)+LOG_FOUR_PLUS_TWO);
            done = true;
        }
    }

    if(done == false) {
        y_square = ro.y * ro.y;
        efield.x = ( log(y_square + x_plus_a * x_plus_a) - log(y_square + x_minus_a * x_minus_a) ) / 2.0;
        efield.y = atan2(ro.y, x_minus_a) - atan2(ro.y, x_plus_a);

        result = DotProd(efield, normal);
// debug
//LogMsg("result %g; y_square %g; efield %g,%g\n", result, y_square, efield.x, efield.y);
    }


	if(divideByLen == true) {
		// divide by length of segment
		return (result / segMod);
	}
	else {
		return result;
	}
}

// Potential generated at an observation point
// by a triangular patch with uniform unit charge
//
// This version uses numerical quadrature based on rule 'rule'
// on the potential integral.
// It is not suitable for observation points very near to the
// source triangle, due to the singularity behavior of the integration
// kernel over the integration surface
//
// r         is the observation point (can also be
//           on boundary or inside the triangle)
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
double CPotential::PotentialNumerical(C3DVector r, C3DVector_float vertexes[3], int rule)
{
	C3DVector side1, side2, side3, normal, x, y, z;
	C3DVector p3d;
	double jacobian, pot, area1, result, point0, point1;
	double vert2d01, vert2d02, vert2d12;
	unsigned int i;

	// compute side vectors of triangle
	side1 = vertexes[1] - vertexes[0];
	side2 = vertexes[2] - vertexes[1];
	side3 = vertexes[0] - vertexes[2];

	// compute unit vector normal to triangle plane (uses vector product)
	// optimized operation
	//normal = CrossProd(side1, side2);
	normal.x = side1.y*side2.z - side1.z*side2.y;
	normal.y = -side1.x*side2.z + side1.z*side2.x;
	normal.z = side1.x*side2.y - side1.y*side2.x;
	z = normal/Mod(normal);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	x = side1/Mod(side1);
	// note that, since x and z are unit vectors, y is too
	// optimized operation - but seems no apparent benefit
	//y = CrossProd(x, z);
	C3D_CROSSPROD(y,x,z)

	// calculates the coordinates of the  triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	vert2d01 = Mod(side1);
	// optimized operation - but seems no apparent benefit
	//vert2d[0].z = -DotProd(side3,x);
	//vert2d[1].z = -DotProd(side3,y);
	vert2d02 = -C3D_DOTPROD(side3,x);
	vert2d12 = -C3D_DOTPROD(side3,y);

	// this is the jacobian of the transform (a determinant)
	jacobian = (vert2d02-vert2d01)*vert2d12 -
				vert2d12*vert2d02;

	result = 0;
	for(i=0; i<m_uiNorder[rule]; i++) {
		// calculates the transformed numerical integration point
		point0 = (vert2d01-vert2d02) * m_clsPoints[rule][i][0] +
					vert2d02 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][1]);
		point1 = vert2d12 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][0] - m_clsPoints[rule][i][1]);

		// calculate 3D point location
		// optimized operation
		//p3d = vertexes[0] + x*point[0] + y*point[1];
		p3d.x = vertexes[0].x + x.x*point0 + y.x*point1;
		p3d.y = vertexes[0].y + x.y*point0 + y.y*point1;
		p3d.z = vertexes[0].z + x.z*point0 + y.z*point1;

		// calculate the potential at given point due to second triangle
		pot = 1/Mod(p3d - r);
		// sum up the result
		result += m_dWeight[rule][i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result *= jacobian / 2.0;

	// area of the triangle
	area1 = Mod(normal)/2.0;
	result /= area1;

	// return result
	return( result );
}

// Electric field component in direction 'normal' at at an observation point
// due to a triangular patch with uniform charge
//
// This version uses numerical quadrature based on rule 'rule'
// on the electric field integral.
// It is not suitable for observation points very near to the
// source triangle, due to the singularity behavior of the integration
// kernel over the integration surface
//
// r         is the observation point (can also be
//           on boundary or inside the triangle)
//
// vertexes  is an array of 3 3D point coordinates,
//           corresponding to the vertexes of triangle
//
// normal    is the normal direction vector on the first triangular panel
//
double CPotential::EnFieldNumerical(C3DVector r, C3DVector_float vertexes[3], C3DVector tnormal, int rule)
{
	C3DVector side1, side2, side3, normal, x, y, z, dist;
	C3DVector p3d;
	double jacobian, pot, area1, result, point0, point1;
	double vert2d01, vert2d02, vert2d12, rdist;
	unsigned int i;

	// compute side vectors of triangle
	side1 = vertexes[1] - vertexes[0];
	side2 = vertexes[2] - vertexes[1];
	side3 = vertexes[0] - vertexes[2];

	// compute unit vector normal to triangle plane (uses vector product)
	// optimized operation
	//normal = CrossProd(side1, side2);
	normal.x = side1.y*side2.z - side1.z*side2.y;
	normal.y = -side1.x*side2.z + side1.z*side2.x;
	normal.z = side1.x*side2.y - side1.y*side2.x;
	z = normal/Mod(normal);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	x = side1/Mod(side1);
	// note that, since x and z are unit vectors, y is too
	// optimized operation - but seems no apparent benefit
	//y = CrossProd(x, z);
	C3D_CROSSPROD(y,x,z)

	// calculates the coordinates of the  triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	vert2d01 = Mod(side1);
	// optimized operation - but seems no apparent benefit
	//vert2d[0].z = -DotProd(side3,x);
	//vert2d[1].z = -DotProd(side3,y);
	vert2d02 = -C3D_DOTPROD(side3,x);
	vert2d12 = -C3D_DOTPROD(side3,y);

	// this is the jacobian of the transform (a determinant)
	jacobian = (vert2d02-vert2d01)*vert2d12 -
				vert2d12*vert2d02;

	result = 0;
	for(i=0; i<m_uiNorder[rule]; i++) {
		// calculates the transformed numerical integration point
		point0 = (vert2d01-vert2d02) * m_clsPoints[rule][i][0] +
					vert2d02 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][1]);
		point1 = vert2d12 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][0] - m_clsPoints[rule][i][1]);

		// calculate 3D point location
		// optimized operation
		//p3d = vertexes[0] + x*point[0] + y*point[1];
		p3d.x = vertexes[0].x + x.x*point0 + y.x*point1;
		p3d.y = vertexes[0].y + x.y*point0 + y.y*point1;
		p3d.z = vertexes[0].z + x.z*point0 + y.z*point1;


		// calculate the electric field component in direction
		// 'normal' at given point due to the triangle patch
		dist = r - p3d;
		rdist = Mod(dist);
		// this is cos((x1-x2),n1) * mod(x1-x2)
		pot = DotProd(dist, tnormal) / (rdist * rdist * rdist);

		// sum up the result
		result += m_dWeight[rule][i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result *= jacobian / 2.0;

	// divide by the area of the triangle
	area1 = Mod(normal)/2.0;
	result /= area1;

	// return result
	return( result );
}

// Mutual cofficient of potential between two triangular patches
// with uniform charge
//
// This version uses numerical quadrature based on Stroud
// formula 19 points, 9th order
// on the outer integral and analytical integration of the potential
// on the inner integral (also to deal with singularity)
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
double CPotential::Mutual_9thOrd_HalfNum(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], bool divideByArea)
{
	C3DVector side1, side2, side3, normal, x, y, z;
	C3DVector p3d;
	C3DVector vert2d[2], transform[2];
//	CLin_Matrix vert2d(2,3), transform;
	// this is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
//	static const double unittranmtx[9] = { 0, 1, 0,
//	                                      1, 0, 0,
//							             -1,-1, 1};
//	CLin_Matrix unittranmtx(3, 3, tranmatrix);
	double jacobian, pot, area1, area2, result, point[2];
	int i;

	// compute side vectors of triangle 1
	side1 = vertexes1[1] - vertexes1[0];
	side2 = vertexes1[2] - vertexes1[1];
	side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	normal = CrossProd(side1, side2);
	z = normal/Mod(normal);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	x = side1/Mod(side1);
	// note that, since x and z are unit vectors, y is too
	y = CrossProd(x, z);

	// numerical 2D quadrature over a triangle
	// uses Stroud 19 points, 9th order
	//
	// considers an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration
	// triangular domain

	// initialize matrices of integration points and weights
	// (in homogeneous coordinates)
	static const double points[19][3] = {
		{0.33333333333333333, 0.33333333333333333, 1.0},
		{0.02063496160252593, 0.4896825191987370, 1.0},
		{0.4896825191987370,  0.02063496160252593, 1.0},
		{0.4896825191987370,  0.4896825191987370, 1.0},
		{0.1258208170141290,  0.4370895914929355, 1.0},
		{0.4370895914929355,  0.1258208170141290, 1.0},
		{0.4370895914929355,  0.4370895914929355, 1.0},
		{0.6235929287619356,  0.1882035356190322, 1.0},
		{0.1882035356190322,  0.6235929287619356, 1.0},
		{0.1882035356190322,  0.1882035356190322, 1.0},
		{0.9105409732110941,  0.04472951339445297, 1.0},
		{0.04472951339445297, 0.9105409732110941, 1.0},
		{0.04472951339445297, 0.04472951339445297, 1.0},
		{0.7411985987844980,  0.03683841205473626, 1.0},
		{0.7411985987844980,  0.22196288916076574, 1.0},
		{0.03683841205473626, 0.7411985987844980, 1.0},
		{0.03683841205473626, 0.22196288916076574, 1.0},
		{0.22196288916076574, 0.7411985987844980, 1.0},
		{0.22196288916076574, 0.03683841205473626, 1.0}};

	static const double weights[19] =
		{0.09713579628279610, 0.03133470022713983, 0.03133470022713983,
		0.03133470022713983, 0.07782754100477543, 0.07782754100477543,
		0.07782754100477543, 0.07964773892720910, 0.07964773892720910,
		0.07964773892720910, 0.02557767565869810, 0.02557767565869810,
		0.02557767565869810, 0.04328353937728940, 0.04328353937728940,
		0.04328353937728940, 0.04328353937728940, 0.04328353937728940,
		0.04328353937728940};

	// calculates the coordinates of the first triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	vert2d[0].y = Mod(side1);
	vert2d[0].z = -DotProd(side3,x);
	vert2d[1].z = -DotProd(side3,y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = vert2d * unittranmtx;
	//
	// where unittranmtx is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)
	transform[0].x = vert2d[0].y-vert2d[0].z;
	transform[0].y = -vert2d[0].z;
	transform[0].z = vert2d[0].z;
	transform[1].x = -vert2d[1].z;
	transform[1].y = -vert2d[1].z;
	transform[1].z = vert2d[1].z;

	// this is the jacobian of the transform (a determinant)
	jacobian = transform[0][0]*transform[1][1]-
				transform[1][0]*transform[0][1];

	result = 0;
	for(i=0; i<19; i++) {
		// calculates the transformed numerical integration point
		point[0] = transform[0].x * points[i][0] +
					transform[0].y * points[i][1] +
					transform[0].z * points[i][2];
		point[1] = transform[1].x * points[i][0] +
					transform[1].y * points[i][1] +
					transform[1].z * points[i][2];

		// calculate 3D point location
		p3d = vertexes1[0] + x*point[0] + y*point[1];
		// calculate the potential at given point due to second triangle
		pot = Potential(p3d, vertexes2);
		// sum up the result
		result = result + weights[i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result = result * jacobian / 2.0;

	// area of the first triangle
	area1 = Mod(normal)/2.0;

	// potential calculated by Potential() is already considering
	// division by area of source panel
	if( divideByArea == true ) {
		result /= area1;
	}
	else {
		// area of second triangle
		//
		// compute side vectors of triangle 2
		side1 = vertexes2[1] - vertexes2[0];
		side2 = vertexes2[2] - vertexes2[1];
		// area (norm of vector product of two sides divided by 2)
		area2 = Mod(CrossProd(side1,side2))/2.0;

		result *= area2;
	}

	// return result
	return( result/FOUR_PI_TIMES_E0 );
}

double CPotential::QMutual_9thOrd_HalfNum(C3DVector vertexes1[4],
						    C3DVector vertexes2[4], bool divideByArea)
{
	C3DVector_float tri11[3], tri12[3], tri21[3], tri22[3];
	double potestim1, potestim2, potestim3, potestim4;
	double area1, area2;

	tri11[0] = vertexes1[0];
	tri11[1] = vertexes1[1];
	tri11[2] = vertexes1[2];

	tri12[0] = vertexes1[2];
	tri12[1] = vertexes1[0];
	tri12[2] = vertexes1[3];

	tri21[0] = vertexes2[0];
	tri21[1] = vertexes2[1];
	tri21[2] = vertexes2[2];

	tri22[0] = vertexes2[2];
	tri22[1] = vertexes2[0];
	tri22[2] = vertexes2[3];

	potestim1 = Mutual_9thOrd_HalfNum(tri11, tri21, false);
	potestim2 = Mutual_9thOrd_HalfNum(tri11, tri22, false);
	potestim3 = Mutual_9thOrd_HalfNum(tri12, tri21, false);
	potestim4 = Mutual_9thOrd_HalfNum(tri12, tri22, false);

	if(divideByArea == true) {
		area1 = Mod(vertexes1[1] - vertexes1[0]) * Mod(vertexes1[2] - vertexes1[1]);
		area2 = Mod(vertexes2[1] - vertexes2[0]) * Mod(vertexes2[2] - vertexes2[1]);

		return (potestim1+potestim2+potestim3+potestim4)/(area1*area2);
	}
	else {
		return (potestim1+potestim2+potestim3+potestim4);
	}
}

// Mutual cofficient of potential between two triangular patches
// with uniform charge
//
// This version uses numerical quadrature based on rule 'rule'
// on the outer integral and analytical integration of the potential
// on the inner integral (also to deal with singularity)
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
double CPotential::MutualHalfNumerical(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], int rule, bool divideByArea)
{
	C3DVector side1, side2, side3, normal, x, y, z;
	C3DVector p3d;
	//C3DVector vert2d[2], transform[2];
	//CLin_Matrix vert2d(2,3), transform;
	// this is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//static const double unittranmtx[9] = { 0, 1, 0,
	//                                      1, 0, 0,
	//						             -1,-1, 1};
	//CLin_Matrix unittranmtx(3, 3, tranmatrix);
	double jacobian, pot, area1, area2, result, point0, point1;
	double vert2d01, vert2d02, vert2d12;
	unsigned int i;

	// compute side vectors of triangle 1
	side1 = vertexes1[1] - vertexes1[0];
	side2 = vertexes1[2] - vertexes1[1];
	side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	// optimized operation
	//normal = CrossProd(side1, side2);
	normal.x = side1.y*side2.z - side1.z*side2.y;
	normal.y = -side1.x*side2.z + side1.z*side2.x;
	normal.z = side1.x*side2.y - side1.y*side2.x;
	z = normal/Mod(normal);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	x = side1/Mod(side1);
	// note that, since x and z are unit vectors, y is too
	// optimized operation - but seems no apparent benefit
	//y = CrossProd(x, z);
	C3D_CROSSPROD(y,x,z)

	// calculates the coordinates of the first triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	vert2d01 = Mod(side1);
	// optimized operation - but seems no apparent benefit
	//vert2d[0].z = -DotProd(side3,x);
	//vert2d[1].z = -DotProd(side3,y);
	vert2d02 = -C3D_DOTPROD(side3,x);
	vert2d12 = -C3D_DOTPROD(side3,y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = vert2d * unittranmtx;
	//
	// where unittranmtx is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)
	//transform[0].x = vert2d[0].y-vert2d[0].z;
	//transform[0].y = -vert2d[0].z;
	//transform[0].z = vert2d[0].z;
	//transform[1].x = -vert2d[1].z;
	//transform[1].y = -vert2d[1].z;
	//transform[1].z = vert2d[1].z;

	// this is the jacobian of the transform (a determinant)
	//jacobian = transform[0][0]*transform[1][1]-
	//			transform[1][0]*transform[0][1];
	jacobian = (vert2d02-vert2d01)*vert2d12 -
				vert2d12*vert2d02;

	result = 0;
	for(i=0; i<m_uiNorder[rule]; i++) {
		// calculates the transformed numerical integration point
		//point[0] = transform[0].x * m_clsPoints[rule][i][0] +
		//			transform[0].y * m_clsPoints[rule][i][1] +
		//			transform[0].z * m_clsPoints[rule][i][2];
		//point[1] = transform[1].x * m_clsPoints[rule][i][0] +
		//			transform[1].y * m_clsPoints[rule][i][1] +
		//			transform[1].z * m_clsPoints[rule][i][2];
		point0 = (vert2d01-vert2d02) * m_clsPoints[rule][i][0] +
					vert2d02 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][1]);
		point1 = vert2d12 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][0] - m_clsPoints[rule][i][1]);

		// calculate 3D point location
		// optimized operation
		//p3d = vertexes1[0] + x*point[0] + y*point[1];
		p3d.x = vertexes1[0].x + x.x*point0 + y.x*point1;
		p3d.y = vertexes1[0].y + x.y*point0 + y.y*point1;
		p3d.z = vertexes1[0].z + x.z*point0 + y.z*point1;

		// calculate the potential at given point due to second triangle
		pot = PotentialOpt(p3d, vertexes2);
		// sum up the result
		result += m_dWeight[rule][i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result *= jacobian / 2.0;

	// area of the first triangle
	area1 = Mod(normal)/2.0;

	// potential calculated by PotentialOpt() is already considering
	// division by area of source panel
	if( divideByArea == true ) {
		result /= area1;
	}
	else {
		// area of second triangle
		//
		// compute side vectors of triangle 2
		side1 = vertexes2[1] - vertexes2[0];
		side2 = vertexes2[2] - vertexes2[1];
		// area (norm of vector product of two sides divided by 2)
		area2 = Mod(CrossProd(side1,side2))/2.0;

		result *= area2;
	}

	// return result
	return( result/FOUR_PI_TIMES_E0 );
}

// Mutual cofficient of potential between two triangular patches
// with uniform charge
//
// This version uses numerical quadrature based on rule 'rule'
// on the outer integral and analytical integration of the potential
// on the inner integral (also to deal with singularity)
//
// this is an optimized version of the MutualHalfNumerical() routine,
// using also the quantities already calculated as CAutoPanel
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
// Overloaded version passing panels as arguments, instead of vertexes
double CPotential::MutualHalfNumerical(CAutoPanel panel1,
						  CAutoPanel panel2, int rule, bool divideByArea)
{
	C3DVector side1, side2, side3, x, y;
	C3DVector p3d;
	//C3DVector vert2d[2], transform[2];
	//CLin_Matrix vert2d(2,3), transform;
	// this is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//static const double unittranmtx[9] = { 0, 1, 0,
	//                                      1, 0, 0,
	//						             -1,-1, 1};
	//CLin_Matrix unittranmtx(3, 3, tranmatrix);
	double jacobian, pot, result, point0, point1;
	double vert2d01, vert2d02, vert2d12;
	unsigned int i;

	// compute side vectors of triangle 1
	side1 = panel1.m_clsVertex[1] - panel1.m_clsVertex[0];
	side3 = panel1.m_clsVertex[0] - panel1.m_clsVertex[2];

	// sidelen1, also used below as coordinate of the first triangle in local 2D frame
	vert2d01 = Mod(side1);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
//	x = side1 / panel1.m_dSideLen[0];
	x = side1 / vert2d01;
	// note that, since x and z are unit vectors, y is too
	// optimized operation - but seems no apparent benefit
	//y = CrossProd(x, panel1.m_clsNormal);
	C3D_CROSSPROD(y,x,panel1.m_clsNormal)

	// calculates the coordinates of the first triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	//
	// old assignment, moved above
	//vert2d01 = panel1.m_dSideLen[0];
	// optimized operation - but seems no apparent benefit
	//vert2d[0].z = -DotProd(side3,x);
	//vert2d[1].z = -DotProd(side3,y);
	vert2d02 = -C3D_DOTPROD(side3,x);
	vert2d12 = -C3D_DOTPROD(side3,y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = vert2d * unittranmtx;
	//
	// where unittranmtx is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)
	//transform[0].x = vert2d[0].y-vert2d[0].z;
	//transform[0].y = -vert2d[0].z;
	//transform[0].z = vert2d[0].z;
	//transform[1].x = -vert2d[1].z;
	//transform[1].y = -vert2d[1].z;
	//transform[1].z = vert2d[1].z;

	// this is the jacobian of the transform (a determinant)
	//jacobian = transform[0][0]*transform[1][1]-
	//			transform[1][0]*transform[0][1];
	jacobian = (vert2d02-vert2d01)*vert2d12 -
				vert2d12*vert2d02;

	result = 0;
	for(i=0; i<m_uiNorder[rule]; i++) {
		// calculates the transformed numerical integration point
		// calculates the transformed numerical integration point
		//point[0] = transform[0].x * m_clsPoints[rule][i][0] +
		//			transform[0].y * m_clsPoints[rule][i][1] +
		//			transform[0].z * m_clsPoints[rule][i][2];
		//point[1] = transform[1].x * m_clsPoints[rule][i][0] +
		//			transform[1].y * m_clsPoints[rule][i][1] +
		//			transform[1].z * m_clsPoints[rule][i][2];
		point0 = (vert2d01-vert2d02) * m_clsPoints[rule][i][0] +
					vert2d02 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][1]);
		point1 = vert2d12 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][0] - m_clsPoints[rule][i][1]);

		// calculate 3D point location
		// optimized operation
		//p3d = panel1.m_clsVertex[0] + x*point[0] + y*point[1];
		p3d.x = panel1.m_clsVertex[0].x + x.x*point0 + y.x*point1;
		p3d.y = panel1.m_clsVertex[0].y + x.y*point0 + y.y*point1;
		p3d.z = panel1.m_clsVertex[0].z + x.z*point0 + y.z*point1;

		// calculate the potential at given point due to second triangle
		pot = PotentialOpt(p3d, panel2);
		// sum up the result
		result += m_dWeight[rule][i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result *= jacobian / 2.0;

	// potential calculated by PotentialOpt() is already considering
	// division by area of source panel
	if( divideByArea == true ) {
		result /= panel1.GetDimension();
	}
	else {
		result *= panel2.GetDimension();
	}

	// return result
	return( result/FOUR_PI_TIMES_E0 );
}

// Electric field component in direction normal to one triangular patch
// due to a second triangular patch with uniform charge
//
// This version uses numerical quadrature based on rule 'rule'
// on the outer integral and divided differences of two analytical
// integration on the potential on the inner integral,
// to deal with singularity
//
// this is an optimized version of the MutualHalfNumerical() routine,
// using also the quantities already calculated as CAutoPanel
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
// normal is the normal direction vector on the first triangular panel
//
// h is the distance of the observation points above and below the
// first panel, along the normal, for usage in the divided differences
//
// rule is the integration formula to use, see InitNumerical()
double CPotential::MutualDHalfNumerical(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], C3DVector pnormal, double h, int rule, bool divideByArea)
{
	C3DVector side1, side2, side3, normal, x, y, z;
	C3DVector p3d, evalPlus, evalMinus;
	double potestim_p, potestim_n;
	//C3DVector vert2d[2], transform[2];
	//CLin_Matrix vert2d(2,3), transform;
	// this is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//static const double unittranmtx[9] = { 0, 1, 0,
	//                                      1, 0, 0,
	//						             -1,-1, 1};
	//CLin_Matrix unittranmtx(3, 3, tranmatrix);
	double jacobian, pot, area1, area2, result, point0, point1;
	double vert2d01, vert2d02, vert2d12;
	unsigned int i;

	// compute side vectors of triangle 1
	side1 = vertexes1[1] - vertexes1[0];
	side2 = vertexes1[2] - vertexes1[1];
	side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	// optimized operation
	//normal = CrossProd(side1, side2);
	normal.x = side1.y*side2.z - side1.z*side2.y;
	normal.y = -side1.x*side2.z + side1.z*side2.x;
	normal.z = side1.x*side2.y - side1.y*side2.x;
	z = normal/Mod(normal);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	x = side1/Mod(side1);
	// note that, since x and z are unit vectors, y is too
	// optimized operation - but seems no apparent benefit
	//y = CrossProd(x, z);
	C3D_CROSSPROD(y,x,z)

	// calculates the coordinates of the first triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	vert2d01 = Mod(side1);
	// optimized operation - but seems no apparent benefit
	//vert2d[0].z = -DotProd(side3,x);
	//vert2d[1].z = -DotProd(side3,y);
	vert2d02 = -C3D_DOTPROD(side3,x);
	vert2d12 = -C3D_DOTPROD(side3,y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = vert2d * unittranmtx;
	//
	// where unittranmtx is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)
	//transform[0].x = vert2d[0].y-vert2d[0].z;
	//transform[0].y = -vert2d[0].z;
	//transform[0].z = vert2d[0].z;
	//transform[1].x = -vert2d[1].z;
	//transform[1].y = -vert2d[1].z;
	//transform[1].z = vert2d[1].z;

	// this is the jacobian of the transform (a determinant)
	//jacobian = transform[0][0]*transform[1][1]-
	//			transform[1][0]*transform[0][1];
	jacobian = (vert2d02-vert2d01)*vert2d12 -
				vert2d12*vert2d02;

	result = 0;
	for(i=0; i<m_uiNorder[rule]; i++) {
		// calculates the transformed numerical integration point
		//point[0] = transform[0].x * m_clsPoints[rule][i][0] +
		//			transform[0].y * m_clsPoints[rule][i][1] +
		//			transform[0].z * m_clsPoints[rule][i][2];
		//point[1] = transform[1].x * m_clsPoints[rule][i][0] +
		//			transform[1].y * m_clsPoints[rule][i][1] +
		//			transform[1].z * m_clsPoints[rule][i][2];
		point0 = (vert2d01-vert2d02) * m_clsPoints[rule][i][0] +
					vert2d02 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][1]);
		point1 = vert2d12 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][0] - m_clsPoints[rule][i][1]);

		// calculate 3D point location
		// optimized operation
		//p3d = vertexes1[0] + x*point[0] + y*point[1];
		p3d.x = vertexes1[0].x + x.x*point0 + y.x*point1;
		p3d.y = vertexes1[0].y + x.y*point0 + y.y*point1;
		p3d.z = vertexes1[0].z + x.z*point0 + y.z*point1;

		// calculate electric field at given point due to second triangle,
		// using the divided differences of two potential evaluations

		// evaluation points above and below the target point
		evalPlus = p3d + pnormal*h;
		evalMinus = p3d - pnormal*h;
		// calculate the potentials
		potestim_p = PotentialOpt(evalPlus, vertexes2);
		potestim_n = PotentialOpt(evalMinus, vertexes2);
		// and get the electric field
		pot = (potestim_n - potestim_p) / (2.0*h);

		// sum up the result
		result += m_dWeight[rule][i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result *= jacobian / 2.0;

	// area of the first triangle
	area1 = Mod(normal)/2.0;

	// potential calculated by PotentialOpt() is already considering
	// division by area of source panel
	if( divideByArea == true ) {
		result /= area1;
	}
	else {
		// area of second triangle
		//
		// compute side vectors of triangle 2
		side1 = vertexes2[1] - vertexes2[0];
		side2 = vertexes2[2] - vertexes2[1];
		// area (norm of vector product of two sides divided by 2)
		area2 = Mod(CrossProd(side1,side2))/2.0;

		result *= area2;
	}

	// return result
	return( result/FOUR_PI_TIMES_E0 );
}

// Electric field component in direction normal to one triangular patch
// due to a second triangular patch with uniform charge
//
// This version uses numerical quadrature based on rule 'rule'
// on the outer integral and divided differences of two analytical
// integration on the potential on the inner integral,
// to deal with singularity
//
// this is an optimized version of the MutualHalfNumerical() routine,
// using also the quantities already calculated as CAutoPanel
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
// normal is the normal direction vector on the first triangular panel
//
// h is the distance of the observation points above and below the
// first panel, along the normal, for usage in the divided differences
//
// rule is the integration formula to use, see InitNumerical()
//
// Overloaded version passing panels as arguments, instead of vertexes
double CPotential::MutualDHalfNumerical(CAutoPanel panel1,
						  CAutoPanel panel2, C3DVector pnormal, double h, int rule, bool divideByArea)
{
	C3DVector side1, side2, side3, x, y;
	C3DVector p3d, evalPlus, evalMinus;
	double potestim_p, potestim_n;
	//C3DVector vert2d[2], transform[2];
	//CLin_Matrix vert2d(2,3), transform;
	// this is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//static const double unittranmtx[9] = { 0, 1, 0,
	//                                      1, 0, 0,
	//						             -1,-1, 1};
	//CLin_Matrix unittranmtx(3, 3, tranmatrix);
	double jacobian, pot, result, point0, point1;
	double vert2d01, vert2d02, vert2d12;
	unsigned int i;

	// compute side vectors of triangle 1
	side1 = panel1.m_clsVertex[1] - panel1.m_clsVertex[0];
	side3 = panel1.m_clsVertex[0] - panel1.m_clsVertex[2];

	// sidelen1, also used below as coordinate of the first triangle in local 2D frame
	vert2d01 = Mod(side1);

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
//	x = side1 / panel1.m_dSideLen[0];
	x = side1 / vert2d01;
	// note that, since x and z are unit vectors, y is too
	// optimized operation - but seems no apparent benefit
	//y = CrossProd(x, panel1.m_clsNormal);
	C3D_CROSSPROD(y,x,panel1.m_clsNormal)

	// calculates the coordinates of the first triangle in local 2D frame;
	// the origin is located on first vertex.
	// Therefore vert2d is an array of 3 2D point coordinates,
	// corresponding to the vertexes of the intergration triangular domain, [x1,y1;x2,y2;x3,y3]
	//
	// old assignment, moved above
	//vert2d01 = panel1.m_dSideLen[0];
	// optimized operation - but seems no apparent benefit
	//vert2d[0].z = -DotProd(side3,x);
	//vert2d[1].z = -DotProd(side3,y);
	vert2d02 = -C3D_DOTPROD(side3,x);
	vert2d12 = -C3D_DOTPROD(side3,y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = vert2d * unittranmtx;
	//
	// where unittranmtx is the transpose of the inverse of [0,1,1;1,0,1;0,0,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)
	//transform[0].x = vert2d[0].y-vert2d[0].z;
	//transform[0].y = -vert2d[0].z;
	//transform[0].z = vert2d[0].z;
	//transform[1].x = -vert2d[1].z;
	//transform[1].y = -vert2d[1].z;
	//transform[1].z = vert2d[1].z;

	// this is the jacobian of the transform (a determinant)
	//jacobian = transform[0][0]*transform[1][1]-
	//			transform[1][0]*transform[0][1];
	jacobian = (vert2d02-vert2d01)*vert2d12 -
				vert2d12*vert2d02;

	result = 0;
	for(i=0; i<m_uiNorder[rule]; i++) {
		// calculates the transformed numerical integration point
		// calculates the transformed numerical integration point
		//point[0] = transform[0].x * m_clsPoints[rule][i][0] +
		//			transform[0].y * m_clsPoints[rule][i][1] +
		//			transform[0].z * m_clsPoints[rule][i][2];
		//point[1] = transform[1].x * m_clsPoints[rule][i][0] +
		//			transform[1].y * m_clsPoints[rule][i][1] +
		//			transform[1].z * m_clsPoints[rule][i][2];
		point0 = (vert2d01-vert2d02) * m_clsPoints[rule][i][0] +
					vert2d02 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][1]);
		point1 = vert2d12 * (m_clsPoints[rule][i][2] - m_clsPoints[rule][i][0] - m_clsPoints[rule][i][1]);

		// calculate 3D point location
		// optimized operation
		//p3d = panel1.m_clsVertex[0] + x*point[0] + y*point[1];
		p3d.x = panel1.m_clsVertex[0].x + x.x*point0 + y.x*point1;
		p3d.y = panel1.m_clsVertex[0].y + x.y*point0 + y.y*point1;
		p3d.z = panel1.m_clsVertex[0].z + x.z*point0 + y.z*point1;

		// calculate electric field at given point due to second triangle,
		// using the divided differences of two potential evaluations

		// evaluation points above and below the target point
		evalPlus = p3d + pnormal*h;
		evalMinus = p3d - pnormal*h;
		// calculate the potentials
		potestim_p = PotentialOpt(evalPlus, panel2);
		potestim_n = PotentialOpt(evalMinus, panel2);
		// and get the electric field
		pot = (potestim_n - potestim_p) / (2.0*h);

		// sum up the result
		result += m_dWeight[rule][i]*pot;
	}

	// do not forget the jacobian!
	// note that we must divide by the unit triangle area (could be included
	// in the weights; it is a side effect of how, in Stroud, the transformed
	// coordinates are calculated)
	result *= jacobian / 2.0;

	// potential calculated by PotentialOpt() is already considering
	// division by area of source panel
	if( divideByArea == true ) {
		result /= panel1.GetDimension();
	}
	else {
		result *= panel2.GetDimension();
	}

	// return result
	return( result/FOUR_PI_TIMES_E0 );
}

// Mutual cofficient of potential between two triangular patches
// with uniform charge
//
// This version is fully numerical and uses Strang and Fix, formula #1,
// 3 points, 2th order, for numerical quadrature (from Stroud)
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
double CPotential::Mutual_2thOrd_FullNum(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], bool divideByArea)
{
	C3DVector t1side1, t1side2, t1side3, t1normal, t1x, t1y, t1z;
	C3DVector t2side1, t2side2, t2side3, t2normal, t2x, t2y, t2z;
	C2DVector t1local2D[3], t2local2D[3];
	C3DVector t1p3d[POTENTIAL_PTS_2TH_ORDER], t2p3d[POTENTIAL_PTS_2TH_ORDER];
	C3DVector t1transform[2], t2transform[2];
	double mod_t1side1, mod_t2side1, mod_t1normal, mod_t2normal;
	double t1jacobian, t2jacobian, pot, result, point_x, point_y;
	int i, j, k;

	// numerical 2D quadrature over a triangle
	// uses Stroud 3 points, 2th order
	//
	// initialize matrices of integration points and weights
	// (in homogeneous coordinates)
	static const double points[POTENTIAL_PTS_2TH_ORDER][3] = {
		{0.66666666666666666666666666666667, 0.16666666666666666666666666666667, 1.0},
		{0.16666666666666666666666666666667, 0.66666666666666666666666666666667, 1.0},
		{0.16666666666666666666666666666667, 0.16666666666666666666666666666667, 1.0}};

	static const double weights2[POTENTIAL_W2_2TH_ORDER] =
		{ 0.11111111111111111111111111111111, 0.11111111111111111111111111111111, 0.11111111111111111111111111111111,
		  0.11111111111111111111111111111111, 0.11111111111111111111111111111111, 0.11111111111111111111111111111111,
		  0.11111111111111111111111111111111, 0.11111111111111111111111111111111, 0.11111111111111111111111111111111};

	//
	// calculate local coordinate frame for triangle 1
	//

	// compute side vectors of triangle
	t1side1 = vertexes1[1] - vertexes1[0];
	t1side2 = vertexes1[2] - vertexes1[1];
	t1side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t1normal = CrossProd(t1side1, t1side2);
	mod_t1normal = Mod(t1normal);
	t1z = t1normal / mod_t1normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t1side1 = Mod(t1side1);
	t1x = t1side1 / mod_t1side1;
	// note that, since x and z are unit vectors, y is too
	t1y = CrossProd(t1x, t1z);

	//
	// calculate local coordinate frame for triangle 2
	//

	// compute side vectors of triangle
	t2side1 = vertexes2[1] - vertexes2[0];
	t2side2 = vertexes2[2] - vertexes2[1];
	t2side3 = vertexes2[0] - vertexes2[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t2normal = CrossProd(t2side1, t2side2);
	mod_t2normal = Mod(t2normal);
	t2z = t2normal / mod_t2normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t2side1 = Mod(t2side1);
	t2x = t2side1 / mod_t2side1;
	// note that, since x and z are unit vectors, y is too
	t2y = CrossProd(t2x, t2z);


	// calculates the coordinates of the triangles in local 2D frame;
	// the origin is located on first vertex.
	// Remark: it is assumed that C2DVertex coordinates defualt to zero

	t1local2D[1].x = mod_t1side1;
	t1local2D[2].x = -DotProd(t1side3, t1x);
	t1local2D[2].y = -DotProd(t1side3, t1y);

	t2local2D[1].x = mod_t2side1;
	t2local2D[2].x = -DotProd(t2side3, t2x);
	t2local2D[2].y = -DotProd(t2side3, t2y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = local2D * unittranmtx;
	//
	// where unittranmtx is the inverse of [0,1,0;1,0,0;1,1,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)

	// calculate transform for triangle 1 and 2
	// (only significative values are calculated, so last line of transform mtx
	// is trivial, moreover t1local2D[0] is (0.0, 0.0) and t1local2D[1].y is 0.0)

	t1transform[0].x = t1local2D[1].x - t1local2D[2].x;
	t1transform[0].y = -t1local2D[2].x;
	t1transform[0].z = t1local2D[2].x;
	t1transform[1].x = -t1local2D[2].y;
	t1transform[1].y = -t1local2D[2].y;
	t1transform[1].z = t1local2D[2].y;

	t2transform[0].x = t2local2D[1].x - t2local2D[2].x;
	t2transform[0].y = -t2local2D[2].x;
	t2transform[0].z = t2local2D[2].x;
	t2transform[1].x = -t2local2D[2].y;
	t2transform[1].y = -t2local2D[2].y;
	t2transform[1].z = t2local2D[2].y;


	// this is the jacobian of the transforms (a determinant)

	t1jacobian = t1transform[0].x*t1transform[1].y-
				t1transform[1].x*t1transform[0].y;

	t2jacobian = t2transform[0].x*t2transform[1].y-
				t2transform[1].x*t2transform[0].y;


	// calculate the transformed 3D points for triangle 1 and 2
	for(i=0; i<POTENTIAL_PTS_2TH_ORDER; i++) {

		// triangle 1

		// calculates the transformed 2D numerical integration point
		point_x = t1transform[0].x * points[i][0] +
			 	  t1transform[0].y * points[i][1] +
			 	  t1transform[0].z * points[i][2];
		point_y = t1transform[1].x * points[i][0] +
				  t1transform[1].y * points[i][1] +
				  t1transform[1].z * points[i][2];

		// calculate 3D point location
		t1p3d[i] = vertexes1[0] + t1x*point_x + t1y*point_y;

		// triangle 2

		// calculates the transformed 2D numerical integration point
		point_x = t2transform[0].x * points[i][0] +
			 	  t2transform[0].y * points[i][1] +
			 	  t2transform[0].z * points[i][2];
		point_y = t2transform[1].x * points[i][0] +
				  t2transform[1].y * points[i][1] +
				  t2transform[1].z * points[i][2];

		// calculate 3D point location
		t2p3d[i] = vertexes2[0] + t2x*point_x + t2y*point_y;
	}

	// calculate the integral
	result = 0;
	for(i=0, k=0; i<POTENTIAL_PTS_2TH_ORDER; i++) {
		for(j=0; j<POTENTIAL_PTS_2TH_ORDER; j++) {
			// calculate the kernel value between given integration points
			pot = 1/Mod(t1p3d[i] - t2p3d[j]);
			// sum up the result
			result = result + weights2[k]*pot;
			// increment weight counter
			k++;
		}
	}

	// do not forget the jacobian!
	// note that we should multiply by the unit triangle area, i.e. 1/2
	// (could be included in the weights; it is a side effect of how, in Stroud,
	// the transformed coordinates are calculated); however we do NOT perform
	// this operation here, since we should then divide by the two triangles area,
	// which include a divide by 2 operation, so they compensate
	result = result * t1jacobian * t2jacobian;

	// The modulus of the triangle normal is 2 * area of triangle
	if( divideByArea == true )
		result = result / (mod_t1normal*mod_t2normal);

	// return result
	return( result / FOUR_PI_TIMES_E0 );
}

// Electric field component in direction normal to one triangular patch
// due to a second triangular patch with uniform charge
//
// This version is fully numerical and uses Strang and Fix, formula #1,
// 3 points, 2th order, for numerical quadrature (from Stroud)
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
// normal is the normal direction vector on the first triangular panel
//
double CPotential::MutualD_2thOrd_FullNum(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], C3DVector normal, bool divideByArea)
{
	C3DVector t1side1, t1side2, t1side3, t1normal, t1x, t1y, t1z;
	C3DVector t2side1, t2side2, t2side3, t2normal, t2x, t2y, t2z;
	C2DVector t1local2D[3], t2local2D[3];
	C3DVector t1p3d[POTENTIAL_PTS_2TH_ORDER], t2p3d[POTENTIAL_PTS_2TH_ORDER];
	C3DVector t1transform[2], t2transform[2];
	C3DVector dist;
	double mod_t1side1, mod_t2side1, mod_t1normal, mod_t2normal;
	double t1jacobian, t2jacobian, pot, result, point_x, point_y;
	double rdist;
	int i, j, k;
	bool isfullfield;

	// if required, calculate full electrical field and not only
	// the normal electrical fields
	if(normal.x == 1.0 && normal.y == 1.0 && normal.z == 1.0)
		isfullfield = true;
	else
		isfullfield = false;

	// numerical 2D quadrature over a triangle
	// uses Stroud 3 points, 2th order
	//
	// initialize matrices of integration points and weights
	// (in homogeneous coordinates)
	static const double points[POTENTIAL_PTS_2TH_ORDER][3] = {
		{0.66666666666666666666666666666667, 0.16666666666666666666666666666667, 1.0},
		{0.16666666666666666666666666666667, 0.66666666666666666666666666666667, 1.0},
		{0.16666666666666666666666666666667, 0.16666666666666666666666666666667, 1.0}};

	static const double weights2[POTENTIAL_W2_2TH_ORDER] =
		{ 0.11111111111111111111111111111111, 0.11111111111111111111111111111111, 0.11111111111111111111111111111111,
		  0.11111111111111111111111111111111, 0.11111111111111111111111111111111, 0.11111111111111111111111111111111,
		  0.11111111111111111111111111111111, 0.11111111111111111111111111111111, 0.11111111111111111111111111111111};

	//
	// calculate local coordinate frame for triangle 1
	//

	// compute side vectors of triangle
	t1side1 = vertexes1[1] - vertexes1[0];
	t1side2 = vertexes1[2] - vertexes1[1];
	t1side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t1normal = CrossProd(t1side1, t1side2);
	mod_t1normal = Mod(t1normal);
	t1z = t1normal / mod_t1normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t1side1 = Mod(t1side1);
	t1x = t1side1 / mod_t1side1;
	// note that, since x and z are unit vectors, y is too
	t1y = CrossProd(t1x, t1z);

	//
	// calculate local coordinate frame for triangle 2
	//

	// compute side vectors of triangle
	t2side1 = vertexes2[1] - vertexes2[0];
	t2side2 = vertexes2[2] - vertexes2[1];
	t2side3 = vertexes2[0] - vertexes2[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t2normal = CrossProd(t2side1, t2side2);
	mod_t2normal = Mod(t2normal);
	t2z = t2normal / mod_t2normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t2side1 = Mod(t2side1);
	t2x = t2side1 / mod_t2side1;
	// note that, since x and z are unit vectors, y is too
	t2y = CrossProd(t2x, t2z);


	// calculates the coordinates of the triangles in local 2D frame;
	// the origin is located on first vertex.
	// Remark: it is assumed that C2DVertex coordinates defualt to zero

	t1local2D[1].x = mod_t1side1;
	t1local2D[2].x = -DotProd(t1side3, t1x);
	t1local2D[2].y = -DotProd(t1side3, t1y);

	t2local2D[1].x = mod_t2side1;
	t2local2D[2].x = -DotProd(t2side3, t2x);
	t2local2D[2].y = -DotProd(t2side3, t2y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = local2D * unittranmtx;
	//
	// where unittranmtx is the inverse of [0,1,0;1,0,0;1,1,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)

	// calculate transform for triangle 1 and 2
	// (only significative values are calculated, so last line of transform mtx
	// is trivial, moreover t1local2D[0] is (0.0, 0.0) and t1local2D[1].y is 0.0)

	t1transform[0].x = t1local2D[1].x - t1local2D[2].x;
	t1transform[0].y = -t1local2D[2].x;
	t1transform[0].z = t1local2D[2].x;
	t1transform[1].x = -t1local2D[2].y;
	t1transform[1].y = -t1local2D[2].y;
	t1transform[1].z = t1local2D[2].y;

	t2transform[0].x = t2local2D[1].x - t2local2D[2].x;
	t2transform[0].y = -t2local2D[2].x;
	t2transform[0].z = t2local2D[2].x;
	t2transform[1].x = -t2local2D[2].y;
	t2transform[1].y = -t2local2D[2].y;
	t2transform[1].z = t2local2D[2].y;


	// this is the jacobian of the transforms (a determinant)

	t1jacobian = t1transform[0].x*t1transform[1].y-
				t1transform[1].x*t1transform[0].y;

	t2jacobian = t2transform[0].x*t2transform[1].y-
				t2transform[1].x*t2transform[0].y;


	// calculate the transformed 3D points for triangle 1 and 2
	for(i=0; i<POTENTIAL_PTS_2TH_ORDER; i++) {

		// triangle 1

		// calculates the transformed 2D numerical integration point
		point_x = t1transform[0].x * points[i][0] +
			 	  t1transform[0].y * points[i][1] +
			 	  t1transform[0].z * points[i][2];
		point_y = t1transform[1].x * points[i][0] +
				  t1transform[1].y * points[i][1] +
				  t1transform[1].z * points[i][2];

		// calculate 3D point location
		t1p3d[i] = vertexes1[0] + t1x*point_x + t1y*point_y;

		// triangle 2

		// calculates the transformed 2D numerical integration point
		point_x = t2transform[0].x * points[i][0] +
			 	  t2transform[0].y * points[i][1] +
			 	  t2transform[0].z * points[i][2];
		point_y = t2transform[1].x * points[i][0] +
				  t2transform[1].y * points[i][1] +
				  t2transform[1].z * points[i][2];

		// calculate 3D point location
		t2p3d[i] = vertexes2[0] + t2x*point_x + t2y*point_y;
	}

	// calculate the integral
	result = 0;
	for(i=0, k=0; i<POTENTIAL_PTS_2TH_ORDER; i++) {
		for(j=0; j<POTENTIAL_PTS_2TH_ORDER; j++) {
			// calculate the kernel value between given integration points:
			// calculate the electric field component in direction
			// 'normal' at given point due to the triangle patch
			dist = t1p3d[i] - t2p3d[j];
			rdist = Mod(dist);
			if( isfullfield == false ) {
				// this is cos((x1-x2),n1) * mod(x1-x2)
				pot = DotProd(dist, normal) / (rdist * rdist * rdist);
			}
			else {
				pot = 1.0 / (rdist * rdist);
			}

			// sum up the result
			result = result + weights2[k]*pot;
			// increment weight counter
			k++;
		}
	}

	// do not forget the jacobian!
	// note that we should multiply by the unit triangle area, i.e. 1/2
	// (could be included in the weights; it is a side effect of how, in Stroud,
	// the transformed coordinates are calculated); however we do NOT perform
	// this operation here, since we should then divide by the two triangles area,
	// which include a divide by 2 operation, so they compensate
	result = result * t1jacobian * t2jacobian;

	// The modulus of the triangle normal is 2 * area of triangle
	if( divideByArea == true )
		result = result / (mod_t1normal*mod_t2normal);

	// return result
	return( result / FOUR_PI_TIMES_E0 );
}

// Mutual cofficient of potential between two triangular patches
// with uniform charge
//
// This version is fully numerical and uses Strang and Fix, formula #3,
// 4 points, 3th order, for numerical quadrature (from Stroud)
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
double CPotential::Mutual_3thOrd_FullNum(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], bool divideByArea)
{
	C3DVector t1side1, t1side2, t1side3, t1normal, t1x, t1y, t1z;
	C3DVector t2side1, t2side2, t2side3, t2normal, t2x, t2y, t2z;
	C2DVector t1local2D[3], t2local2D[3];
	C3DVector t1p3d[POTENTIAL_PTS_3TH_ORDER], t2p3d[POTENTIAL_PTS_3TH_ORDER];
	C3DVector t1transform[2], t2transform[2];
	double mod_t1side1, mod_t2side1, mod_t1normal, mod_t2normal;
	double t1jacobian, t2jacobian, pot, result, point_x, point_y;
	int i, j, k;

	// numerical 2D quadrature over a triangle
	// uses Stroud 4 points, 3th order
	//
	// initialize matrices of integration points and weights
	// (in homogeneous coordinates)
	static const double points[POTENTIAL_PTS_3TH_ORDER][3] = {
		{0.33333333333333333333333333333333, 0.33333333333333333333333333333333, 1.0},
		{0.6,                                0.2,                                1.0},
		{0.2,                                0.6,                                1.0},
		{0.2,                                0.2,                                1.0}};

	static const double weights2[POTENTIAL_W2_3TH_ORDER] =
		{ 0.31640625, -0.29296875,                         -0.29296875,                         -0.29296875,
		 -0.29296875,  0.27126736111111111111111111111111,  0.27126736111111111111111111111111,  0.27126736111111111111111111111111,
		 -0.29296875,  0.27126736111111111111111111111111,  0.27126736111111111111111111111111,  0.27126736111111111111111111111111,
		 -0.29296875,  0.27126736111111111111111111111111,  0.27126736111111111111111111111111,  0.27126736111111111111111111111111};

	//
	// calculate local coordinate frame for triangle 1
	//

	// compute side vectors of triangle
	t1side1 = vertexes1[1] - vertexes1[0];
	t1side2 = vertexes1[2] - vertexes1[1];
	t1side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t1normal = CrossProd(t1side1, t1side2);
	mod_t1normal = Mod(t1normal);
	t1z = t1normal / mod_t1normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t1side1 = Mod(t1side1);
	t1x = t1side1 / mod_t1side1;
	// note that, since x and z are unit vectors, y is too
	t1y = CrossProd(t1x, t1z);

	//
	// calculate local coordinate frame for triangle 2
	//

	// compute side vectors of triangle
	t2side1 = vertexes2[1] - vertexes2[0];
	t2side2 = vertexes2[2] - vertexes2[1];
	t2side3 = vertexes2[0] - vertexes2[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t2normal = CrossProd(t2side1, t2side2);
	mod_t2normal = Mod(t2normal);
	t2z = t2normal / mod_t2normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t2side1 = Mod(t2side1);
	t2x = t2side1 / mod_t2side1;
	// note that, since x and z are unit vectors, y is too
	t2y = CrossProd(t2x, t2z);


	// calculates the coordinates of the triangles in local 2D frame;
	// the origin is located on first vertex.
	// Remark: it is assumed that C2DVertex coordinates defualt to zero

	t1local2D[1].x = mod_t1side1;
	t1local2D[2].x = -DotProd(t1side3, t1x);
	t1local2D[2].y = -DotProd(t1side3, t1y);

	t2local2D[1].x = mod_t2side1;
	t2local2D[2].x = -DotProd(t2side3, t2x);
	t2local2D[2].y = -DotProd(t2side3, t2y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = local2D * unittranmtx;
	//
	// where unittranmtx is the inverse of [0,1,0;1,0,0;1,1,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)

	// calculate transform for triangle 1 and 2
	// (only significative values are calculated, so last line of transform mtx
	// is trivial, moreover t1local2D[0] is (0.0, 0.0) and t1local2D[1].y is 0.0)

	t1transform[0].x = t1local2D[1].x - t1local2D[2].x;
	t1transform[0].y = -t1local2D[2].x;
	t1transform[0].z = t1local2D[2].x;
	t1transform[1].x = -t1local2D[2].y;
	t1transform[1].y = -t1local2D[2].y;
	t1transform[1].z = t1local2D[2].y;

	t2transform[0].x = t2local2D[1].x - t2local2D[2].x;
	t2transform[0].y = -t2local2D[2].x;
	t2transform[0].z = t2local2D[2].x;
	t2transform[1].x = -t2local2D[2].y;
	t2transform[1].y = -t2local2D[2].y;
	t2transform[1].z = t2local2D[2].y;


	// this is the jacobian of the transforms (a determinant)

	t1jacobian = t1transform[0].x*t1transform[1].y-
				t1transform[1].x*t1transform[0].y;

	t2jacobian = t2transform[0].x*t2transform[1].y-
				t2transform[1].x*t2transform[0].y;


	// calculate the transformed 3D points for triangle 1 and 2
	for(i=0; i<POTENTIAL_PTS_3TH_ORDER; i++) {

		// triangle 1

		// calculates the transformed 2D numerical integration point
		point_x = t1transform[0].x * points[i][0] +
			 	  t1transform[0].y * points[i][1] +
			 	  t1transform[0].z * points[i][2];
		point_y = t1transform[1].x * points[i][0] +
				  t1transform[1].y * points[i][1] +
				  t1transform[1].z * points[i][2];

		// calculate 3D point location
		t1p3d[i] = vertexes1[0] + t1x*point_x + t1y*point_y;

		// triangle 2

		// calculates the transformed 2D numerical integration point
		point_x = t2transform[0].x * points[i][0] +
			 	  t2transform[0].y * points[i][1] +
			 	  t2transform[0].z * points[i][2];
		point_y = t2transform[1].x * points[i][0] +
				  t2transform[1].y * points[i][1] +
				  t2transform[1].z * points[i][2];

		// calculate 3D point location
		t2p3d[i] = vertexes2[0] + t2x*point_x + t2y*point_y;
	}

	// calculate the integral
	result = 0;
	for(i=0, k=0; i<POTENTIAL_PTS_3TH_ORDER; i++) {
		for(j=0; j<POTENTIAL_PTS_3TH_ORDER; j++) {
			// calculate the kernel value between given integration points
			pot = 1/Mod(t1p3d[i] - t2p3d[j]);
			// sum up the result
			result = result + weights2[k]*pot;
			// increment weight counter
			k++;
		}
	}

	// do not forget the jacobian!
	// note that we should multiply by the unit triangle area, i.e. 1/2
	// (could be included in the weights; it is a side effect of how, in Stroud,
	// the transformed coordinates are calculated); however we do NOT perform
	// this operation here, since we should then divide by the two triangles area,
	// which include a divide by 2 operation, so they compensate
	result = result * t1jacobian * t2jacobian;

	// The modulus of the triangle normal is 2 * area of triangle
	if( divideByArea == true )
		result /= mod_t1normal*mod_t2normal;

	// return result
	return( result / FOUR_PI_TIMES_E0 );
}

// Mutual cofficient of potential between two triangular patches
// with uniform charge
//
// This version is fully numerical and uses Strang and Fix, formula #5,
// 6 points, 4th order, for numerical quadrature (from Stroud)
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
double CPotential::Mutual_4thOrd_FullNum(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], bool divideByArea)
{
	C3DVector t1side1, t1side2, t1side3, t1normal, t1x, t1y, t1z;
	C3DVector t2side1, t2side2, t2side3, t2normal, t2x, t2y, t2z;
	C2DVector t1local2D[3], t2local2D[3];
	C3DVector t1p3d[POTENTIAL_PTS_4TH_ORDER], t2p3d[POTENTIAL_PTS_4TH_ORDER];
	C3DVector t1transform[2], t2transform[2];
	double mod_t1side1, mod_t2side1, mod_t1normal, mod_t2normal;
	double t1jacobian, t2jacobian, pot, result, point_x, point_y;
	int i, j, k;

	// numerical 2D quadrature over a triangle
	// uses Stroud 6 points, 4th order
	//
	// initialize matrices of integration points and weights
	// (in homogeneous coordinates)
	static const double points[POTENTIAL_PTS_4TH_ORDER][3] = {
		{0.816847572980459, 0.091576213509771, 1.0},
		{0.091576213509771, 0.816847572980459, 1.0},
		{0.091576213509771, 0.091576213509771, 1.0},
		{0.108103018168070, 0.445948490915965, 1.0},
		{0.445948490915965, 0.108103018168070, 1.0},
		{0.445948490915965, 0.445948490915965, 1.0}};

	static const double weights2[POTENTIAL_W2_4TH_ORDER] =
		{ 0.012089385932845641681938923684, 0.012089385932845641681938923684, 0.012089385932845641681938923684, 0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.024561195285594988334146524542,
	      0.012089385932845641681938923684, 0.012089385932845641681938923684, 0.012089385932845641681938923684, 0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.024561195285594988334146524542,
	      0.012089385932845641681938923684, 0.012089385932845641681938923684, 0.012089385932845641681938923684, 0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.024561195285594988334146524542,
          0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.049899334607075270538656916121, 0.049899334607075270538656916121, 0.049899334607075270538656916121,
          0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.049899334607075270538656916121, 0.049899334607075270538656916121, 0.049899334607075270538656916121,
          0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.024561195285594988334146524542, 0.049899334607075270538656916121, 0.049899334607075270538656916121, 0.049899334607075270538656916121};

	//
	// calculate local coordinate frame for triangle 1
	//

	// compute side vectors of triangle
	t1side1 = vertexes1[1] - vertexes1[0];
	t1side2 = vertexes1[2] - vertexes1[1];
	t1side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t1normal = CrossProd(t1side1, t1side2);
	mod_t1normal = Mod(t1normal);
	t1z = t1normal / mod_t1normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t1side1 = Mod(t1side1);
	t1x = t1side1 / mod_t1side1;
	// note that, since x and z are unit vectors, y is too
	t1y = CrossProd(t1x, t1z);

	//
	// calculate local coordinate frame for triangle 2
	//

	// compute side vectors of triangle
	t2side1 = vertexes2[1] - vertexes2[0];
	t2side2 = vertexes2[2] - vertexes2[1];
	t2side3 = vertexes2[0] - vertexes2[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t2normal = CrossProd(t2side1, t2side2);
	mod_t2normal = Mod(t2normal);
	t2z = t2normal / mod_t2normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t2side1 = Mod(t2side1);
	t2x = t2side1 / mod_t2side1;
	// note that, since x and z are unit vectors, y is too
	t2y = CrossProd(t2x, t2z);


	// calculates the coordinates of the triangles in local 2D frame;
	// the origin is located on first vertex.
	// Remark: it is assumed that C2DVertex coordinates defualt to zero

	t1local2D[1].x = mod_t1side1;
	t1local2D[2].x = -DotProd(t1side3, t1x);
	t1local2D[2].y = -DotProd(t1side3, t1y);

	t2local2D[1].x = mod_t2side1;
	t2local2D[2].x = -DotProd(t2side3, t2x);
	t2local2D[2].y = -DotProd(t2side3, t2y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = local2D * unittranmtx;
	//
	// where unittranmtx is the inverse of [0,1,0;1,0,0;1,1,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)

	// calculate transform for triangle 1 and 2
	// (only significative values are calculated, so last line of transform mtx
	// is trivial, moreover t1local2D[0] is (0.0, 0.0) and t1local2D[1].y is 0.0)

	t1transform[0].x = t1local2D[1].x - t1local2D[2].x;
	t1transform[0].y = -t1local2D[2].x;
	t1transform[0].z = t1local2D[2].x;
	t1transform[1].x = -t1local2D[2].y;
	t1transform[1].y = -t1local2D[2].y;
	t1transform[1].z = t1local2D[2].y;

	t2transform[0].x = t2local2D[1].x - t2local2D[2].x;
	t2transform[0].y = -t2local2D[2].x;
	t2transform[0].z = t2local2D[2].x;
	t2transform[1].x = -t2local2D[2].y;
	t2transform[1].y = -t2local2D[2].y;
	t2transform[1].z = t2local2D[2].y;


	// this is the jacobian of the transforms (a determinant)

	t1jacobian = t1transform[0].x*t1transform[1].y-
				t1transform[1].x*t1transform[0].y;

	t2jacobian = t2transform[0].x*t2transform[1].y-
				t2transform[1].x*t2transform[0].y;


	// calculate the transformed 3D points for triangle 1 and 2
	for(i=0; i<POTENTIAL_PTS_4TH_ORDER; i++) {

		// triangle 1

		// calculates the transformed 2D numerical integration point
		point_x = t1transform[0].x * points[i][0] +
			 	  t1transform[0].y * points[i][1] +
			 	  t1transform[0].z * points[i][2];
		point_y = t1transform[1].x * points[i][0] +
				  t1transform[1].y * points[i][1] +
				  t1transform[1].z * points[i][2];

		// calculate 3D point location
		t1p3d[i] = vertexes1[0] + t1x*point_x + t1y*point_y;

		// triangle 2

		// calculates the transformed 2D numerical integration point
		point_x = t2transform[0].x * points[i][0] +
			 	  t2transform[0].y * points[i][1] +
			 	  t2transform[0].z * points[i][2];
		point_y = t2transform[1].x * points[i][0] +
				  t2transform[1].y * points[i][1] +
				  t2transform[1].z * points[i][2];

		// calculate 3D point location
		t2p3d[i] = vertexes2[0] + t2x*point_x + t2y*point_y;
	}

	// calculate the integral
	result = 0;
	for(i=0, k=0; i<POTENTIAL_PTS_4TH_ORDER; i++) {
		for(j=0; j<POTENTIAL_PTS_4TH_ORDER; j++) {
			// calculate the kernel value between given integration points
			pot = 1/Mod(t1p3d[i] - t2p3d[j]);
			// sum up the result
			result = result + weights2[k]*pot;
			// increment weight counter
			k++;
		}
	}

	// do not forget the jacobian!
	// note that we should multiply by the unit triangle area, i.e. 1/2
	// (could be included in the weights; it is a side effect of how, in Stroud,
	// the transformed coordinates are calculated); however we do NOT perform
	// this operation here, since we should then divide by the two triangles area,
	// which include a divide by 2 operation, so they compensate
	result = result * t1jacobian * t2jacobian;

	// The modulus of the triangle normal is 2 * area of triangle
	if( divideByArea == true )
		result = result / (mod_t1normal*mod_t2normal);

	// return result
	return( result / FOUR_PI_TIMES_E0 );
}

// Mutual cofficient of potential between two triangular patches
// with uniform charge
//
// This version is fully numerical
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
// rule is the integration formula to use, see InitNumerical()
double CPotential::MutualFullNumerical(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], int rule, bool divideByArea)
{
	C3DVector t1side1, t1side2, t1side3, t1normal, t1x, t1y, t1z;
	C3DVector t2side1, t2side2, t2side3, t2normal, t2x, t2y, t2z;
	C2DVector t1local2D[3], t2local2D[3];
	C3DVector t1p3d[POTENTIAL_PTS_MAX], t2p3d[POTENTIAL_PTS_MAX];
	C3DVector t1transform[2], t2transform[2], dist;
	double mod_t1side1, mod_t2side1, mod_t1normal, mod_t2normal;
	double t1jacobian, t2jacobian, pot, result, point_x, point_y;
	unsigned int i, j, k;

	//
	// calculate local coordinate frame for triangle 1
	//

	// compute side vectors of triangle
	t1side1 = vertexes1[1] - vertexes1[0];
	t1side2 = vertexes1[2] - vertexes1[1];
	t1side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t1normal = CrossProd(t1side1, t1side2);
	mod_t1normal = Mod(t1normal);
	t1z = t1normal / mod_t1normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t1side1 = Mod(t1side1);
	t1x = t1side1 / mod_t1side1;
	// note that, since x and z are unit vectors, y is too
	t1y = CrossProd(t1x, t1z);

	//
	// calculate local coordinate frame for triangle 2
	//

	// compute side vectors of triangle
	t2side1 = vertexes2[1] - vertexes2[0];
	t2side2 = vertexes2[2] - vertexes2[1];
	t2side3 = vertexes2[0] - vertexes2[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t2normal = CrossProd(t2side1, t2side2);
	mod_t2normal = Mod(t2normal);
	t2z = t2normal / mod_t2normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t2side1 = Mod(t2side1);
	t2x = t2side1 / mod_t2side1;
	// note that, since x and z are unit vectors, y is too
	t2y = CrossProd(t2x, t2z);


	// calculates the coordinates of the triangles in local 2D frame;
	// the origin is located on first vertex.
	// Remark: it is assumed that C2DVertex coordinates defualt to zero

	t1local2D[1].x = mod_t1side1;
	t1local2D[2].x = -DotProd(t1side3, t1x);
	t1local2D[2].y = -DotProd(t1side3, t1y);

	t2local2D[1].x = mod_t2side1;
	t2local2D[2].x = -DotProd(t2side3, t2x);
	t2local2D[2].y = -DotProd(t2side3, t2y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = local2D * unittranmtx;
	//
	// where unittranmtx is the inverse of [0,1,0;1,0,0;1,1,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)

	// calculate transform for triangle 1 and 2
	// (only significative values are calculated, so last line of transform mtx
	// is trivial, moreover t1local2D[0] is (0.0, 0.0) and t1local2D[1].y is 0.0)

	t1transform[0].x = t1local2D[1].x - t1local2D[2].x;
	t1transform[0].y = -t1local2D[2].x;
	t1transform[0].z = t1local2D[2].x;
	t1transform[1].x = -t1local2D[2].y;
	t1transform[1].y = -t1local2D[2].y;
	t1transform[1].z = t1local2D[2].y;

	t2transform[0].x = t2local2D[1].x - t2local2D[2].x;
	t2transform[0].y = -t2local2D[2].x;
	t2transform[0].z = t2local2D[2].x;
	t2transform[1].x = -t2local2D[2].y;
	t2transform[1].y = -t2local2D[2].y;
	t2transform[1].z = t2local2D[2].y;


	// this is the jacobian of the transforms (a determinant)

	t1jacobian = t1transform[0].x*t1transform[1].y-
				t1transform[1].x*t1transform[0].y;

	t2jacobian = t2transform[0].x*t2transform[1].y-
				t2transform[1].x*t2transform[0].y;


	// calculate the transformed 3D points for triangle 1 and 2
	for(i=0; i<m_uiNorder[rule]; i++) {

		// triangle 1

		// calculates the transformed 2D numerical integration point
		point_x = t1transform[0].x * m_clsPoints[rule][i][0] +
			 	  t1transform[0].y * m_clsPoints[rule][i][1] +
			 	  t1transform[0].z * m_clsPoints[rule][i][2];
		point_y = t1transform[1].x * m_clsPoints[rule][i][0] +
				  t1transform[1].y * m_clsPoints[rule][i][1] +
				  t1transform[1].z * m_clsPoints[rule][i][2];

		// calculate 3D point location
		t1p3d[i] = vertexes1[0] + t1x*point_x + t1y*point_y;

		// triangle 2

		// calculates the transformed 2D numerical integration point
		point_x = t2transform[0].x * m_clsPoints[rule][i][0] +
			 	  t2transform[0].y * m_clsPoints[rule][i][1] +
			 	  t2transform[0].z * m_clsPoints[rule][i][2];
		point_y = t2transform[1].x * m_clsPoints[rule][i][0] +
				  t2transform[1].y * m_clsPoints[rule][i][1] +
				  t2transform[1].z * m_clsPoints[rule][i][2];

		// calculate 3D point location
		t2p3d[i] = vertexes2[0] + t2x*point_x + t2y*point_y;
	}

	// calculate the integral
	result = 0;
	for(i=0, k=0; i<m_uiNorder[rule]; i++) {
		for(j=0; j<m_uiNorder[rule]; j++) {
			// calculate the kernel value between given integration points
			pot = 1/Mod(t1p3d[i] - t2p3d[j]);
			// sum up the result
			result = result + m_dWeight[rule][i]*m_dWeight[rule][j]*pot;
			// increment weight counter
			k++;
		}
	}

	// do not forget the jacobian!
	result = result * t1jacobian * t2jacobian;

	// The modulus of the triangle normal is 2 * area of triangle
	if( divideByArea == true )
		result = result / (mod_t1normal*mod_t2normal);

	// return result
	return( result / FOUR_PI_TIMES_E0 );
}

// Electric field component in direction normal to one triangular patch
// due to a second triangular patch with uniform charge
//
// This version is fully numerical
//
// vertexes1 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the first triangle
//
// vertexes2 is an array of 3 3D point coordinates,
//           corresponding to the vertexes of the second triangle
//
// normal is the normal direction vector on the first triangular panel
//
// rule is the integration formula to use, see InitNumerical()
//
double CPotential::MutualDFullNumerical(C3DVector_float vertexes1[3],
						  C3DVector_float vertexes2[3], C3DVector normal, int rule, bool divideByArea)
{
	C3DVector t1side1, t1side2, t1side3, t1normal, t1x, t1y, t1z;
	C3DVector t2side1, t2side2, t2side3, t2normal, t2x, t2y, t2z;
	C2DVector t1local2D[3], t2local2D[3];
	C3DVector t1p3d[POTENTIAL_PTS_MAX], t2p3d[POTENTIAL_PTS_MAX];
	C3DVector t1transform[2], t2transform[2], dist;
	double mod_t1side1, mod_t2side1, mod_t1normal, mod_t2normal;
	double t1jacobian, t2jacobian, pot, result, point_x, point_y, rdist;
	unsigned int i, j, k;

	//
	// calculate local coordinate frame for triangle 1
	//

	// compute side vectors of triangle
	t1side1 = vertexes1[1] - vertexes1[0];
	t1side2 = vertexes1[2] - vertexes1[1];
	t1side3 = vertexes1[0] - vertexes1[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t1normal = CrossProd(t1side1, t1side2);
	mod_t1normal = Mod(t1normal);
	t1z = t1normal / mod_t1normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t1side1 = Mod(t1side1);
	t1x = t1side1 / mod_t1side1;
	// note that, since x and z are unit vectors, y is too
	t1y = CrossProd(t1x, t1z);

	//
	// calculate local coordinate frame for triangle 2
	//

	// compute side vectors of triangle
	t2side1 = vertexes2[1] - vertexes2[0];
	t2side2 = vertexes2[2] - vertexes2[1];
	t2side3 = vertexes2[0] - vertexes2[2];

	// compute unit vector normal to triangle plane (uses vector product)
	t2normal = CrossProd(t2side1, t2side2);
	mod_t2normal = Mod(t2normal);
	t2z = t2normal / mod_t2normal;

	// compute local x, y axes versors; x is along side1,
	// y is perp. to side1 in counter-clockwise direction
	mod_t2side1 = Mod(t2side1);
	t2x = t2side1 / mod_t2side1;
	// note that, since x and z are unit vectors, y is too
	t2y = CrossProd(t2x, t2z);


	// calculates the coordinates of the triangles in local 2D frame;
	// the origin is located on first vertex.
	// Remark: it is assumed that C2DVertex coordinates defualt to zero

	t1local2D[1].x = mod_t1side1;
	t1local2D[2].x = -DotProd(t1side3, t1x);
	t1local2D[2].y = -DotProd(t1side3, t1y);

	t2local2D[1].x = mod_t2side1;
	t2local2D[2].x = -DotProd(t2side3, t2x);
	t2local2D[2].y = -DotProd(t2side3, t2y);

	// calculates the coordinate transformation needed
	// to map the first triangle to a unit triangle
	// Basically, this is the following operation
	//
	//	static const double unittranmtx[9] = { 0, 1, 0,
	//	                                      1, 0, 0,
	//							             -1,-1, 1};
	// 	transform = local2D * unittranmtx;
	//
	// where unittranmtx is the inverse of [0,1,0;1,0,0;1,1,1],
	// the matrix representing the unit triangle
	//
	// Explanation: we are in 2D, but we use homogeneous coordinates.
	// The mapping is a linear transformation, i.e. to map a point (xt, yt, 1)
	// from a unit triangle to a point (x0, y0, 1) on a general triangle we use:
	// x0 = a*xt + b*yt + c
	// y0 = d*xt + e*yt + f
	// Therefore to map the corners of a unit triangle to the corners of a general
	// triangle we have:
	// [ x0 x1 x2 ]    [ a b c ]   [ 0 1 0 ]
	// [ y0 y1 y2 ]  = [ d e f ] * [ 1 0 0 ]
	// [  1  1  1 ]    [ 0 0 1 ]   [ 1 1 1 ]
	// where the last matrix is the matrix representing the vertex of the unit
	// triangle in homogeneous coordinates
	// We would like to know the matrix Trans = [a b c; d e f; 0 0 1] to be able to map
	// integration points in the unit triangle to integration points in our triangle.
	// The matrix is therefore Trans = MyTri * inv(UnitTri)

	// calculate transform for triangle 1 and 2
	// (only significative values are calculated, so last line of transform mtx
	// is trivial, moreover t1local2D[0] is (0.0, 0.0) and t1local2D[1].y is 0.0)

	t1transform[0].x = t1local2D[1].x - t1local2D[2].x;
	t1transform[0].y = -t1local2D[2].x;
	t1transform[0].z = t1local2D[2].x;
	t1transform[1].x = -t1local2D[2].y;
	t1transform[1].y = -t1local2D[2].y;
	t1transform[1].z = t1local2D[2].y;

	t2transform[0].x = t2local2D[1].x - t2local2D[2].x;
	t2transform[0].y = -t2local2D[2].x;
	t2transform[0].z = t2local2D[2].x;
	t2transform[1].x = -t2local2D[2].y;
	t2transform[1].y = -t2local2D[2].y;
	t2transform[1].z = t2local2D[2].y;


	// this is the jacobian of the transforms (a determinant)

	t1jacobian = t1transform[0].x*t1transform[1].y-
				t1transform[1].x*t1transform[0].y;

	t2jacobian = t2transform[0].x*t2transform[1].y-
				t2transform[1].x*t2transform[0].y;


	// calculate the transformed 3D points for triangle 1 and 2
	for(i=0; i<m_uiNorder[rule]; i++) {

		// triangle 1

		// calculates the transformed 2D numerical integration point
		point_x = t1transform[0].x * m_clsPoints[rule][i][0] +
			 	  t1transform[0].y * m_clsPoints[rule][i][1] +
			 	  t1transform[0].z * m_clsPoints[rule][i][2];
		point_y = t1transform[1].x * m_clsPoints[rule][i][0] +
				  t1transform[1].y * m_clsPoints[rule][i][1] +
				  t1transform[1].z * m_clsPoints[rule][i][2];

		// calculate 3D point location
		t1p3d[i] = vertexes1[0] + t1x*point_x + t1y*point_y;

		// triangle 2

		// calculates the transformed 2D numerical integration point
		point_x = t2transform[0].x * m_clsPoints[rule][i][0] +
			 	  t2transform[0].y * m_clsPoints[rule][i][1] +
			 	  t2transform[0].z * m_clsPoints[rule][i][2];
		point_y = t2transform[1].x * m_clsPoints[rule][i][0] +
				  t2transform[1].y * m_clsPoints[rule][i][1] +
				  t2transform[1].z * m_clsPoints[rule][i][2];

		// calculate 3D point location
		t2p3d[i] = vertexes2[0] + t2x*point_x + t2y*point_y;
	}

	// calculate the integral
	result = 0;
	for(i=0, k=0; i<m_uiNorder[rule]; i++) {
		for(j=0; j<m_uiNorder[rule]; j++) {
			// calculate the kernel value between given integration points:
			// calculate the electric field component in direction
			// 'normal' at given point due to the triangle patch
			dist = t1p3d[i] - t2p3d[j];
			rdist = Mod(dist);
			// this is cos((x1-x2),n1) * mod(x1-x2)
			pot = DotProd(dist, normal) / (rdist * rdist * rdist);

			// sum up the result
			result = result + m_dWeight[rule][i]*m_dWeight[rule][j]*pot;
			// increment weight counter
			k++;
		}
	}

	// do not forget the jacobian!
	result = result * t1jacobian * t2jacobian;

	// The modulus of the triangle normal is 2 * area of triangle
	if( divideByArea == true )
		result = result / (mod_t1normal*mod_t2normal);

	// return result
	return( result / FOUR_PI_TIMES_E0 );
}
