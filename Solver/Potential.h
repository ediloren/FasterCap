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


// Potential.h : CPotential class header file
//

#if !defined(AFX_POTENTIAL_H__E89AAD21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_POTENTIAL_H__E89AAD21_5486_11D5_9282_04F014C10000__INCLUDED_

#include "SolverGlobal.h"

// includes for LinALg
#include "Geometry/Vector3D.h"
#include "Geometry/Vector2D.h"

// includes for CAutoPanel class
// to pass panels directly as arguments
#include "AutoPanel.h"

#define PI_TIMES_E0			2.781625139939e-11
#define TWO_PI_TIMES_E0	    5.563250279878e-11
#define FOUR_PI_TIMES_E0	1.1126500559756e-10
#define TWO_TIMES_E0		1.7708375634e-11
#define E0                  8.854187817E-12
#define PI					3.1415926535897932384626433832795
#define PI_HALF     		1.5707963267948966192313216916398

// tolerance for distances in potential calculations (e.g. to check if a point lies on a line etc.)
#define POTENTIAL_GEO_EPS 1E-6
// numerical tolerance for singularity check
#define POTENTIAL_TOL_EPS 1E-12

// max number of points for the n-th order quadrature
#define POTENTIAL_PTS_MAX			144
// max number of quadrature rules
#define POTENTIAL_RULE_MAX			20

class CPotential
{
	class CPotParam
	{
	public:
		C3DVector_float *vertex1, *vertex2;
		C3DVector r, n, rho;
		C3DVector u, Rminus, Rplus, P0v;
		double d, P0, lplus, lminus, R0;
	};

public:
	CPotential();
//	double Auto(C3DVector vertexes[3]);
	double Auto(C3DVector_float vertexes[3]);
	double Auto(CAutoPanel &panel);
	double QAuto(C3DVector vertexes[4]);
//	double AutoNumerical(C3DVector vertexes[3], int rule, int formula, bool divideByArea);
	double AutoNumerical(C3DVector_float vertexes[3], int rule, int formula, bool divideByArea);
//	double InsidePotential(C3DVector point, C3DVector vertexes[3]);
	double InsidePotential(C3DVector point, C3DVector_float vertexes[3]);
//	double CornerPotential(C3DVector vertexes[3], bool divideByArea = true);
	double CornerPotential(C3DVector_float vertexes[3], bool divideByArea = true);
//	double Potential(C3DVector r, C3DVector vertexes[3]);
	double Potential(C3DVector r, C3DVector_float vertexes[3]);
//	double PotentialOpt(C3DVector r, C3DVector vertexes[3], bool divideByArea = true);
	double PotentialOpt(C3DVector r, C3DVector_float vertexes[3], bool divideByArea = true);
	double PotentialOpt(C3DVector r, CAutoPanel panel, bool divideByArea = true);
    double PotentialOpt(C2DVector r, C2DVector_float vertexes[2], bool divideByLen = true);
//	double PotentialNumerical(C3DVector r, C3DVector vertexes[3], int rule);
	double PotentialNumerical(C3DVector r, C3DVector_float vertexes[3], int rule);
//	double EnFieldNumerical(C3DVector r, C3DVector vertexes[3], C3DVector normal, int rule);
	double EnFieldNumerical(C3DVector r, C3DVector_float vertexes[3], C3DVector normal, int rule);
    double EnField(C2DVector r, C2DVector_float vertexes[2], C2DVector pnormal, bool divideByLen = true);
//	double Mutual_9thOrd_HalfNum(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], bool divideByArea = true);
	double Mutual_9thOrd_HalfNum(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], bool divideByArea = true);
	double QMutual_9thOrd_HalfNum(C3DVector vertexes1[4], C3DVector vertexes2[4],
					bool divideByArea = true);
//	double MutualHalfNumerical(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], int rule, bool divideByArea = true);
	double MutualHalfNumerical(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], int rule, bool divideByArea = true);
	double MutualHalfNumerical(CAutoPanel panel1,
				  CAutoPanel panel2, int rule, bool divideByArea = true);
	double MutualDHalfNumerical(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3],C3DVector normal, double h, int rule, bool divideByArea = true);
	double MutualDHalfNumerical(CAutoPanel panel1,
						  CAutoPanel panel2, C3DVector normal, double h, int rule, bool divideByArea = true);
//	double Mutual_2thOrd_FullNum(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], bool divideByArea = true);
	double Mutual_2thOrd_FullNum(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], bool divideByArea = true);
//	double MutualD_2thOrd_FullNum(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], C3DVector normal, bool divideByArea = true);
	double MutualD_2thOrd_FullNum(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], C3DVector normal, bool divideByArea = true);
//	double Mutual_3thOrd_FullNum(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], bool divideByArea = true);
	double Mutual_3thOrd_FullNum(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], bool divideByArea = true);
//	double Mutual_4thOrd_FullNum(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], bool divideByArea = true);
	double Mutual_4thOrd_FullNum(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], bool divideByArea = true);
//	double MutualFullNumerical(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], int rule, bool divideByArea = true);
	double MutualFullNumerical(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], int rule, bool divideByArea = true);
//	double MutualDFullNumerical(C3DVector vertexes1[3],
//				  C3DVector vertexes2[3], C3DVector normal, int rule, bool divideByArea = true);
	double MutualDFullNumerical(C3DVector_float vertexes1[3],
				  C3DVector_float vertexes2[3], C3DVector normal, int rule, bool divideByArea = true);

protected:
	void InitNumerical();
	void Sidecontrib(CPotParam *param);

	C3DVector m_clsPoints[POTENTIAL_RULE_MAX][POTENTIAL_PTS_MAX];
	double m_dWeight[POTENTIAL_RULE_MAX][POTENTIAL_PTS_MAX];
	unsigned int m_uiNorder[POTENTIAL_RULE_MAX];
};


#endif //!defined(AFX_POTENTIAL_H__E89AAD21_5486_11D5_9282_04F014C10000__INCLUDED_)
