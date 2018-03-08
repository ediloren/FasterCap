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


// AutoPanel.cpp : autopanel class
//
// Basic panel class, to be used in AutoRefine


//#include "stdafx.h"

#include "AutoPanel.h"

// link with FasterCap main frame
#include "../FasterCapGlobal.h"

#define AUTOPANEL_2_DIV_FOURTHROOT_27   0.87738267530166164054614593453133
#define AUTOPANEL_2_DIV_FOURTHROOT_3    1.5196713713031850946623755013091
// -1/2+sqr(2)/4
#define AUTOPANEL_MINUS_0_5_PLUS_SQR2_DIV_4     -0.14644660940672623779957781894758
// 1/2+sqr(2)/4
#define AUTOPANEL_0_5_PLUS_SQR2_DIV_4           0.85355339059327376220042218105242
// sqr(2)
#define AUTOPANEL_SQR2                          1.4142135623730950488016887242097


///////////////////////
// CAutoPanel
///////////////////////

double CAutoPanel::CalcPanelGeomPar()
{
	C3DVector side[3];
	double normalMod, sidelen[3], cosmin;
	unsigned char a, b, c;

	// compute panel geometrical parameters
	//

	side[0] = m_clsVertex[1] - m_clsVertex[0];
	side[1] = m_clsVertex[2] - m_clsVertex[1];
	side[2] = m_clsVertex[0] - m_clsVertex[2];

	sidelen[0] = Mod(side[0]);
	sidelen[1] = Mod(side[1]);
	sidelen[2] = Mod(side[2]);

	normalMod = CalculateNormal(m_clsNormal);

	// compute panel area
	m_dDimension = normalMod / 2.0;

	// compute centroid
	m_clsCentroid = (m_clsVertex[0] + m_clsVertex[1] + m_clsVertex[2]) / 3.0;


	// compare and store max relative side sizes, and get index of min side
	b = MaxSide(sidelen[0], sidelen[1], sidelen[2]);

	//
	// return cos(min angle), for thin triangle check
	//

	// calculate a and c over circular reference to the sides
	if(b == 0) {
		a = 2;
		c = 1;
	}
	else if (b == 1) {
		a = 0;
		c = 2;
	}
	else {
		a = 1;
		c = 0;
	}
	// calculate cos(min angle) using cosine rule (angles are opposite to correspondent sides)
	//  b² = a² + c² -  2ac cosB
	// for reference, there is also sine rule (not used here):
	// a / sin(A) = b / sin(B) = c / sin(C)
	cosmin = (sidelen[a]*sidelen[a] + sidelen[c]*sidelen[c] - sidelen[b]*sidelen[b]) / (2.0 * sidelen[a] * sidelen[c]);

	return cosmin;
}


char CAutoPanel::MaxSide(double lside0, double lside1, double lside2)
{
	char minside;

	// find out max lenght panel sides, and index of min side

	if( (lside0 >= lside1) && (lside0 >= lside2) ) {
		m_ucMaxSide = 0;
		m_dMaxSideLen = lside0;
		if(lside1 >= lside2) {
			minside = 2;
		}
		else {
			minside = 1;
		}
	}
	else if( (lside1 >= lside0) && (lside1 >= lside2) ) {
		m_ucMaxSide = 1;
		m_dMaxSideLen = lside1;
		if(lside0 >= lside2) {
			minside = 2;
		}
		else {
			minside = 0;
		}
	}
	else {
		m_ucMaxSide = 2;
		m_dMaxSideLen = lside2;
		if(lside0 >= lside1) {
			minside = 1;
		}
		else {
			minside = 0;
		}
	}

	return minside;
}

// normal to be updated is passed as an argument, as in some cases this is not
// the 'm_clsNormal' member of the CAutoPanel class (see GetGeoNormal() )
double CAutoPanel::CalculateNormal(C3DVector_float &normal)
{
	double normalMod;
	unsigned char i, j;


	// compute unit vector normal to triangle plane
	// using Newell's method (see e.g. graphics gems III, V.5)
	// Plane normal direction is so that N x point1-point3
	// is directed inside the polygon
	//
	//         3                    O
	//         |\                   O
	//    N    | \                  O
	//     \   |  \                 O
	//      \  |   \                O
	//       \ |    \               O
	//        \|     \              O
	//         1------2             O
	//
	// Side benefit: module of the normal is 2xArea (like for cross(side1, side2) )
	// This is because the Newell method
	// calculates the normal from the areas of the projections of the polygon on the three
	// cartesian planes. The projected areas are computed as the sum of the signed areas
	// of the trapezoidal regions enclosed between the each edge of the projected polygon
	// and its projection onto the cartesian axes.
	// The 2x coefficient is only because, to speed up the calculation, the trapezoidal
	// area formula that would require division by 2 is not used, in favor of a simple
	// multiplication with no division by 2, as only the proportions are important for the
	// polygon normal calculation.

	// clear the normal vector
	normal.pos(0.0, 0.0, 0.0);
	// scan all points
	for(i = 0; i < 3; i++) {

		// compute cyclic index of next vertex
		j = i+1;
		if(j>2) {
			j = 0;
		}

		// calculate components of normal vector, using the
		// fact that the area of the trapezoids projected on the
		// xy, xz, yz planes is proportional to the z, y, x
		// components of the normal
		normal.x += (m_clsVertex[i].y - m_clsVertex[j].y) * (m_clsVertex[i].z + m_clsVertex[j].z);
		normal.y += (m_clsVertex[i].z - m_clsVertex[j].z) * (m_clsVertex[i].x + m_clsVertex[j].x);
		normal.z += (m_clsVertex[i].x - m_clsVertex[j].x) * (m_clsVertex[i].y + m_clsVertex[j].y);
	}
	// no need to check the actual length of the vector, since if the points are almost collinear,
	// the check below on cosmin will return a thin angle detection
	normalMod = normal.Mod();
	if(normalMod > AUTOPANEL_EPS) {
		normal /= normalMod;
	}

	return normalMod;
}

void CAutoPanel::MakeSuperPanel(CAutoPanel *leftSubPanel, CAutoPanel *rightSubPanel)
{
	C3DVector lNormal, rNormal, normal, cmToCs, xversor, yversor;
	double modnormal, R, a;

	// if subpanels are not coherent
	if( ((rightSubPanel->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL &&
	        (leftSubPanel->m_ucType & AUTOPANEL_IS_DIEL) != AUTOPANEL_IS_DIEL) ||
	        ((rightSubPanel->m_ucType & AUTOPANEL_IS_DIEL) != AUTOPANEL_IS_DIEL &&
	         (leftSubPanel->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) ) {

		ASSERT(true);
		ErrMsg("Internal error: dielectric subpanels %x, %x not coherent in type\n",
		       rightSubPanel, leftSubPanel);
	}

	// define the correct type
	if((rightSubPanel->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL &&
	        (leftSubPanel->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
		m_ucType |= (AUTOPANEL_IS_SUPER_NODE | AUTOPANEL_IS_DIEL | AUTOPANEL_OUTPERM_NORMAL_DIR);
	}
	else {
		m_ucType |= AUTOPANEL_IS_SUPER_NODE;
	}

	// not a leaf any more
	SetNotLeaf();
	m_pLeft = leftSubPanel;
	m_pRight = rightSubPanel;

#ifdef DEBUG_DUMP_BASIC
	m_pCond = leftSubPanel->m_pCond;
	m_pParent = this;
	m_pParent = this;
#endif

	// construction of the new super panel starts from the centroid.
	// This is the center of mass, i.e. the weighted centroid.
	// c3 = (c1*A1 + c2*A2) / (A1+A2)
	// (so when we calculate c4 = (c4*A4 + c3*A3) / (A4+A3) = (c4*A4 + ((c1*A1 + c2*A2) / (A1+A2))*(A1+A2)) / (A4+A1+A2) =
	// = (c4*A4 + c1*A1 + c2*A2) / (A4+A1+A2) this is also weighted at next level)
	//
	// compute super panel area
	m_dDimension = leftSubPanel->m_dDimension + rightSubPanel->m_dDimension;
	// compute centroid
	m_clsCentroid = (leftSubPanel->GetCentroid() * leftSubPanel->m_dDimension + rightSubPanel->GetCentroid() * rightSubPanel->m_dDimension) / m_dDimension;

	// now calculate the normal direction as averaged panel normal
	// Remark: do not normalize, otherwise when building higher level super nodes,
	// average will be wrong in respect to lower level nodes / leaves
	// (e.g. two panels of area 1/2 and one panel of area 1 with opposite orientation)
	//
	// first check if there was a reverted normal
	if( (leftSubPanel->m_ucType & AUTOPANEL_OUTPERM_NORMAL_DIR) == AUTOPANEL_OUTPERM_NORMAL_DIR) {
		lNormal = leftSubPanel->m_clsNormal;
	}
	else {
		lNormal = -leftSubPanel->m_clsNormal;
	}
	if( (rightSubPanel->m_ucType & AUTOPANEL_OUTPERM_NORMAL_DIR) == AUTOPANEL_OUTPERM_NORMAL_DIR) {
		rNormal = rightSubPanel->m_clsNormal;
	}
	else {
		rNormal = -rightSubPanel->m_clsNormal;
	}
	m_clsNormal = (lNormal * leftSubPanel->m_dDimension + rNormal * rightSubPanel->m_dDimension) / m_dDimension;
	// now check the length of the normal (could be very small or even zero for opposite normals;
	// in this case, considering the limit of two vectors that move aligning with opposite directions
	// along a line, the supporting plane of the super triangle contains the line
	modnormal = Mod(m_clsNormal);
	if(modnormal < AUTOPANEL_EPS) {
		// must pick a direction to fix one of the infinite possible planes with the same support
		// of the line formed by the normals. Let's solve:
		// DotProd(v1, v2) = x1*x2+y1*y2+z1*z2 = 0.0
		// so if x1 != 0,  x2 = (-y1*y2-z1*z2) / x1 and we assume y2=1 z2=1
		// etc for y1 and z1
		if(fabs(lNormal.x) > AUTOPANEL_EPS ) {
			normal.x = (-lNormal.y -lNormal.z) / lNormal.x;
			normal.y = 1.0f;
			normal.z = 1.0f;
		}
		else if(fabs(lNormal.y) > AUTOPANEL_EPS ) {
			normal.y = (-lNormal.x -lNormal.z) / lNormal.y;
			normal.x = 1.0f;
			normal.z = 1.0f;
		}
		else if(fabs(lNormal.z) > AUTOPANEL_EPS ) {
			normal.z = (-lNormal.x -lNormal.y) / lNormal.z;
			normal.x = 1.0f;
			normal.y = 1.0f;
		}
		else {
			// this can happen if both left and right subpanels are super panels, and by chance
			// both of them have zero dielectric normals.
			// in this case, compute the normal from the geometric normal (only chance)
			normal = leftSubPanel->GetGeoNormal() * leftSubPanel->m_dDimension + rightSubPanel->GetGeoNormal() * rightSubPanel->m_dDimension;
			normal.Normalize();
			// if still no chance (also geometric normals are opposite), take an arbitrary direction
			if(Mod(normal) < AUTOPANEL_EPS) {
				normal.x = 0.0f;
				normal.y = 0.0f;
				normal.z = 1.0f;
			}
		}
	}
	else {
		normal = m_clsNormal / modnormal;
	}

	// and then the vertexes
	//

	// we want the pointed side of the triangle directed towards the smallest of the two triangles
	//
	// so find first the plane containing the normal and the segment connecting the center of the super panel
	// to the center of the smallest of the two (left, right) panels
	//
	// this is the vector connecting the super panel centroid to the centroid of the smallest area panel
	if(leftSubPanel->m_dDimension < rightSubPanel->m_dDimension) {
		cmToCs = leftSubPanel->GetCentroid() - m_clsCentroid;
	}
	else {
		cmToCs = rightSubPanel->GetCentroid() - m_clsCentroid;
	}

	// now, it can happen that this vector and the normal are parallel
	// (e.g. when there are parallel, aligned panels)
	// in this case it does not matter where the super panel triangle points to,
	// just define an arbitrary vector, as long as different from the normal;
	// so rotate the normal 45 degrees along any angle
	// [cos(te)cos(psi)  -cos(fi)sin(psi)+sin(fi)sin(te)cos(psi)  sin(fi)sin(psi)+cos(fi)sin(te)cos(psi);
	//  cos(te)sin(psi)   cos(fi)cos(psi)+sin(fi)sin(te)sin(psi) -sin(fi)cos(psi)+cos(fi)sin(te)sin(psi);
	//  -sin(te)          sin(fi)cos(te)                          cos(fi)cos(te)                         ]
	// that is
	// [1/2     -1/2+sqr(2)/4   1/2+sqr(2)/4;
	//  1/2      1/2+sqr(2)/4  -1/2+sqr(2)/4;
	//  -sqr(2)  1/2            1/2]
	if( Mod(CrossProd(normal, cmToCs)) < AUTOPANEL_EPS) {
		cmToCs.x = normal.x * 0.5 + normal.y * (AUTOPANEL_MINUS_0_5_PLUS_SQR2_DIV_4) + normal.z * AUTOPANEL_0_5_PLUS_SQR2_DIV_4;
		cmToCs.y = normal.x * 0.5 + normal.y * AUTOPANEL_0_5_PLUS_SQR2_DIV_4 + normal.z * (AUTOPANEL_MINUS_0_5_PLUS_SQR2_DIV_4);
		cmToCs.z = normal.x * (-AUTOPANEL_SQR2) + normal.y * 0.5 + normal.z * 0.5;
	}

	// now we have the two vectors defining the plane perpendicular to the new panel
	// and contaning the centroid as well as the reference vertex.
	// We'll center a local coordinate system on the centroid; the local x versor
	// is the unit vector perpendicular to the plane defined by the two vectors
	// (remember that these two vectors does NOT usually span a 90 deg angle)
	xversor = CrossProd(cmToCs, normal);
	xversor.Normalize();
	// now y versor is simply perpendicular to both normal and xversor
	yversor = CrossProd(normal, xversor);
	// should already be a unit vector, but in any case..
	yversor.Normalize();

	// Let's find the position of the reference vertex. For an equilateral triangle
	// of area A, the radius of the circumcircle is R = 2*sqrt(A) / (27^(1/4))
	// (this is because R = a / sqrt(3) and A = sqrt(3)*a^2/4 where a is the triangle side)
	R = AUTOPANEL_2_DIV_FOURTHROOT_27 * sqrt(m_dDimension);
	m_clsVertex[0] = m_clsCentroid + yversor * R;
	// while the remaining two vertexes are chosed such as that
	// N x (m_clsVertex[0] - m_clsVertex[2]) is pointing inside the triangle
	// (in line with the choice we have done when calculating a panel normal).
	// Note that centroid to base side of the triangle distance is R / 2
	a = AUTOPANEL_2_DIV_FOURTHROOT_3 * sqrt(m_dDimension);
	m_clsVertex[1] = m_clsCentroid - yversor * (R/2.0) - xversor * (a/2.0);
	m_clsVertex[2] = m_clsCentroid - yversor * (R/2.0) + xversor * (a/2.0);

	// and finally mark max side (any of the three, the triangle is equilateral)
	m_ucMaxSide = 0;
	m_dMaxSideLen = a;
}

int CAutoPanel::Subdivide()
{
	C3DVector midpoint;
	char maxSide;
	double newsidelen, halfsidelen, leftSideLen[3], rightSideLen[3];

	if(g_bFCContinue == false) {
		return FC_USER_BREAK;
	}

	// if not a leaf, the panel is already subdivided
	if(IsLeaf() != true) {
		return FC_NORMAL_END;
	}

	// create new child panels (newly created panels are leaves by default)
	// SAFENEW_RET(TYPE, VAR, MEM)
	SAFENEW_RET(CAutoPanel, m_pLeft, g_clsMemUsage.m_ulPanelsMem)
	SAFENEW_RET(CAutoPanel, m_pRight, g_clsMemUsage.m_ulPanelsMem)


#ifdef DEBUG_DUMP_BASIC
	GetLeftChild()->m_pParent = this;
	GetRightChild()->m_pParent = this;
#endif

	maxSide = m_ucMaxSide;

	if(maxSide == 0) {

		midpoint = (m_clsVertex[1] + m_clsVertex[0]) / 2.0;

		newsidelen = Mod(m_clsVertex[2] - midpoint);
		halfsidelen = m_dMaxSideLen / 2.0;

		GetLeftChild()->m_clsVertex[0] = m_clsVertex[0];
		GetLeftChild()->m_clsVertex[1] = midpoint;
		GetLeftChild()->m_clsVertex[2] = m_clsVertex[2];

		leftSideLen[0] = halfsidelen;
		leftSideLen[1] = newsidelen;
		leftSideLen[2] = Mod(m_clsVertex[0] - m_clsVertex[2]);

		GetRightChild()->m_clsVertex[0] = midpoint;
		GetRightChild()->m_clsVertex[1] = m_clsVertex[1];
		GetRightChild()->m_clsVertex[2] = m_clsVertex[2];

		rightSideLen[0] = halfsidelen;
		rightSideLen[1] = Mod(m_clsVertex[2] - m_clsVertex[1]);
		rightSideLen[2] = newsidelen;
	}
	else if(maxSide == 1) {

		midpoint = (m_clsVertex[2] + m_clsVertex[1])/2;

		newsidelen = Mod(m_clsVertex[0] - midpoint);
		halfsidelen = m_dMaxSideLen / 2.0;

		GetLeftChild()->m_clsVertex[0] = m_clsVertex[1];
		GetLeftChild()->m_clsVertex[1] = midpoint;
		GetLeftChild()->m_clsVertex[2] = m_clsVertex[0];

		leftSideLen[0] = halfsidelen;
		leftSideLen[1] = newsidelen;
		leftSideLen[2] = Mod(m_clsVertex[1] - m_clsVertex[0]);

		GetRightChild()->m_clsVertex[0] = midpoint;
		GetRightChild()->m_clsVertex[1] = m_clsVertex[2];
		GetRightChild()->m_clsVertex[2] = m_clsVertex[0];

		rightSideLen[0] = halfsidelen;
		rightSideLen[1] = Mod(m_clsVertex[0] - m_clsVertex[2]);
		rightSideLen[2] = newsidelen;
	}
	else {

		midpoint = (m_clsVertex[0] + m_clsVertex[2])/2;

		newsidelen = Mod(m_clsVertex[1] - midpoint);
		halfsidelen = m_dMaxSideLen / 2.0;

		GetLeftChild()->m_clsVertex[0] = m_clsVertex[2];
		GetLeftChild()->m_clsVertex[1] = midpoint;
		GetLeftChild()->m_clsVertex[2] = m_clsVertex[1];

		leftSideLen[0] = halfsidelen;
		leftSideLen[1] = newsidelen;
		leftSideLen[2] = Mod(m_clsVertex[2] - m_clsVertex[1]);

		GetRightChild()->m_clsVertex[0] = midpoint;
		GetRightChild()->m_clsVertex[1] = m_clsVertex[0];
		GetRightChild()->m_clsVertex[2] = m_clsVertex[1];

		rightSideLen[0] = halfsidelen;
		rightSideLen[1] = Mod(m_clsVertex[1] - m_clsVertex[0]);
		rightSideLen[2] = newsidelen;
	}

	// copy panel type
	GetLeftChild()->m_ucType = m_ucType;
	GetRightChild()->m_ucType = m_ucType;

	// copy diel constant index
	GetLeftChild()->m_ucDielIndex = m_ucDielIndex;
	GetRightChild()->m_ucDielIndex = m_ucDielIndex;

	// not a leaf any more (remark: modifies m_ucType !)
	SetNotLeaf();

#ifdef DEBUG_DUMP_BASIC
	GetLeftChild()->m_pCond = m_pCond;
	GetRightChild()->m_pCond = m_pCond;

	GetLeftChild()->m_iLevel = panel->m_iLevel + 1;
	GetRightChild()->m_iLevel = panel->m_iLevel + 1;
#endif

	GetLeftChild()->m_dDimension = m_dDimension / 2.0;
	GetRightChild()->m_dDimension = m_dDimension / 2.0;

	GetLeftChild()->m_clsNormal = m_clsNormal;
	GetRightChild()->m_clsNormal = m_clsNormal;

	// extending the already computed charge (if any)
//	GetLeftChild()->m_dCharge = panel->m_dCharge / 2.0;
//	GetRightChild()->m_dCharge = panel->m_dCharge / 2.0;
	// extending the already computed charge density (if any)
	GetLeftChild()->m_dCharge = m_dCharge;
	GetRightChild()->m_dCharge = m_dCharge;
	// and copy the potential
	GetLeftChild()->m_dPotential = m_dPotential;
	GetRightChild()->m_dPotential = m_dPotential;

	// compute centroid
	GetLeftChild()->m_clsCentroid = (GetLeftChild()->m_clsVertex[0] + GetLeftChild()->m_clsVertex[1] + GetLeftChild()->m_clsVertex[2]) / 3.0;
	GetRightChild()->m_clsCentroid = (GetRightChild()->m_clsVertex[0] + GetRightChild()->m_clsVertex[1] + GetRightChild()->m_clsVertex[2]) / 3.0;

	GetLeftChild()->MaxSide(leftSideLen[0], leftSideLen[1], leftSideLen[2]);
	GetRightChild()->MaxSide(rightSideLen[0], rightSideLen[1], rightSideLen[2]);

	return FC_NORMAL_END;
}

void CAutoPanel::ErrorPrintCoords()
{
	ErrMsg("       panel %lx corner coordinates are: (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", this,
	       m_clsVertex[0].x, m_clsVertex[0].y, m_clsVertex[0].z,
	       m_clsVertex[1].x, m_clsVertex[1].y, m_clsVertex[1].z,
	       m_clsVertex[2].x, m_clsVertex[2].y, m_clsVertex[2].z);
	ErrMsg("       centroid is : (%g,%g,%g)\n",
	       GetCentroid().x, GetCentroid().y, GetCentroid().z);
}

///////////////////////
// CAutoQPanel
///////////////////////


double CAutoQPanel::CalcPanelGeomPar()
{
	C3DVector side[4];
	double normalMod, sidelen[4], cosmin, cosangle;
	int i;

	// compute panel geometrical parameters
	//

	side[0] = m_clsQVertex[1] - m_clsQVertex[0];
	side[1] = m_clsQVertex[2] - m_clsQVertex[1];
	side[2] = m_clsQVertex[3] - m_clsQVertex[2];
	side[3] = m_clsQVertex[0] - m_clsQVertex[3];

	sidelen[0] = Mod(side[0]);
	sidelen[1] = Mod(side[1]);
	sidelen[2] = Mod(side[2]);
	sidelen[3] = Mod(side[3]);

	normalMod = CalculateNormal(m_clsNormal);

	// compute panel area
	m_dDimension = normalMod / 2.0;

	// compute centroid
	m_clsCentroid = (m_clsQVertex[0] + m_clsQVertex[1] + m_clsQVertex[2] + m_clsQVertex[3]) / 4.0;


	// compare and store max relative side sizes
	MaxSide(sidelen);

	//
	// return cos(min angle), for thin quadrilateral check
	//

	// calculate cos(min angle)
	cosmin = 1.0;
	for(i=0; i<4; i++) {
        cosangle = DotProd(side[i], side[(i+1)%4]) / (sidelen[i] * sidelen[(i+1)%4]);
        if(cosangle < cosmin) {
            cosmin = cosangle;
        }
	}
	return cosmin;
}

char CAutoQPanel::MaxSide(double lside[4])
{
	unsigned char minside;
	int i;

	// find out max lenght panel sides, and index of min side
	m_ucMaxSide = 0;
	minside = 0;
    for (i=1; i<4; i++) {
        if(lside[i] >= lside[m_ucMaxSide]) {
            m_ucMaxSide = i;
        }
        if(lside[i] < lside[minside]) {
            minside = i;
        }
    }
    m_dMaxSideLen = lside[m_ucMaxSide];

	return minside;
}

double CAutoQPanel::CalculateNormal(C3DVector_float &normal)
{
	double normalMod;
	unsigned char i, j;


	// compute unit vector normal of the polygon
	// using Newell's method (see e.g. graphics gems III, V.5)
	// Plane normal direction is so that N x point1-point3
	// is directed inside the polygon.
	// For nonplanar polygons, Newell’s method computes a “best-fit” normal
	//
	//         3                    O
	//         |\                   O
	//    N    | \                  O
	//     \   |  \                 O
	//      \  |   \                O
	//       \ |    \               O
	//        \|     \              O
	//         1------2             O
	//
	// Side benefit: module of the normal is 2xArea. This is because the Newell method
	// calculates the normal from the areas of the projections of the polygon on the three
	// cartesian planes. The projected areas are computed as the sum of the signed areas
	// of the trapezoidal regions enclosed between the each edge of the projected polygon
	// and its projection onto the cartesian axes.
	// The 2x coefficient is only because, to speed up the calculation, the trapezoidal
	// area formula that would require division by 2 is not used, in favor of a simple
	// multiplication with no division by 2, as only the proportions are important for the
	// polygon normal calculation.

	// clear the normal vector
	normal.pos(0.0, 0.0, 0.0);
	// scan all points
	for(i = 0; i < 4; i++) {

		// compute cyclic index of next vertex
		j = i+1;
		if(j>3) {
			j = 0;
		}

		// calculate components of normal vector, using the
		// fact that the area of the trapezoids projected on the
		// xy, xz, yz planes is proportional to the z, y, x
		// components of the normal
		normal.x += (m_clsQVertex[i].y - m_clsQVertex[j].y) * (m_clsQVertex[i].z + m_clsQVertex[j].z);
		normal.y += (m_clsQVertex[i].z - m_clsQVertex[j].z) * (m_clsQVertex[i].x + m_clsQVertex[j].x);
		normal.z += (m_clsQVertex[i].x - m_clsQVertex[j].x) * (m_clsQVertex[i].y + m_clsQVertex[j].y);
	}
	// no need to check the actual length of the vector, since if the points are almost collinear,
	// the check below on cosmin will return a thin angle detection
	normalMod = normal.Mod();
	if(normalMod > AUTOPANEL_EPS) {
		normal /= normalMod;
	}

	return normalMod;
}

void CAutoQPanel::ErrorPrintCoords()
{
	ErrMsg("       Qpanel %lx corner coordinates are: (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", this,
	       m_clsQVertex[0].x, m_clsQVertex[0].y, m_clsQVertex[0].z,
	       m_clsQVertex[1].x, m_clsQVertex[1].y, m_clsQVertex[1].z,
	       m_clsQVertex[2].x, m_clsQVertex[2].y, m_clsQVertex[2].z,
	       m_clsQVertex[3].x, m_clsQVertex[3].y, m_clsQVertex[3].z);
//	ErrMsg("       centroid is : (%g,%g,%g)\n",
//	       GetCentroid().x, GetCentroid().y, GetCentroid().z);
}

