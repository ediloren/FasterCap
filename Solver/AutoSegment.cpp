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


// Basic 2D segment class, to be used in AutoRefine
// E. Di Lorenzo, 2013/02/02

//#include "stdafx.h"

#include "AutoSegment.h"

// link with FasterCap main frame
#include "../FasterCapGlobal.h"

void CAutoSegment::CalcSegmentGeomPar()
{
	C2DVector segment;

	// compute panel geometrical parameters
	//

	segment = m_clsVertex[1] - m_clsVertex[0];

    // compute segment length
    m_dDimension = Mod(segment);
    // compute centroid
	m_clsCentroid = (m_clsVertex[0] + m_clsVertex[1]) / 2.0;

	// compute normal
    CalculateNormal(m_clsNormal);
}

double CAutoSegment::CalculateNormal(C2DVector_float &normal)
{
	double normalMod;

    // y is normal to the segment and the smallest angle
    // is in ccw direction from x
	normal.x = m_clsVertex[0].y - m_clsVertex[1].y;
	normal.y = m_clsVertex[1].x - m_clsVertex[0].x;

	// clear the normal vector
	normalMod = normal.Mod();
	if(normalMod > AUTOPANEL_EPS) {
		normal /= normalMod;
	}

	return normalMod;
}

void CAutoSegment::MakeSuperSegment(CAutoSegment *leftSubPanel, CAutoSegment *rightSubPanel)
{
	C2DVector lNormal, rNormal, normal, segment;
	double modnormal;

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

	// construction of the new super segment starts from the centroid.
	// This is the center of mass, i.e. the weighted centroid.
	// c3 = (c1*A1 + c2*A2) / (A1+A2)
	// (so when we calculate c4 = (c4*A4 + c3*A3) / (A4+A3) = (c4*A4 + ((c1*A1 + c2*A2) / (A1+A2))*(A1+A2)) / (A4+A1+A2) =
	// = (c4*A4 + c1*A1 + c2*A2) / (A4+A1+A2) this is also weighted at next level)
	//
	// compute super panel area
	m_dDimension = leftSubPanel->m_dDimension + rightSubPanel->m_dDimension;
	// compute centroid
	m_clsCentroid = (leftSubPanel->GetCentroid() * leftSubPanel->m_dDimension + rightSubPanel->GetCentroid() * rightSubPanel->m_dDimension) / m_dDimension;

	// now calculate the normal direction as averaged segment normal
	// Remark: do not normalize, otherwise when building higher level super segments,
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
	// along a line, the supporting line of the super segment contains the normal
	modnormal = Mod(m_clsNormal);
	if(modnormal < AUTOPANEL_EPS) {
		// must pick a direction of one of the two normals; let's take the first
		normal = lNormal;
	}
	else {
		normal = m_clsNormal / modnormal;
	}

	// and then the vertexes
	//
	// first a unit vector along the segment
    segment.x = normal.y;
    segment.y = -normal.x;
    // then the vertexes
    m_clsVertex[0] = m_clsCentroid - segment * (m_dDimension / 2.0);
    m_clsVertex[1] = m_clsCentroid + segment * (m_dDimension / 2.0);
}

int CAutoSegment::Subdivide()
{
	C2DVector midpoint;

	if(g_bFCContinue == false) {
		return FC_USER_BREAK;
	}

	// if not a leaf, the panel is already subdivided
	if(IsLeaf() != true) {
		return FC_NORMAL_END;
	}

	// create new child panels (newly created panels are leaves by default)
	// SAFENEW_RET(TYPE, VAR, MEM)
	SAFENEW_RET(CAutoSegment, m_pLeft, g_clsMemUsage.m_ulPanelsMem)
	SAFENEW_RET(CAutoSegment, m_pRight, g_clsMemUsage.m_ulPanelsMem)


#ifdef DEBUG_DUMP_BASIC
	GetLeftChild()->m_pParent = this;
	GetRightChild()->m_pParent = this;
#endif

    // geometric parameters

    midpoint = (m_clsVertex[1] + m_clsVertex[0]) / 2.0;

    GetLeftChild()->m_clsVertex[0] = m_clsVertex[0];
	GetLeftChild()->m_clsVertex[1] = midpoint;

	GetRightChild()->m_clsVertex[0] = midpoint;
	GetRightChild()->m_clsVertex[1] = m_clsVertex[1];

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
	GetLeftChild()->m_clsCentroid = (GetLeftChild()->m_clsVertex[0] + GetLeftChild()->m_clsVertex[1]) / 2.0;
	GetRightChild()->m_clsCentroid = (GetRightChild()->m_clsVertex[0] + GetRightChild()->m_clsVertex[1]) / 2.0;

	return FC_NORMAL_END;
}

// scale the segment geometrical data by the 'scale' factor
void CAutoSegment::Scale(double scale)
{
	m_dDimension *= scale;
	m_clsVertex[0] *= scale;
	m_clsVertex[1] *= scale;
	m_clsCentroid *= scale;

}

void CAutoSegment::ErrorPrintCoords()
{
	ErrMsg("       segment %lx end point coordinates are: (%g,%g) (%g,%g)\n", this,
	       m_clsVertex[0].x, m_clsVertex[0].y,
	       m_clsVertex[1].x, m_clsVertex[1].y);
	ErrMsg("       centroid is : (%g,%g)\n",
	       GetCentroid().x, GetCentroid().y);
}

