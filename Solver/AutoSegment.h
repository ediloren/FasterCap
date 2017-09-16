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


// AutoSegment.h : 2D segments definitions
// E. Di Lorenzo, 2013/02/02

#if !defined(AFX_AUTOSEGMENT_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_AUTOSEGMENT_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_

#include "SolverGlobal.h"

#include <deque>
#include <string>

// for base class
#include "AutoElement.h"

// includes for vector class
#include "Geometry/Vector2D.h"

// for common defines
#include "AutoPanel.h"

using namespace std;

class CAutoConductor;

// panel structure used by CAutoRefine
class CAutoSegment : public CAutoElement
{
public:
	void CalcSegmentGeomPar();
	void MakeSuperSegment(CAutoSegment *leftSubPanel, CAutoSegment *rightSubPanel);
    int Subdivide();
    void Scale(double scale);

    //
	// virtual functions implementation

	void ErrorPrintCoords();
    unsigned char GetClass()
    {
        return AUTOELEMENT_SEGMENT;
    }

    //
    // inlines

	inline C2DVector_float& GetCentroid() {
		return 	m_clsCentroid;
	}
	inline CAutoSegment *GetLeftChild() {
		return 	(CAutoSegment*)m_pLeft;
	}
	inline CAutoSegment *GetRightChild() {
		return 	(CAutoSegment*)m_pRight;
	}

    // geometric normal (i.e. calculated from the segment vertexes) can be
    // different from the normal used for the dielectric calculations in two cases:
    // Case 1:
    // in case of super-panels. In particular, if a super panel is created
    // from two normal-opposite segments, dielectric normal is zero,
    // even if the geometric normal is the one of the actual segment
    // Case 2:
    // for diel normal in opposide direction w.r.t the geometric normal
    // (i.e. panel->m_ucType not flagged AUTOPANEL_OUTPERM_NORMAL_DIR)
    //
    // 'm_clsNormal' stores the gemetric normal for standard segments,
    // while for super segments it stores the averaged normal from sub-segments
	inline C2DVector GetGeoNormal()
	{
		C2DVector_float unitvector;

		if(m_ucType & AUTOPANEL_IS_SUPER_NODE) {
			// calculate on the fly
            CalculateNormal(unitvector);
		}
		else {
			unitvector = m_clsNormal;
		}

		return unitvector;
	}

	inline C2DVector GetDielNormal()
	{
		double modnormal;
		C2DVector unitvector;

		if(m_ucType & AUTOPANEL_IS_SUPER_NODE) {
			// must normalize since super panel normals are not unit vectors
			// TBC warning: could normalize all super panel normals
			// after the end of hierarchical super panel build pass, to avoid calculating every time
			modnormal = Mod(m_clsNormal);
			if(modnormal < AUTOPANEL_EPS) {
				unitvector = C2DVector(0,0);
			}
			else {
				unitvector = m_clsNormal / modnormal;
			}
		}
		else {
			unitvector = m_clsNormal;
		}

        if( (m_ucType & AUTOPANEL_OUTPERM_NORMAL_DIR) != AUTOPANEL_OUTPERM_NORMAL_DIR ) {
            unitvector.Invert();
        }

		return unitvector;
	}

	inline double GetLength()
	{
		return m_dDimension;
	}

/*
	// for saving memory of CAutoPanel class structure

	inline C3DVector GetCentroid() {

		if( (m_ucType & AUTOPANEL_IS_SUPER_NODE) == AUTOPANEL_IS_SUPER_NODE) {
			// if super panel, compute centroid from bbox
			return (m_clsVertex[0] + m_clsVertex[1]) / 2.0;
		}
		else {
			return (m_clsVertex[0] + m_clsVertex[1] + m_clsVertex[2]) / 3.0;
		}
	}
*/

	C2DVector_float m_clsVertex[2];
	C2DVector_float m_clsNormal, m_clsCentroid;

protected:

    double CalculateNormal(C2DVector_float &normal);

};

typedef std::deque<CAutoSegment*>  StlAutoSegmentDeque;


#endif //!defined(AFX_AUTOSEGMENT_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
