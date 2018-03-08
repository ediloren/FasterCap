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


// AutoPanel.h : autopanel class header file
//

#if !defined(AFX_AUTOPANEL_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_AUTOPANEL_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_

#include "SolverGlobal.h"

#include <deque>
#include <string>

// for base class
#include "AutoElement.h"

// includes for vector class
#include "Geometry/Vector3D.h"

using namespace std;

class CAutoConductor;

// panel structure used by CAutoRefine
class CAutoPanel : public CAutoElement
{
public:
    // virtual functions
	virtual double CalcPanelGeomPar();

	// standard functions
	char MaxSide(double lside0, double lside1, double lside2);
	void MakeSuperPanel(CAutoPanel *leftSubPanel, CAutoPanel *rightSubPanel);
    int Subdivide();

    //
	// virtual functions implementation

	void ErrorPrintCoords();
    unsigned char GetClass()
    {
        return AUTOELEMENT_PANEL;
    }

    //
    // inlines

	inline C3DVector_float& GetCentroid() {
		return 	m_clsCentroid;
	}

	inline CAutoPanel *GetLeftChild() {
		return 	(CAutoPanel*)m_pLeft;
	}

	inline CAutoPanel *GetRightChild() {
		return 	(CAutoPanel*)m_pRight;
	}

	inline double GetMaxSideLen()
	{
		return m_dMaxSideLen;
	}

    // geometric normal (i.e. calculated from the panel vertexes) can be
    // different from the normal used for the dielectric calculations in two cases:
    // Case 1:
    // in case of super-panels. In particular, if a super panel is created
    // from two normal-opposite triangles, dielectric normal is zero,
    // even if the geometric normal is the one of the actual triangular panel
    // Case 2:
    // for diel normal in opposide direction w.r.t the geometric normal
    // (i.e. panel->m_ucType not flagged AUTOPANEL_OUTPERM_NORMAL_DIR)
    //
    // 'm_clsNormal' stores the gemetric normal for standard panels,
    // while for super panels it stores the averaged normal from sub-panels
	inline C3DVector GetGeoNormal()
	{
		C3DVector_float unitvector;

		if(m_ucType & AUTOPANEL_IS_SUPER_NODE) {
			// calculate on the fly
            CalculateNormal(unitvector);
		}
		else {
			unitvector = m_clsNormal;
		}

		return unitvector;
	}

	inline C3DVector GetDielNormal()
	{
		double modnormal;
		C3DVector unitvector;

		if(m_ucType & AUTOPANEL_IS_SUPER_NODE) {
			// must normalize since super panel normals are not unit vectors
			// TBC warning: could normalize all super panel normals
			// after the end of hierarchical super panel build pass, to avoid calculation every time
			modnormal = Mod(m_clsNormal);
			if(modnormal < AUTOPANEL_EPS) {
				unitvector = C3DVector(0,0,0);
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

	inline double GetArea()
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

	inline double GetMaxSideLen()
	{
		double lside0, lside1, lside2;
		C3DVector side1, side2, side3;

		if( (m_ucType & AUTOPANEL_IS_SUPER_NODE) == AUTOPANEL_IS_SUPER_NODE) {
			return m_dSideLen[0];
		}
		else {
			side1 = m_clsVertex[1] - m_clsVertex[0];
			side2 = m_clsVertex[2] - m_clsVertex[1];
			side3 = m_clsVertex[0] - m_clsVertex[2];

			//		lside0 = Mod(m_clsVertex[1] - m_clsVertex[0]);
			//		lside1 = Mod(m_clsVertex[2] - m_clsVertex[1]);
			//		lside2 = Mod(m_clsVertex[0] - m_clsVertex[2]);

			lside0 = Mod(side1);
			lside1 = Mod(side2);
			lside2 = Mod(side3);

			// find out max and mid lenght panel sides

			if( (lside0 >= lside1) && (lside0 >= lside2) ) {
				return lside0;
			}
			else if( (lside1 >= lside0) && (lside1 >= lside2) ) {
				return lside1;
			}
			else {
				return lside2;
			}
		}
	}
*/

	double m_dMaxSideLen;
	unsigned char m_ucMaxSide;
	C3DVector_float m_clsVertex[3];
	C3DVector_float m_clsNormal, m_clsCentroid;

protected:
    // virtual functions
    virtual double CalculateNormal(C3DVector_float &normal);

};

// Quadrilateral panel is derived from standard triangular panel
// The QPanel is used only for dumping input geometry; all other operations
// rely on the standard triangular panels CAutoPanel.
// In this case, adding 'm_clsQVertex' on top of 'm_clsVertex' array
// (wasting memory) is not a big issue, as the number of CAutoQPanel
// variables will be limited by definiton (no split will occur)
class CAutoQPanel : public CAutoPanel
{
public:
	// standard functions
	char MaxSide(double lside[4]);

    //
	// virtual functions implementation

	void ErrorPrintCoords();
    unsigned char GetClass()
    {
        return AUTOELEMENT_QPANEL;
    }
	double CalcPanelGeomPar();

	C3DVector_float m_clsQVertex[4];

protected:
    // virtual functions implementation
    double CalculateNormal(C3DVector_float &normal);
};

typedef std::deque<CAutoPanel*>  StlAutoPanelDeque;


#endif //!defined(AFX_AUTOPANEL_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
