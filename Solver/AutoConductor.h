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


// AutoConductor.h : CAutoConductor class definitions
// E. Di Lorenzo, 2013/02/02

#if !defined(AFX_AUTOCONDUCTOR_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_AUTOCONDUCTOR_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_

#include "AutoPanel.h"
#include "AutoSegment.h"
// includes for bbox
#include "Geometry/Geometry3D.h"
#include "Geometry/Geometry2D.h"

#include <deque>

using namespace std;

#define AUTOCONDUCTOR_MAX_NAME_LEN				512

class CInteraction
{
public:
	CAutoPanel* m_clsPanel;
	double m_dPotCoeff;

	CInteraction() {}
	CInteraction(CAutoPanel *intPanel, double intPotCoeff) : m_clsPanel(intPanel), m_dPotCoeff(intPotCoeff) { }
};

typedef C3DList<CInteraction> InteractionC3DList;

// conductor structure used by CAutoRefine
class CAutoConductor
{
public:
	CAutoConductor();
	CAutoConductor(char *name, bool isdiel, double outpermRe, double outpermIm);
	CAutoConductor(char *name, bool isdiel, double outpermRe, double outpermIm,
					double inpermRe, double inpermIm, C3DVector dielrefpoint);
	CAutoConductor(char *name, bool isdiel, double outpermRe, double outpermIm,
					double inpermRe, double inpermIm, C2DVector dielrefpoint);
    void Scale(const double scale);

	char m_sName[AUTOCONDUCTOR_MAX_NAME_LEN];
	// panels in 3D, segments in 2D
    StlAutoPanelDeque m_stlPanels;
    StlAutoSegmentDeque m_stlSegments;

	unsigned long m_ulLeafPanelNum;
	unsigned long m_ulInputPanelNum;
	bool m_bIsDiel;
	// diel constants for dielectrics
	double m_dOutperm[2], m_dInperm[2];
	// diel constants for conductors
	double m_dSurfOutperm[AUTOPANEL_MAX_DIEL_NUM][2];
	unsigned char m_ucMaxSurfOutperm;
	// 'm_clsDielRefxDPoint' is always on 'm_dOutperm' side
    C3DVector m_clsDielRef3DPoint;
	// conductor bounding box
    C3DBBox m_cls3DBbox;
	// top-level panel (either real panel or superpanel)
	union {
	    CAutoPanel *m_pTopPanel;
	    CAutoSegment *m_pTopSegment;
	    CAutoElement *m_pTopElement;
	} m_uTopElement;

protected:

};

typedef std::deque<CAutoConductor*>  StlAutoCondDeque;


#endif //!defined(AFX_AUTOCONDUCTOR_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
