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


// CAutoConductor class, to be used in AutoRefine
// E. Di Lorenzo, 2013/02/02

//#include "stdafx.h"

#include "AutoConductor.h"

// link with FasterCap main frame
#include "../FasterCapGlobal.h"


CAutoConductor::CAutoConductor()
{
	m_sName[0] = '\0';
	m_ulInputPanelNum = 0;
	m_ulLeafPanelNum = 0;
	m_uTopElement.m_pTopPanel = NULL;
	m_ucMaxSurfOutperm = 0;
}

CAutoConductor::CAutoConductor(char *name, bool isdiel, double outpermRe, double outpermIm)
{
	strcpy(m_sName, name);
	m_ulInputPanelNum = 0;
	m_ulLeafPanelNum = 0;
	m_bIsDiel = isdiel;
	m_dOutperm[0] = outpermRe;
	m_dOutperm[1] = outpermIm;
	m_dInperm[0] = 0.0;
	m_dInperm[1] = 0.0;
	m_clsDielRef3DPoint = C3DVector(0,0,0);
	m_uTopElement.m_pTopPanel = NULL;
	m_ucMaxSurfOutperm = 0;
}

CAutoConductor::CAutoConductor(char *name, bool isdiel, double outpermRe, double outpermIm,
                               double inpermRe, double inpermIm, C3DVector dielrefpoint)
{
	strcpy(m_sName, name);
	m_ulInputPanelNum = 0;
	m_ulLeafPanelNum = 0;
	m_bIsDiel = isdiel;
	m_dOutperm[0] = outpermRe;
	m_dOutperm[1] = outpermIm;
	m_dInperm[0] = inpermRe;
	m_dInperm[1] = inpermIm;
	m_clsDielRef3DPoint = dielrefpoint;
	m_uTopElement.m_pTopPanel = NULL;
	m_ucMaxSurfOutperm = 0;
}

CAutoConductor::CAutoConductor(char *name, bool isdiel, double outpermRe, double outpermIm,
                               double inpermRe, double inpermIm, C2DVector dielrefpoint)
{
	strcpy(m_sName, name);
	m_ulInputPanelNum = 0;
	m_ulLeafPanelNum = 0;
	m_bIsDiel = isdiel;
	m_dOutperm[0] = outpermRe;
	m_dOutperm[1] = outpermIm;
	m_dInperm[0] = inpermRe;
	m_dInperm[1] = inpermIm;
	// trick of storing the dielectric 2D ref point inside the 3D ref point,
	// nulling the third coordinate
	m_clsDielRef3DPoint.x = dielrefpoint.x;
	m_clsDielRef3DPoint.y = dielrefpoint.y;
	m_clsDielRef3DPoint.z = 0.0;
	m_uTopElement.m_pTopPanel = NULL;
	m_ucMaxSurfOutperm = 0;
}

// scale conductor
void CAutoConductor::Scale(const double scale)
{
    m_clsDielRef3DPoint *= scale;
    // Remark: formally we should scale the bbox as well, however no need,
    // since this will be calculated AFTER the panels are scaled
}
