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


// AutoSegment.h : 2D segments & 3D patch base class definition
// E. Di Lorenzo, 2013/02/05

#if !defined(AFX_AUTOELEMENT_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_AUTOELEMENT_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_

#include "SolverGlobal.h"

#include <deque>
#include <string>

using namespace std;

// constant used to identify the class of the derived classes
// throught the virtual function GetClass()
#define AUTOELEMENT_SEGMENT     1
#define AUTOELEMENT_PANEL       2
#define AUTOELEMENT_QPANEL      4

// constants used in CAutoPanel type 'm_ucType'
#define AUTOPANEL_IS_LEAF					(unsigned char)1
#define AUTOPANEL_IS_DIEL					(unsigned char)2
#define AUTOPANEL_OUTPERM_NORMAL_DIR		(unsigned char)4
#define AUTOPANEL_OUTPERM_ELEMENT_LEVEL		(unsigned char)8
#define AUTOPANEL_IS_SUPER_NODE				(unsigned char)128

// maximum number of hierarchies supported by the panel structure
// (e.g. for hierarchical preconditioning)
#define AUTOPANEL_MAX_NUM_OF_HIERARCHIES	(unsigned char)2

// maximum number of different dielectric constants for the mediums surrounding
// the different surfaces belonging to a single conductor
#define AUTOPANEL_MAX_DIEL_NUM				256

// EPS used to avoid divisions by zero, and other use
#define AUTOPANEL_EPS 1E-12

// panel structure used by CAutoRefine
class CAutoElement
{
public:

    CAutoElement();
    virtual ~CAutoElement() {}

	inline bool IsLeaf()
	{
		return (m_ucType & AUTOPANEL_IS_LEAF);
	}

	inline void SetNotLeaf()
	{
		m_ucType &= ~(AUTOPANEL_IS_LEAF);
	}

	inline double GetDimension()
	{
		return m_dDimension;
	}

    void InitElementsTree();
    void DeleteElementsTree();

    // virtual functions (pure virtual)
	virtual void ErrorPrintCoords() = 0;
	virtual unsigned char GetClass() = 0;

	unsigned char m_ucType;
	unsigned long m_ulLinkIndexStart[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	unsigned long m_ulLinkIndexEnd[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	unsigned long m_lNumOfChildren;
	CAutoElement *m_pLeft, *m_pRight;
	double m_dCharge, m_dPotential;
	long m_lIndex[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	unsigned char m_ucDielIndex;
	double m_dDimension;


#ifdef DEBUG_DUMP_BASIC
	int m_iLevel;
	CAutoPanel *m_pParent;
	CAutoConductor *m_pCond;
#endif

};

#endif //!defined(AFX_AUTOELEMENT_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
