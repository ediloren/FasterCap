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


// 2D segments & 3D patch base class, to be used in AutoRefine
// E. Di Lorenzo, 2013/02/05

#include "AutoElement.h"

// link with FasterCap main frame
#include "../FasterCapGlobal.h"

CAutoElement::CAutoElement()
{
	unsigned char i;

	// initialize panel structure

	m_ucType = AUTOPANEL_IS_LEAF;

	for(i=0; i<AUTOPANEL_MAX_NUM_OF_HIERARCHIES; i++) {
		m_ulLinkIndexStart[i] = 0;
		m_ulLinkIndexEnd[i] = 0;
	}
	m_pLeft = m_pRight = NULL;
	m_dCharge = 1.0;

#ifdef DEBUG_DUMP_BASIC
	m_pParent = NULL;
	m_iLevel = 0;
#endif

}

void CAutoElement::InitElementsTree()
{
    unsigned char i;

	// reset link counter
	for(i=0; i<AUTOPANEL_MAX_NUM_OF_HIERARCHIES; i++) {
		m_ulLinkIndexStart[i] = 0;
		m_ulLinkIndexEnd[i] = 0;
	}

	if(IsLeaf() != true) {
		// call recursively the routine for each child
		m_pLeft->InitElementsTree();
		m_pRight->InitElementsTree();
	}
}

void CAutoElement::DeleteElementsTree()
{
	// visit the tree

	if (IsLeaf() != true) {
		// scan the subtree, only after delete panel
		m_pLeft->DeleteElementsTree();
		m_pRight->DeleteElementsTree();
	}

	// then delete panel
	delete this;
}

