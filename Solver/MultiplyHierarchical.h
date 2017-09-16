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


// MultiplyHierarchical.h : hierarchical multiplication class header file
//

#if !defined(AFX_MULTIPLYHIERARCHICAL_H__ENRI1945_5487_55D2_9284_04F014CF5600__INCLUDED_)
#define AFX_MULTIPLYHIERARCHICAL_H__ENRI1945_5487_55D2_9284_04F014CF5600__INCLUDED_

#include "SolverGlobal.h"

#include "Autorefine.h"

// includes for LinALg
#include "LinAlgebra/Vect.h"

#define MULTHIER_MAX_RECURS_DEPTH	128


class CMultHier : public CAutoRefine
{

public:
	CMultHier();
	~CMultHier();
	int AllocateMemory();
	void DeallocateMemory();
	int MultiplyMatByVec_fast(CLin_Vector *v, CLin_Vector *q);
	void CopyChargesToVec(CLin_Vector *q);
	void CopyVecToCharges(CLin_Vector *q);
	void InitFlatLinks();

	unsigned long m_ulFirstCondElemIndex;

protected:
	void ComputePanelCharges_fast();
	int ComputePanelPotentials_2fast();
	void ComputeLeafPotentials_fast();
	void CopyPanelCharges();
	void CopyCharges(CAutoElement* panel);
	void CopyVecCharges();
	void CopyVec(CAutoElement* panel);

	CLin_Range m_clsChargeVect, m_clsPotVect;
	long m_dIndex;
	CAutoElement *m_clsRecursVec[MULTHIER_MAX_RECURS_DEPTH];

};

#endif //!defined(AFX_MULTIPLYHIERARCHICAL_H__ENRI1945_5487_55D2_9284_04F014CF5600__INCLUDED_)
