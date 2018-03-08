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


// AutomationHelper.h : definition of MS Win Automation helper functions
// Enrico Di Lorenzo 2013/03/14

// These classes are wrappers for Automation variables.
// The main ideas are taken from MS Win SDK 2003, files
// afxdisp.h, afxole.inl, olevar.cpp
// but here the classes are re-implemented from scratch

#ifndef AUTOMATIONHELPER_DEFS
#define AUTOMATIONHELPER_DEFS

// definition of wxVariantData needed by the oleutils.h include
#include <wx/variant.h>
// for COM Automation
#include <windows.h>
#include <objbase.h>
#include <initguid.h>
#include <wx/msw/ole/oleutils.h>

// for wxASSERT
#include <wx/debug.h>

#define ASSERT wxASSERT


///////////////////////////////////
// AutoBSTR class (BSTR wrapper)
///////////////////////////////////

class AutoBSTR
{
public:
	AutoBSTR();
	~AutoBSTR();

	void Clear();

	const AutoBSTR& operator=(const BSTR& BSTRsrc);
	const AutoBSTR& operator=(const char *textSrc);
	operator BSTR();
	operator LPBSTR();

private:
    BSTR    m_clsBSTR;
};

inline AutoBSTR::AutoBSTR() { m_clsBSTR = NULL; }
inline AutoBSTR::~AutoBSTR() { if(m_clsBSTR != NULL) {SysFreeString(m_clsBSTR);} }
inline void AutoBSTR::Clear() { if(m_clsBSTR != NULL) {SysFreeString(m_clsBSTR);} }
inline AutoBSTR::operator BSTR() { return m_clsBSTR; }
inline AutoBSTR::operator LPBSTR() { return &m_clsBSTR; }


///////////////////////////////////
// AutoVariant class (VARIANT wrapper)
///////////////////////////////////

class AutoVariant : public VARIANT
{
public:
	AutoVariant();
	~AutoVariant();

	void Clear();
	void Attach(VARIANT& varSrc);
	VARIANT Detach();

	const AutoVariant& operator=(const VARIANT& varSrc);
	const AutoVariant& operator=(const AutoVariant& varSrc);
	const AutoVariant& operator=(const char* pCharSrc);
	const AutoVariant& operator=(double dblSrc);
	operator LPVARIANT();
};

inline AutoVariant::AutoVariant() { VariantInit(this); }
inline AutoVariant::~AutoVariant() { VariantClear(this); }
inline void AutoVariant::Clear() { VariantClear(this); }
inline AutoVariant::operator LPVARIANT() { return this; }


///////////////////////////////////
// AutoSafeArray class (SAFEARRAY wrapper)
///////////////////////////////////

class AutoSafeArray : public VARIANT
{
public:
	AutoSafeArray();
	~AutoSafeArray();

	void Clear();
	void Create(VARTYPE vtSrc, DWORD dwDims, DWORD* rgElements);
	void Create(VARTYPE vtSrc, DWORD dwDims, SAFEARRAYBOUND* rgsabounds);
	void PutElement(long* rgIndices, void* pvData);
	DWORD GetElemSize();

	operator LPVARIANT();
};

// helper function
void SafeArrayInit(AutoSafeArray* psa);

inline AutoSafeArray::AutoSafeArray()	{ SafeArrayInit(this); vt = VT_EMPTY; }
inline AutoSafeArray::~AutoSafeArray()	{ Clear(); }
inline void AutoSafeArray::Clear()	{ VariantClear(this); }
inline DWORD AutoSafeArray::GetElemSize() { return SafeArrayGetElemsize(parray); }
inline AutoSafeArray::operator LPVARIANT()	{ return this; }



#endif // AUTOMATIONHELPER_DEFS
