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


//  MS Win Automation helper classes implementation
//
//   Enrico Di Lorenzo, 2013/02/14
//

#include "AutomationHelper.h"

#include <exception>
using namespace std;

///////////////////////////////////
// AutoBSTR class (BSTR wrapper)
///////////////////////////////////

const AutoBSTR& AutoBSTR::operator=(const BSTR& BSTRsrc)
{
    if(m_clsBSTR != NULL) {
        SysFreeString(m_clsBSTR);
    }

    m_clsBSTR = SysAllocString(BSTRsrc);

	return *this;
}

const AutoBSTR& AutoBSTR::operator=(const char *textSrc)
{
    DWORD len;

    if(m_clsBSTR != NULL) {
        SysFreeString(m_clsBSTR);
    }

    // get the length (in wide chars) of the Unicode buffer needed to
    // convert 'textSrc' to Unicode.
    len = MultiByteToWideChar(CP_ACP, 0, textSrc, -1, 0, 0);
    // Allocate the Unicode buffer. SysAllocStringLen() will also allocate
    // the space for a terminating 0 short, and the unsigned long count.
    // SysAllocStringLen() will fill in the unsigned long with the value of
    // "len * sizeof(wchar_t)" and then return a pointer to the third short in
    // the buffer it allocates
    m_clsBSTR = SysAllocStringLen(0, len);
    // Convert 'textSrc' to Unicode in the allocated buffer
    MultiByteToWideChar(CP_ACP, 0, textSrc, -1, m_clsBSTR, len);

	return *this;
}


/////////////////////////
// AutoVariant class
/////////////////////////

const AutoVariant& AutoVariant::operator=(const VARIANT& varSrc)
{
    VariantCopy(this, (LPVARIANT)&varSrc);

	return *this;
}

const AutoVariant& AutoVariant::operator=(const AutoVariant& varSrc)
{
	VariantCopy(this, (LPVARIANT)&varSrc);

	return *this;
}

const AutoVariant& AutoVariant::operator=(const char* pCharSrc)
{
    DWORD len;

	// free the base VARIANT
	Clear();

	vt = VT_BSTR;
	// if the input C string is empty, store the same
	if (pCharSrc == NULL) {
		bstrVal = NULL;
	}
	else
	{
        // get the length (in wide chars) of the Unicode buffer we'll need to
        // convert pCharSrc to Unicode.
        len = MultiByteToWideChar(CP_ACP, 0, pCharSrc, -1, 0, 0);
        // allocate a Unicode buffer with the required length. SysAllocStringLen() will
        // also allocate the space for a terminating 0 short, and the unsigned long count.
        // SysAllocStringLen() will fill in the unsigned long with the value of
        // "len * sizeof(wchar_t)" and then will return a pointer to the third short in
        // the buffer it allocates.
        bstrVal = SysAllocStringLen(0, len);
        // convert pCharSrc to Unicode and put it in the buffer allocated by SysAllocStringLen()
        // remark: the assumption here is that SysFreeString() will be called by the AutoVariant
        // class destructor, via VariantClear()
        MultiByteToWideChar(CP_ACP, 0, pCharSrc, -1, bstrVal, len);
	}
	return *this;
}

const AutoVariant& AutoVariant::operator=(double dblSrc)
{
	// free the previous VARIANT, if of a different type
	if (vt != VT_R8)
	{
		Clear();
		vt = VT_R8;
	}

	dblVal = dblSrc;

	return *this;
}


/////////////////////////
// AutoSafeArray class
/////////////////////////

void AutoSafeArray::Create(VARTYPE vtSrc, DWORD dwDims, DWORD* rgElements)
{
    DWORD dwIndex;
    SAFEARRAYBOUND* rgsaBounds;

	ASSERT(rgElements != NULL);

	// Allocate and fill proxy array of bounds (with lower bound of zero)
	rgsaBounds = new SAFEARRAYBOUND[dwDims];
	ASSERT(rgsaBounds != NULL);

	for (dwIndex = 0; dwIndex < dwDims; dwIndex++)
	{
		// Assume lower bound is 0 and fill in element count
		rgsaBounds[dwIndex].lLbound = 0;
		rgsaBounds[dwIndex].cElements = rgElements[dwIndex];
	}

	try	{
		Create(vtSrc, dwDims, rgsaBounds);
	}
    catch( bad_alloc& )	{
		// Must free up memory
		delete[] rgsaBounds;
		rgsaBounds = NULL;
	}
    catch(...)	{
		// Must free up memory
		delete[] rgsaBounds;
		rgsaBounds = NULL;
	}

    if(rgsaBounds != NULL) {
        delete[] rgsaBounds;
    }
}

void AutoSafeArray::Create(VARTYPE vtSrc, DWORD dwDims, SAFEARRAYBOUND* rgsabound)
{
	ASSERT(dwDims > 0);
	ASSERT(rgsabound != NULL);

	// Validate the VARTYPE for SafeArrayCreate call
	ASSERT(!(vtSrc & VT_ARRAY));
	ASSERT(!(vtSrc & VT_BYREF));
	ASSERT(!(vtSrc & VT_VECTOR));
	ASSERT(vtSrc != VT_EMPTY);
	ASSERT(vtSrc != VT_NULL);

	// Free up the AutoSafeArray
	Clear();

    // parray is a member of VARIANT base class union structure
	parray = SafeArrayCreate(vtSrc, dwDims, rgsabound);
	ASSERT(parray != NULL);

	vt = (unsigned short) (vtSrc | VT_ARRAY);
}

void AutoSafeArray::PutElement(long* rgIndices, void* pvData)
{
	SafeArrayPutElement(parray, rgIndices, pvData);
}

// helper function
void SafeArrayInit(AutoSafeArray* psa)
{
	memset(psa, 0, sizeof(*psa));
}
