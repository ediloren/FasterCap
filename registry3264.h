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

 // Purpose:   class derived from wxWidgets wxRegKey to manage 64 bits vs. 32 bits keys

#include "wx_pch.h"

// valid only for MS Win
#ifdef __WXMSW__

#include <wx/msw/registry.h>

#define KEY_REGKEY3264_DEFAULT  (0)

class wxRegKey3264 : public wxRegKey
{
public:
  wxRegKey3264(const wxString& strKey, REGSAM mode = KEY_REGKEY3264_DEFAULT);
  wxRegKey3264(const wxRegKey3264& keyParent, const wxString& strKey, REGSAM mode = KEY_REGKEY3264_DEFAULT);
  wxRegKey3264() : key3264access(KEY_REGKEY3264_DEFAULT) {}
  bool  Open(AccessMode mode = Write);
  bool  Create(bool bOkIfExists = true);
  bool  Exists() const;
  bool  HasSubKey(const wxChar *szKey) const;
  void  SetKey32bits();

protected:
  bool  KeyExists(WXHKEY hRootKey, const wxChar *szKey);

  REGSAM  key3264access;
};

#endif //__WXMSW__
