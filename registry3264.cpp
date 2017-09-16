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

//   class derived from wxWidgets wxRegKey to manage 64 bits vs. 32 bits keys

#include "wx_pch.h"

// valid only for MS Win
#ifdef __WXMSW__

#include "registry3264.h"

wxRegKey3264::wxRegKey3264(const wxString& strKey, REGSAM mode)
{
  key3264access = mode;
  ::wxRegKey(strKey);
}

// parent is a normal regkey
wxRegKey3264::wxRegKey3264(const wxRegKey3264& keyParent, const wxString& strKey, REGSAM mode)
{
  key3264access = mode;
  ::wxRegKey(keyParent, strKey);
}

void wxRegKey3264::SetKey32bits()
{
  key3264access = KEY_WOW64_32KEY;
}

// opens key (it's not an error to call Open() on an already opened key)
bool wxRegKey3264::Open(AccessMode mode)
{
    REGSAM mode3264;

    if ( IsOpened() )
    {
        if ( mode <= m_mode )
            return true;

        // we had been opened in read mode but now must be reopened in write
        Close();
    }

    HKEY tmpKey;

    mode3264 = (Read ? KEY_READ : KEY_ALL_ACCESS) | key3264access;

    m_dwLastError = ::RegOpenKeyEx
                    (
                        (HKEY) m_hRootKey,
                        m_strKey,
                        RESERVED,
                        mode3264,
                        &tmpKey
                    );

    if ( m_dwLastError != ERROR_SUCCESS )
    {
        wxLogSysError(m_dwLastError, _("Can't open registry key '%s'"),
                      (const char*)GetName());
        return false;
    }

    m_hKey = (WXHKEY) tmpKey;
    m_mode = mode;

    return true;
}

// returns true if the key exists
bool wxRegKey3264::Exists() const
{
  // opened key has to exist, try to open it if not done yet
  return IsOpened() ? true : KeyExists(m_hRootKey, m_strKey);
}

// returns true if given subkey exists
bool wxRegKey3264::HasSubKey(const wxChar *szKey) const
{
  // this function should be silent, so suppress possible messages from Open()
  wxLogNull nolog;

  if ( !CONST_CAST Open(Read) )
    return false;

  return KeyExists(m_hKey, szKey);
}

// creates key, failing if it exists and !bOkIfExists
bool wxRegKey3264::Create(bool bOkIfExists)
{
  // check for existence only if asked (i.e. order is important!)
  if ( !bOkIfExists && Exists() )
    return false;

  if ( IsOpened() )
    return true;

  HKEY tmpKey;

  DWORD disposition;
  m_dwLastError = RegCreateKeyEx((HKEY) m_hRootKey, m_strKey,
      NULL, // reserved
      NULL, // class string
      0,
      KEY_READ | KEY_WRITE | key3264access,
      NULL,
      &tmpKey,
      &disposition);

  if ( m_dwLastError != ERROR_SUCCESS ) {
    wxLogSysError(m_dwLastError, _("Can't create registry key '%s'"),
                  (const char*)GetName());
    return false;
  }
  else
  {
    m_hKey = (WXHKEY) tmpKey;
    return true;
  }
}

bool wxRegKey3264::KeyExists(WXHKEY hRootKey, const wxChar *szKey)
{
    // don't close this key itself for the case of empty szKey!
    if ( wxIsEmpty(szKey) )
        return true;

    HKEY hkeyDummy;
    if ( ::RegOpenKeyEx
         (
            (HKEY)hRootKey,
            szKey,
            RESERVED,
            KEY_READ | key3264access,  // we might not have enough rights for rw access
            &hkeyDummy
         ) == ERROR_SUCCESS )
    {
        ::RegCloseKey(hkeyDummy);

        return true;
    }

    return false;
}

#endif
