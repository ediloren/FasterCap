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


#ifndef LICENSEBOX_H
#define LICENSEBOX_H

#ifndef WX_PRECOMP
	//(*HeadersPCH(LicenseBox)
	#include <wx/sizer.h>
	#include <wx/stattext.h>
	#include <wx/textctrl.h>
	#include <wx/panel.h>
	#include <wx/dialog.h>
	//*)
#endif
//(*Headers(LicenseBox)
//*)

class LicenseBox: public wxDialog
{
	public:

		LicenseBox(wxWindow* parent,wxWindowID id=wxID_ANY);
		virtual ~LicenseBox();

		//(*Declarations(LicenseBox)
		wxStaticText* StaticText_Header2;
		wxStaticText* StaticText_Header3;
		wxPanel* Panel1;
		wxStaticText* StaticText_Header1;
		wxTextCtrl* TextCtrl_License;
		//*)

	protected:

		//(*Identifiers(LicenseBox)
		static const long ID_STATICTEXT_HEADER1;
		static const long ID_STATICTEXT_HEADER2;
		static const long ID_STATICTEXT_HEADER3;
		static const long ID_TEXTCTRL_LICENSE;
		static const long ID_PANEL1;
		//*)

	private:

		//(*Handlers(LicenseBox)
		//*)

		DECLARE_EVENT_TABLE()
};

#endif
