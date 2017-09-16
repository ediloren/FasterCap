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


#ifndef ABOUTBOX_H
#define ABOUTBOX_H

#ifndef WX_PRECOMP
	//(*HeadersPCH(AboutBox)
	#include <wx/sizer.h>
	#include <wx/stattext.h>
	#include <wx/panel.h>
	#include <wx/statbmp.h>
	#include <wx/button.h>
	#include <wx/dialog.h>
	//*)
#endif
//(*Headers(AboutBox)
//*)

class AboutBox: public wxDialog
{
	public:

		AboutBox(wxWindow* parent,wxWindowID id=wxID_ANY,const wxPoint& pos=wxDefaultPosition,const wxSize& size=wxDefaultSize);
		virtual ~AboutBox();

		//(*Declarations(AboutBox)
		wxButton* Button_ShowLicense;
		wxStaticText* StaticText_Version;
		wxStaticBitmap* StaticBitmap1;
		wxPanel* Panel1;
		wxStaticText* StaticText_WebSite;
		wxStaticText* StaticText_Copyright;
		//*)

	protected:

		//(*Identifiers(AboutBox)
		static const long ID_STATICBITMAP1;
		static const long ID_ABOUTBOX_VERSION;
		static const long ID_ABOUTBOX_COPYRIGHT;
		static const long ID_ABOUTBOX_WEBSITE;
		static const long ID_ABOUTBOX_SHOWLICENSE;
		static const long ID_PANEL1;
		//*)

	private:

		//(*Handlers(AboutBox)
		void OnButton_ShowLicenseClick(wxCommandEvent& event);
		//*)

		DECLARE_EVENT_TABLE()
};

#endif
