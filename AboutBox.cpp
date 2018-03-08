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

#include "wx_pch.h"
#include "AboutBox.h"

#ifndef WX_PRECOMP
	//(*InternalHeadersPCH(AboutBox)
	#include <wx/intl.h>
	#include <wx/string.h>
	//*)
#endif
//(*InternalHeaders(AboutBox)
//*)

#include "FasterCapGlobal.h"

#include "LicenseBox.h"

// include the big application icon
#include "res/FasterCap_32x32.xpm"

//(*IdInit(AboutBox)
const long AboutBox::ID_STATICBITMAP1 = wxNewId();
const long AboutBox::ID_ABOUTBOX_VERSION = wxNewId();
const long AboutBox::ID_ABOUTBOX_COPYRIGHT = wxNewId();
const long AboutBox::ID_ABOUTBOX_WEBSITE = wxNewId();
const long AboutBox::ID_ABOUTBOX_SHOWLICENSE = wxNewId();
const long AboutBox::ID_PANEL1 = wxNewId();
//*)

BEGIN_EVENT_TABLE(AboutBox,wxDialog)
	//(*EventTable(AboutBox)
	//*)
END_EVENT_TABLE()

AboutBox::AboutBox(wxWindow* parent,wxWindowID ,const wxPoint& ,const wxSize& )
{
	//(*Initialize(AboutBox)
	wxBoxSizer* BoxSizer2;
	wxBoxSizer* BoxSizer1;
	wxFlexGridSizer* FlexGridSizer1;
	wxBoxSizer* BoxSizer3;
	wxStdDialogButtonSizer* StdDialogButtonSizer1;

	Create(parent, wxID_ANY, _("About"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("wxID_ANY"));
	BoxSizer1 = new wxBoxSizer(wxHORIZONTAL);
	Panel1 = new wxPanel(this, ID_PANEL1, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL, _T("ID_PANEL1"));
	FlexGridSizer1 = new wxFlexGridSizer(2, 2, 0, 0);
	StaticBitmap1 = new wxStaticBitmap(Panel1, ID_STATICBITMAP1, wxNullBitmap, wxDefaultPosition, wxDefaultSize, wxSIMPLE_BORDER, _T("ID_STATICBITMAP1"));
	FlexGridSizer1->Add(StaticBitmap1, 1, wxALL|wxALIGN_TOP|wxALIGN_CENTER_HORIZONTAL, 5);
	BoxSizer2 = new wxBoxSizer(wxVERTICAL);
	StaticText_Version = new wxStaticText(Panel1, ID_ABOUTBOX_VERSION, _("FasterCap Version X.X.X"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_ABOUTBOX_VERSION"));
	BoxSizer2->Add(StaticText_Version, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	StaticText_Copyright = new wxStaticText(Panel1, ID_ABOUTBOX_COPYRIGHT, _("Copyright (C) 2018 FastFieldSolvers S.R.L."), wxDefaultPosition, wxDefaultSize, 0, _T("ID_ABOUTBOX_COPYRIGHT"));
	BoxSizer2->Add(StaticText_Copyright, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	StaticText_WebSite = new wxStaticText(Panel1, ID_ABOUTBOX_WEBSITE, _("http://www.fastfieldsolvers.com"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_ABOUTBOX_WEBSITE"));
	BoxSizer2->Add(StaticText_WebSite, 1, wxALL|wxEXPAND|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer1->Add(BoxSizer2, 1, wxALL|wxEXPAND|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer1->Add(-1,-1,1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	BoxSizer3 = new wxBoxSizer(wxHORIZONTAL);
	Button_ShowLicense = new wxButton(Panel1, ID_ABOUTBOX_SHOWLICENSE, _("Show License"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_ABOUTBOX_SHOWLICENSE"));
	BoxSizer3->Add(Button_ShowLicense, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	StdDialogButtonSizer1 = new wxStdDialogButtonSizer();
	StdDialogButtonSizer1->AddButton(new wxButton(Panel1, wxID_OK, wxEmptyString));
	StdDialogButtonSizer1->Realize();
	BoxSizer3->Add(StdDialogButtonSizer1, 1, wxALL|wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer1->Add(BoxSizer3, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	Panel1->SetSizer(FlexGridSizer1);
	FlexGridSizer1->Fit(Panel1);
	FlexGridSizer1->SetSizeHints(Panel1);
	BoxSizer1->Add(Panel1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	SetSizer(BoxSizer1);
	BoxSizer1->Fit(this);
	BoxSizer1->SetSizeHints(this);

	Connect(ID_ABOUTBOX_SHOWLICENSE,wxEVT_COMMAND_BUTTON_CLICKED,(wxObjectEventFunction)&AboutBox::OnButton_ShowLicenseClick);
	//*)

    // update with correct global values
	StaticText_Version->SetLabel(wxT(FCG_HEADER_VERSION));
	StaticText_Copyright->SetLabel(wxT(FCG_HEADER_COPYRIGHT));
	StaticText_WebSite->SetLabel(wxT(FCG_HEADER_WEBSITE));
	StaticBitmap1->SetBitmap(wxBitmap(FasterCap_32x32_xpm));
    //
	// and rearrange the dimension, since the static texts / images etc. can be smaller or larger
	//
	FlexGridSizer1->Fit(Panel1);
	// can/must call SetSizeHints() because 'Panel1' is derived from wxWindow,
	// while the other object in the hierarchy (e.g. BoxSizer) etc. derive from wxObject only
	FlexGridSizer1->SetSizeHints(Panel1);
	BoxSizer1->Fit(this);
	BoxSizer1->SetSizeHints(this);
}

AboutBox::~AboutBox()
{
	//(*Destroy(AboutBox)
	//*)
}

void AboutBox::OnButton_ShowLicenseClick(wxCommandEvent& )
{
    LicenseBox licenseDlg(this);
    licenseDlg.ShowModal();
}
