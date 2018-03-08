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
#include "LicenseBox.h"

#ifndef WX_PRECOMP
	//(*InternalHeadersPCH(LicenseBox)
	#include <wx/intl.h>
	#include <wx/string.h>
	//*)
#endif

#include <wx/textfile.h>
#include <wx/filename.h>

// for #defines
#include "FasterCapGlobal.h"

//(*InternalHeaders(LicenseBox)
#include <wx/button.h>
//*)

//(*IdInit(LicenseBox)
const long LicenseBox::ID_STATICTEXT_HEADER1 = wxNewId();
const long LicenseBox::ID_STATICTEXT_HEADER2 = wxNewId();
const long LicenseBox::ID_STATICTEXT_HEADER3 = wxNewId();
const long LicenseBox::ID_TEXTCTRL_LICENSE = wxNewId();
const long LicenseBox::ID_PANEL1 = wxNewId();
//*)

BEGIN_EVENT_TABLE(LicenseBox,wxDialog)
	//(*EventTable(LicenseBox)
	//*)
END_EVENT_TABLE()

LicenseBox::LicenseBox(wxWindow* parent,wxWindowID )
{
    wxTextFile file;
    wxFileName fileName;
    unsigned int i;

	//(*Initialize(LicenseBox)
	wxBoxSizer* BoxSizer2;
	wxBoxSizer* BoxSizer1;
	wxFlexGridSizer* FlexGridSizer1;
	wxStdDialogButtonSizer* StdDialogButtonSizer1;

	Create(parent, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("wxID_ANY"));
	BoxSizer1 = new wxBoxSizer(wxHORIZONTAL);
	Panel1 = new wxPanel(this, ID_PANEL1, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL, _T("ID_PANEL1"));
	FlexGridSizer1 = new wxFlexGridSizer(3, 1, 0, 0);
	BoxSizer2 = new wxBoxSizer(wxVERTICAL);
	StaticText_Header1 = new wxStaticText(Panel1, ID_STATICTEXT_HEADER1, _("FasterCap Demo License"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT_HEADER1"));
	BoxSizer2->Add(StaticText_Header1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticText_Header2 = new wxStaticText(Panel1, ID_STATICTEXT_HEADER2, _("Please contact FastFieldSolvers S.r.l. if you wish to purchase a commercial license"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT_HEADER2"));
	BoxSizer2->Add(StaticText_Header2, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticText_Header3 = new wxStaticText(Panel1, ID_STATICTEXT_HEADER3, _("http://www.fastfieldsolvers.com"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT_HEADER3"));
	BoxSizer2->Add(StaticText_Header3, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer1->Add(BoxSizer2, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	TextCtrl_License = new wxTextCtrl(Panel1, ID_TEXTCTRL_LICENSE, _("License"), wxDefaultPosition, wxSize(400,300), wxTE_AUTO_SCROLL|wxTE_MULTILINE|wxTE_WORDWRAP, wxDefaultValidator, _T("ID_TEXTCTRL_LICENSE"));
	FlexGridSizer1->Add(TextCtrl_License, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StdDialogButtonSizer1 = new wxStdDialogButtonSizer();
	StdDialogButtonSizer1->AddButton(new wxButton(Panel1, wxID_OK, _("I accept")));
	StdDialogButtonSizer1->Realize();
	FlexGridSizer1->Add(StdDialogButtonSizer1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	Panel1->SetSizer(FlexGridSizer1);
	FlexGridSizer1->Fit(Panel1);
	FlexGridSizer1->SetSizeHints(Panel1);
	BoxSizer1->Add(Panel1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	SetSizer(BoxSizer1);
	BoxSizer1->Fit(this);
	BoxSizer1->SetSizeHints(this);
	//*)

    StaticText_Header1->SetLabel(wxT(FCG_LIC_STD_HEADER1));
    StaticText_Header2->SetLabel(wxT(FCG_HEADER_COPYRIGHT));
    StaticText_Header3->SetLabel(wxT(FCG_HEADER_WEBSITE));

	// load license text

    TextCtrl_License->Clear();

    if ( file.Open((Globals::GetApp())->GetLicenseTextPath()) ) {
        for (i = 0; i < file.GetLineCount(); i++) {
            TextCtrl_License->AppendText(file[i]);
            TextCtrl_License->AppendText(wxT("\n"));
        }
        file.Close();
    }
    else {
        TextCtrl_License->SetDefaultStyle(wxTextAttr(*wxRED));
        TextCtrl_License->AppendText("ERROR - License text cannot be read\n");
        TextCtrl_License->AppendText("You are NOT authorized to use the software\n");
    }

    //
	// and rearrange the dimension, since the static texts / images etc. can be smaller or larger
	//
	FlexGridSizer1->Fit(Panel1);
	// can/must call SetSizeHints() because 'Panel1' is derived from wxWindow,
	// while the other object in the hierarchy (e.g. BoxSizer) etc. derive from wxObject only
	FlexGridSizer1->SetSizeHints(Panel1);
	BoxSizer1->Fit(this);
	BoxSizer1->SetSizeHints(this);

    // go to beginning of text, first by unselecting all the text, then moving to first pos
    // remark: must use SetFocus() before, or it does not unselect the text (selecting when
    // getting focus is standard wxWidgets behavior, must explicitly unselect)
    TextCtrl_License->SetFocus();
    TextCtrl_License->ShowPosition(0);
    TextCtrl_License->SetSelection(0, 0);
}

LicenseBox::~LicenseBox()
{
	//(*Destroy(LicenseBox)
	//*)
}

