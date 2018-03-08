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
#include "RunDialog.h"

#ifndef WX_PRECOMP
	//(*InternalHeadersPCH(RunDialog)
	#include <wx/intl.h>
	#include <wx/string.h>
	//*)
#endif
//(*InternalHeaders(RunDialog)
//*)

//(*IdInit(RunDialog)
const long RunDialog::ID_STATICTEXT1 = wxNewId();
const long RunDialog::ID_INPUT_FILE_NAME = wxNewId();
const long RunDialog::ID_BROWSE_INPUT_FILE = wxNewId();
const long RunDialog::ID_AUTO_SETTINGS = wxNewId();
const long RunDialog::ID_STATICTEXT2 = wxNewId();
const long RunDialog::ID_AUTOMATIC_ERROR = wxNewId();
const long RunDialog::ID_AUTO_PRECOND = wxNewId();
const long RunDialog::ID_STATICTEXT7 = wxNewId();
const long RunDialog::ID_MESH_EPS = wxNewId();
const long RunDialog::ID_STATICTEXT8 = wxNewId();
const long RunDialog::ID_GMRES_TOL = wxNewId();
const long RunDialog::ID_STATICTEXT3 = wxNewId();
const long RunDialog::ID_EPS_RATIO = wxNewId();
const long RunDialog::ID_STATICTEXT4 = wxNewId();
const long RunDialog::ID_CURVATURE = wxNewId();
const long RunDialog::ID_STATICTEXT5 = wxNewId();
const long RunDialog::ID_AUTOREFINE_OUT_OF_CORE = wxNewId();
const long RunDialog::ID_GALERKIN_SCHEME = wxNewId();
const long RunDialog::ID_PRECOND_JACOBI = wxNewId();
const long RunDialog::ID_PRECOND_TWO_LEVELS = wxNewId();
const long RunDialog::ID_STATICTEXT6 = wxNewId();
const long RunDialog::ID_DIM_PRECOND_TWO_LEVELS = wxNewId();
const long RunDialog::ID_OUTPUT_GEOMETRY = wxNewId();
const long RunDialog::ID_DUMP_GEOMETRY = wxNewId();
const long RunDialog::ID_OUTPUT_CHARGE = wxNewId();
const long RunDialog::ID_DUMP_RESIDUAL = wxNewId();
const long RunDialog::ID_DUMP_TIME_MEM = wxNewId();
const long RunDialog::ID_VERBOSE_OUTPUT = wxNewId();
const long RunDialog::ID_OUTPUT_CAPACITANCE = wxNewId();
const long RunDialog::ID_PANEL1 = wxNewId();
//*)

BEGIN_EVENT_TABLE(RunDialog,wxDialog)
	//(*EventTable(RunDialog)
	//*)
END_EVENT_TABLE()

RunDialog::RunDialog(wxWindow* parent,wxWindowID )
{
	//(*Initialize(RunDialog)
	wxStaticBoxSizer* StaticBoxSizer2;
	wxFlexGridSizer* FlexGridSizer4;
	wxStaticBoxSizer* StaticBoxSizer4;
	wxFlexGridSizer* FlexGridSizer3;
	wxFlexGridSizer* FlexGridSizer5;
	wxFlexGridSizer* FlexGridSizer9;
	wxFlexGridSizer* FlexGridSizer2;
	wxFlexGridSizer* FlexGridSizer7;
	wxStaticBoxSizer* StaticBoxSizer3;
	wxGridSizer* GridSizer1;
	wxFlexGridSizer* FlexGridSizer8;
	wxBoxSizer* BoxSizer1;
	wxFlexGridSizer* FlexGridSizer6;
	wxStaticBoxSizer* StaticBoxSizer1;
	wxFlexGridSizer* FlexGridSizer1;
	wxStaticBoxSizer* StaticBoxSizer5;
	wxStdDialogButtonSizer* StdDialogButtonSizer1;

	Create(parent, wxID_ANY, _("Run Menu"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE, _T("wxID_ANY"));
	BoxSizer1 = new wxBoxSizer(wxHORIZONTAL);
	Panel1 = new wxPanel(this, ID_PANEL1, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL, _T("ID_PANEL1"));
	FlexGridSizer1 = new wxFlexGridSizer(0, 1, 0, 0);
	FlexGridSizer2 = new wxFlexGridSizer(0, 3, 0, 0);
	StaticText1 = new wxStaticText(Panel1, ID_STATICTEXT1, _("Input file name"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT1"));
	FlexGridSizer2->Add(StaticText1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	TextCtrl_InputFileName = new wxTextCtrl(Panel1, ID_INPUT_FILE_NAME, wxEmptyString, wxDefaultPosition, wxSize(360,21), 0, wxDefaultValidator, _T("ID_INPUT_FILE_NAME"));
	FlexGridSizer2->Add(TextCtrl_InputFileName, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	Button_BrowseInputFile = new wxButton(Panel1, ID_BROWSE_INPUT_FILE, _("Browse"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_BROWSE_INPUT_FILE"));
	FlexGridSizer2->Add(Button_BrowseInputFile, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer1->Add(FlexGridSizer2, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticBoxSizer1 = new wxStaticBoxSizer(wxHORIZONTAL, Panel1, _("Automatic Settings"));
	FlexGridSizer3 = new wxFlexGridSizer(0, 3, 0, 0);
	CheckBox_AutoSettings = new wxCheckBox(Panel1, ID_AUTO_SETTINGS, _("Automatically calculate settings"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_AUTO_SETTINGS"));
	CheckBox_AutoSettings->SetValue(false);
	FlexGridSizer3->Add(CheckBox_AutoSettings, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	StaticText2 = new wxStaticText(Panel1, ID_STATICTEXT2, _("Stop when relative error is lower than (-a)"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT2"));
	FlexGridSizer3->Add(StaticText2, 1, wxLEFT|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 40);
	TextCtrl_Automatic_Error = new wxTextCtrl(Panel1, ID_AUTOMATIC_ERROR, wxEmptyString, wxDefaultPosition, wxSize(50,21), 0, wxDefaultValidator, _T("ID_AUTOMATIC_ERROR"));
	FlexGridSizer3->Add(TextCtrl_Automatic_Error, 1, wxALL|wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_Auto_Precond = new wxCheckBox(Panel1, ID_AUTO_PRECOND, _("Automatic preconditioner usage (-ap)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_AUTO_PRECOND"));
	CheckBox_Auto_Precond->SetValue(false);
	FlexGridSizer3->Add(CheckBox_Auto_Precond, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	StaticBoxSizer1->Add(FlexGridSizer3, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	FlexGridSizer1->Add(StaticBoxSizer1, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticBoxSizer5 = new wxStaticBoxSizer(wxHORIZONTAL, Panel1, _("Manual Settings"));
	FlexGridSizer9 = new wxFlexGridSizer(0, 4, 0, 0);
	StaticText7 = new wxStaticText(Panel1, ID_STATICTEXT7, _("Mesh relative refinement value (-m)"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT7"));
	FlexGridSizer9->Add(StaticText7, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	TextCtrl_MeshEps = new wxTextCtrl(Panel1, ID_MESH_EPS, wxEmptyString, wxDefaultPosition, wxSize(50,21), 0, wxDefaultValidator, _T("ID_MESH_EPS"));
	FlexGridSizer9->Add(TextCtrl_MeshEps, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticText8 = new wxStaticText(Panel1, ID_STATICTEXT8, _("Gmres iteration tolerance (-t)"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT8"));
	FlexGridSizer9->Add(StaticText8, 1, wxLEFT|wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL, 60);
	TextCtrl_GmresTol = new wxTextCtrl(Panel1, ID_GMRES_TOL, wxEmptyString, wxDefaultPosition, wxSize(58,21), 0, wxDefaultValidator, _T("ID_GMRES_TOL"));
	FlexGridSizer9->Add(TextCtrl_GmresTol, 1, wxALL|wxALIGN_RIGHT|wxALIGN_CENTER_VERTICAL, 5);
	StaticBoxSizer5->Add(FlexGridSizer9, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	FlexGridSizer1->Add(StaticBoxSizer5, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticBoxSizer2 = new wxStaticBoxSizer(wxHORIZONTAL, Panel1, _("Common Settings"));
	FlexGridSizer4 = new wxFlexGridSizer(0, 1, 0, 0);
	FlexGridSizer5 = new wxFlexGridSizer(0, 4, 0, 0);
	StaticText3 = new wxStaticText(Panel1, ID_STATICTEXT3, _("Direct potential interaction coefficient to\nmesh refinement ratio (-d)"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT3"));
	FlexGridSizer5->Add(StaticText3, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	TextCtrl_EPS_Ratio = new wxTextCtrl(Panel1, ID_EPS_RATIO, wxEmptyString, wxDefaultPosition, wxSize(50,21), 0, wxDefaultValidator, _T("ID_EPS_RATIO"));
	FlexGridSizer5->Add(TextCtrl_EPS_Ratio, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticText4 = new wxStaticText(Panel1, ID_STATICTEXT4, _("Mesh curvature coefficient (-mc)"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT4"));
	FlexGridSizer5->Add(StaticText4, 1, wxLEFT|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 20);
	TextCtrl_Curvature = new wxTextCtrl(Panel1, ID_CURVATURE, wxEmptyString, wxDefaultPosition, wxSize(50,21), 0, wxDefaultValidator, _T("ID_CURVATURE"));
	FlexGridSizer5->Add(TextCtrl_Curvature, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer4->Add(FlexGridSizer5, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	FlexGridSizer6 = new wxFlexGridSizer(0, 3, 0, 0);
	StaticText5 = new wxStaticText(Panel1, ID_STATICTEXT5, _("Out-Of-Core free memory to link memory \ncondition (0 = don\'t go OOC) (-f)"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT5"));
	FlexGridSizer6->Add(StaticText5, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	TextCtrl_Autorefine_Out_of_Core = new wxTextCtrl(Panel1, ID_AUTOREFINE_OUT_OF_CORE, wxEmptyString, wxDefaultPosition, wxSize(50,21), 0, wxDefaultValidator, _T("ID_AUTOREFINE_OUT_OF_CORE"));
	FlexGridSizer6->Add(TextCtrl_Autorefine_Out_of_Core, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_GalerkinScheme = new wxCheckBox(Panel1, ID_GALERKIN_SCHEME, _("Use Galerkin scheme (-g)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_GALERKIN_SCHEME"));
	CheckBox_GalerkinScheme->SetValue(false);
	FlexGridSizer6->Add(CheckBox_GalerkinScheme, 1, wxLEFT|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 50);
	FlexGridSizer4->Add(FlexGridSizer6, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	StaticBoxSizer3 = new wxStaticBoxSizer(wxHORIZONTAL, Panel1, _("Preconditioner Settings"));
	FlexGridSizer7 = new wxFlexGridSizer(0, 1, 0, 0);
	CheckBox_Precond_Jacobi = new wxCheckBox(Panel1, ID_PRECOND_JACOBI, _("Use Jacobi Preconditioner (-pj)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_PRECOND_JACOBI"));
	CheckBox_Precond_Jacobi->SetValue(false);
	FlexGridSizer7->Add(CheckBox_Precond_Jacobi, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer8 = new wxFlexGridSizer(0, 3, 0, 0);
	CheckBox_Precond_TwoLevels = new wxCheckBox(Panel1, ID_PRECOND_TWO_LEVELS, _("Use Two-levels Preconditioner"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_PRECOND_TWO_LEVELS"));
	CheckBox_Precond_TwoLevels->SetValue(false);
	FlexGridSizer8->Add(CheckBox_Precond_TwoLevels, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticText6 = new wxStaticText(Panel1, ID_STATICTEXT6, _("Two-levels Preconditioner dimension (-ps)"), wxDefaultPosition, wxDefaultSize, 0, _T("ID_STATICTEXT6"));
	FlexGridSizer8->Add(StaticText6, 1, wxLEFT|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 30);
	TextCtrl_Dim_Precond_Two_Levels = new wxTextCtrl(Panel1, ID_DIM_PRECOND_TWO_LEVELS, wxEmptyString, wxDefaultPosition, wxSize(50,21), 0, wxDefaultValidator, _T("ID_DIM_PRECOND_TWO_LEVELS"));
	FlexGridSizer8->Add(TextCtrl_Dim_Precond_Two_Levels, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	FlexGridSizer7->Add(FlexGridSizer8, 1, wxTOP|wxBOTTOM|wxRIGHT|wxEXPAND|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 0);
	StaticBoxSizer3->Add(FlexGridSizer7, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	FlexGridSizer4->Add(StaticBoxSizer3, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	StaticBoxSizer2->Add(FlexGridSizer4, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	FlexGridSizer1->Add(StaticBoxSizer2, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StaticBoxSizer4 = new wxStaticBoxSizer(wxHORIZONTAL, Panel1, _("General Options"));
	GridSizer1 = new wxGridSizer(0, 2, 0, 0);
	CheckBox_Output_Geometry = new wxCheckBox(Panel1, ID_OUTPUT_GEOMETRY, _("Output refined geometry in FastCap2 format (-o)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_OUTPUT_GEOMETRY"));
	CheckBox_Output_Geometry->SetValue(false);
	GridSizer1->Add(CheckBox_Output_Geometry, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_Dump_Geometry = new wxCheckBox(Panel1, ID_DUMP_GEOMETRY, _("Dump input geometry in FasterCap format and stop (-oi)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_DUMP_GEOMETRY"));
	CheckBox_Dump_Geometry->SetValue(false);
	GridSizer1->Add(CheckBox_Dump_Geometry, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_Output_Charge = new wxCheckBox(Panel1, ID_OUTPUT_CHARGE, _("Dump charge densities in output file (-c)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_OUTPUT_CHARGE"));
	CheckBox_Output_Charge->SetValue(false);
	GridSizer1->Add(CheckBox_Output_Charge, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_Dump_Residual = new wxCheckBox(Panel1, ID_DUMP_RESIDUAL, _("Dump Gmres residual at each iteration (-r)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_DUMP_RESIDUAL"));
	CheckBox_Dump_Residual->SetValue(false);
	GridSizer1->Add(CheckBox_Dump_Residual, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_Dump_Time_Mem = new wxCheckBox(Panel1, ID_DUMP_TIME_MEM, _("Dump detailed time and memory information (-i)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_DUMP_TIME_MEM"));
	CheckBox_Dump_Time_Mem->SetValue(false);
	GridSizer1->Add(CheckBox_Dump_Time_Mem, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_Verbose_Output = new wxCheckBox(Panel1, ID_VERBOSE_OUTPUT, _("Verbose output (-v)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_VERBOSE_OUTPUT"));
	CheckBox_Verbose_Output->SetValue(false);
	GridSizer1->Add(CheckBox_Verbose_Output, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	CheckBox_Output_Capacitance = new wxCheckBox(Panel1, ID_OUTPUT_CAPACITANCE, _("Output capacitance matrix to file (-e)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("ID_OUTPUT_CAPACITANCE"));
	CheckBox_Output_Capacitance->SetValue(false);
	GridSizer1->Add(CheckBox_Output_Capacitance, 1, wxALL|wxALIGN_LEFT|wxALIGN_CENTER_VERTICAL, 5);
	StaticBoxSizer4->Add(GridSizer1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	FlexGridSizer1->Add(StaticBoxSizer4, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	StdDialogButtonSizer1 = new wxStdDialogButtonSizer();
	StdDialogButtonSizer1->AddButton(new wxButton(Panel1, wxID_OK, _("Run")));
	StdDialogButtonSizer1->AddButton(new wxButton(Panel1, wxID_NO, _("Reset")));
	StdDialogButtonSizer1->AddButton(new wxButton(Panel1, wxID_CANCEL, wxEmptyString));
	StdDialogButtonSizer1->Realize();
	FlexGridSizer1->Add(StdDialogButtonSizer1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	Panel1->SetSizer(FlexGridSizer1);
	FlexGridSizer1->Fit(Panel1);
	FlexGridSizer1->SetSizeHints(Panel1);
	BoxSizer1->Add(Panel1, 1, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5);
	SetSizer(BoxSizer1);
	BoxSizer1->Fit(this);
	BoxSizer1->SetSizeHints(this);

	Connect(ID_BROWSE_INPUT_FILE,wxEVT_COMMAND_BUTTON_CLICKED,(wxObjectEventFunction)&RunDialog::OnButton_BrowseInputFileClick);
	//*)

	Connect(wxID_NO,wxEVT_COMMAND_BUTTON_CLICKED,(wxObjectEventFunction)&RunDialog::OnReset);

    //
    // link dialog fields with variables
    // using specific validators will work in next 2.9+ releases of wxWidgets
    //
    //wxFloatingPointValidator<double> myvalidator(&m_clsGlobalVars.m_dEpsRatio, wxNUM_VAL_NO_TRAILING_ZEROES);
	//TextCtrl_EPS_Ratio->SetValidator(myvalidator);

    // reset global vars
    Reset();
    // and copy them in the dialog
	CopyVarsToDialog();
}

RunDialog::~RunDialog()
{
	//(*Destroy(RunDialog)
	//*)
}

void RunDialog::OnButton_BrowseInputFileClick(wxCommandEvent& )
{

    wxFileDialog browseFileDialog(this, wxT("Open"), wxT(""), m_clsGlobalVars.m_sFileIn, wxT("FastCap List or Geometry files (*.lst, *.qui, *.txt)|*.lst;*.qui;*.txt|All files (*.*)|*.*"), wxFD_OPEN|wxFD_FILE_MUST_EXIST);

    if (browseFileDialog.ShowModal() != wxID_CANCEL) {
        TextCtrl_InputFileName->SetValue(browseFileDialog.GetPath());
        m_clsGlobalVars.m_sFileIn = TextCtrl_InputFileName->GetValue();
    }
}

// reset dialog box controls to defaults, except input file name
void RunDialog::OnReset(wxCommandEvent&)
{
    // reset the global vars
	Reset();
    // but keep the previous dialog file name
	m_clsGlobalVars.m_sFileIn = TextCtrl_InputFileName->GetValue();

    CopyVarsToDialog();
}

// reset dialog box controls to defaults
void RunDialog::Reset()
{
	m_clsGlobalVars.Reset();
	// modify defaults for run dialog defaults
	m_clsGlobalVars.m_bAuto = true;
	m_clsGlobalVars.m_bAutoPrecond = true;
	m_clsGlobalVars.m_ucPrecondType = AUTOREFINE_PRECOND_JACOBI;
}

void RunDialog::CopyVarsToDialog()
{

    TextCtrl_InputFileName->SetValue(m_clsGlobalVars.m_sFileIn);

    CheckBox_AutoSettings->SetValue(m_clsGlobalVars.m_bAuto);
    CheckBox_Auto_Precond->SetValue(m_clsGlobalVars.m_bAutoPrecond);
    TextCtrl_Automatic_Error->SetValue(wxString::Format("%g", m_clsGlobalVars.m_dAutoMaxErr));
    TextCtrl_Curvature->SetValue(wxString::Format("%g", m_clsGlobalVars.m_dMeshCurvCoeff));
    TextCtrl_MeshEps->SetValue(wxString::Format("%g", m_clsGlobalVars.m_dMeshEps));
    TextCtrl_EPS_Ratio->SetValue(wxString::Format("%g", m_clsGlobalVars.m_dEpsRatio));
    TextCtrl_Autorefine_Out_of_Core->SetValue(wxString::Format("%g", m_clsGlobalVars.m_dOutOfCoreRatio));
    TextCtrl_GmresTol->SetValue(wxString::Format("%g", m_clsGlobalVars.m_dGmresTol));
	if(m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_JACOBI) {
		CheckBox_Precond_Jacobi->SetValue(true);
	}
	else {
		CheckBox_Precond_Jacobi->SetValue(false);
	}
	if(m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) {
		CheckBox_Precond_TwoLevels->SetValue(true);
	}
	else {
		CheckBox_Precond_TwoLevels->SetValue(false);
	}
	TextCtrl_Dim_Precond_Two_Levels->SetValue(wxString::Format("%u",  m_clsGlobalVars.m_uiSuperPreDim));
	CheckBox_Output_Charge->SetValue(m_clsGlobalVars.m_bOutputCharge);
    CheckBox_Output_Capacitance->SetValue(m_clsGlobalVars.m_bOutputCapMtx);
	CheckBox_Dump_Residual->SetValue(m_clsGlobalVars.m_bDumpResidual);
	CheckBox_Dump_Time_Mem->SetValue(m_clsGlobalVars.m_bDumpTimeMem);
	CheckBox_Verbose_Output->SetValue(m_clsGlobalVars.m_bVerboseOutput);
	if(m_clsGlobalVars.m_cScheme == AUTOREFINE_COLLOCATION) {
		CheckBox_GalerkinScheme->SetValue(false);
    }
	else {
		CheckBox_GalerkinScheme->SetValue(true);
	}
	CheckBox_Output_Geometry->SetValue(m_clsGlobalVars.m_bOutputGeo);
	CheckBox_Dump_Geometry->SetValue(m_clsGlobalVars.m_bDumpInputGeo);
}

void RunDialog::CopyDialogToVars()
{
    unsigned long tmpULong;

    m_clsGlobalVars.m_sFileIn = TextCtrl_InputFileName->GetValue();

    m_clsGlobalVars.m_bAuto = CheckBox_AutoSettings->GetValue();
    m_clsGlobalVars.m_bAutoPrecond = CheckBox_Auto_Precond->GetValue();
    TextCtrl_Automatic_Error->GetValue().ToDouble(&m_clsGlobalVars.m_dAutoMaxErr);
    TextCtrl_Curvature->GetValue().ToDouble(&m_clsGlobalVars.m_dMeshCurvCoeff);
    TextCtrl_MeshEps->GetValue().ToDouble(&m_clsGlobalVars.m_dMeshEps);
    TextCtrl_EPS_Ratio->GetValue().ToDouble(&m_clsGlobalVars.m_dEpsRatio);
    TextCtrl_Autorefine_Out_of_Core->GetValue().ToDouble(&m_clsGlobalVars.m_dOutOfCoreRatio);
    TextCtrl_GmresTol->GetValue().ToDouble(&m_clsGlobalVars.m_dGmresTol);
	if(CheckBox_Precond_Jacobi->GetValue() == true) {
		m_clsGlobalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_JACOBI;
	}
	else {
		m_clsGlobalVars.m_ucPrecondType &= (~AUTOREFINE_PRECOND_JACOBI);
	}
	if(CheckBox_Precond_TwoLevels->GetValue() == true) {
		m_clsGlobalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_SUPER;
	}
	else {
		m_clsGlobalVars.m_ucPrecondType &= (~AUTOREFINE_PRECOND_SUPER);
	}
	TextCtrl_Dim_Precond_Two_Levels->GetValue().ToULong(&tmpULong);
	m_clsGlobalVars.m_uiSuperPreDim = (unsigned int) tmpULong;
	m_clsGlobalVars.m_bOutputCharge = CheckBox_Output_Charge->GetValue();
	m_clsGlobalVars.m_bOutputCapMtx = CheckBox_Output_Capacitance->GetValue();
	m_clsGlobalVars.m_bDumpResidual = CheckBox_Dump_Residual->GetValue();
	m_clsGlobalVars.m_bDumpTimeMem = CheckBox_Dump_Time_Mem->GetValue();
	m_clsGlobalVars.m_bVerboseOutput = CheckBox_Verbose_Output->GetValue();
	if(CheckBox_GalerkinScheme->GetValue() == false) {
		m_clsGlobalVars.m_cScheme = AUTOREFINE_COLLOCATION;
	}
	else {
		m_clsGlobalVars.m_cScheme = AUTOREFINE_GALERKIN;
	}
	m_clsGlobalVars.m_bOutputGeo = CheckBox_Output_Geometry->GetValue();
	m_clsGlobalVars.m_bDumpInputGeo = CheckBox_Dump_Geometry->GetValue();
}

CAutoRefGlobalVars RunDialog::GetGlobalVars()
{
	CopyDialogToVars();

	return m_clsGlobalVars;
}

void RunDialog::SetGlobalVarsFileIn(string fileIn)
{
    m_clsGlobalVars.m_sFileIn = fileIn;

	CopyVarsToDialog();
}
