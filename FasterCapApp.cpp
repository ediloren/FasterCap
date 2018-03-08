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
#include "FasterCapApp.h"

//(*AppHeaders
#include "FasterCapMain.h"
#include <wx/image.h>
//*)
#include <wx/fileconf.h>
#include <wx/wfstream.h>
#include <wx/stdpaths.h>

#include "FasterCapConsole.h"
#include "FasterCapGlobal.h"
#include "Solver/SolveCapacitance.h"

#include <iostream>
#include <stdio.h>
// wxArchiveFSHandler
#include <wx/fs_arc.h>

#ifdef __WXMSW__
#  include <wx/msw/registry.h>
#endif //__WXMSW__

// samples directory
#define FCA_SAMPLES_DIR     "Samples"


// globals
bool g_bIsConsole;


IMPLEMENT_APP(FasterCapApp)
/*
BEGIN_EVENT_TABLE(FasterCapApp, wxApp)
    EVT_IDLE(FasterCapApp::OnIdle)
END_EVENT_TABLE()
*/

bool FasterCapApp::Initialize(int& argc, wxChar **argv)
{
    int i;

    g_bIsConsole = false;

	for(i=0; i<argc; i++) {
		// if started as console application
		if(argv[i][0] == '-') {
			if(argv[i][1] == 'b') {
				g_bIsConsole = true;
			}
		}
	}

    if(g_bIsConsole == false) {
        return wxApp::Initialize(argc, argv);
    }
    else {
        return wxAppConsole::Initialize(argc, argv);
    }
}

void FasterCapApp::CleanUp()
{
    if(g_bIsConsole == false) {
        wxApp::CleanUp();
    }
    else {
        wxAppConsole::CleanUp();
    }
}

MixAppTraits *FasterCapApp::CreateTraits() {
    return new MixAppTraits(g_bIsConsole);
}

bool FasterCapApp::OnInit()
{
	int i;
	bool ret, isAutomated, helpRet;
	wxString inputArgs, errMsg, dirPath;
	CSolveCap solveMain;
	CAutoRefGlobalVars globalVars;
	FasterCapConsole console;
	char bOption;
	wxFileName fileName;

	// init global reference
	Globals::SetApp(this);

	// init vars
	//

	// start as not first use
	m_lIsFirstUse = 0;
	// no help
	m_pHelp = NULL;

	// init values for the registry
	//
	// Sets the name of application's vendor. The name will be used in registry access.
	// (for example, this is the base key for wxConfig; remark, in this case is under HKCU
	// and not HKLM)
	SetVendorName(wxT("FastFieldSolvers"));
	// Sets the name of the application. The name may be used in dialogs
	SetAppName(wxT("FasterCap"));


	// init paths (if Linux, will be null, to possibly default in launch directory)
	m_strBasePath = wxT("");
	m_strFasterCapPath = wxT("");
	m_strSamplePath = wxT("");
	m_strLicenseTextPath = wxT("");

#ifdef __WXMSW__

	//
	// retrieve values from the registry in MS Win version
	//

	wxRegKey keyBase("HKEY_LOCAL_MACHINE\\Software\\FastFieldSolvers");
	wxRegKey keyCommonSet(keyBase, "Settings");
	wxRegKey keyFasterCapSet(keyBase, "FasterCap\\Settings");

	if( keyBase.Exists() == false || keyCommonSet.Exists() == false || keyFasterCapSet.Exists() == false) {
		SysMsg("Cannot find FastFieldSolvers settings in the Registry\nPlease try installing the software again", "Error", wxICON_ERROR);
	}
	else {
		if( keyCommonSet.QueryValue(wxT("Path"), m_strBasePath) == false) {
			SysMsg("Cannot access FastFieldSolvers settings in the Registry\nPlease try installing the software again", "Error", wxICON_ERROR);
		}
		if( keyFasterCapSet.QueryValue(wxT("Path"), m_strFasterCapPath) == false) {
			SysMsg("Cannot access FastFieldSolvers settings in the Registry\nPlease try installing the software again", "Error", wxICON_ERROR);
		}
		if( keyFasterCapSet.QueryValue(wxT("SamplePath"), m_strSamplePath) == false) {
			SysMsg("Cannot access FastFieldSolvers settings in the Registry\nPlease try installing the software again", "Error");
		}
		else {
			// force a file name, otherwise the last directory name is interpreted as a file name
			m_strSamplePath += wxT("\\*.*");
		}
	}


#else // Linux

	// in Linux, the base path is the same dir of the executable
	m_strBasePath = wxPathOnly(wxStandardPaths::Get().GetExecutablePath());
	m_strFasterCapPath = wxPathOnly(wxStandardPaths::Get().GetExecutablePath());

	// must use this constructor, otherwise the last dir name of the path is taken as a file name
	fileName = wxFileName(m_strFasterCapPath, wxT(""));
	fileName.AppendDir(FCA_SAMPLES_DIR);
	// no need of '*'
	//fileName.SetName(wxT("*"));
	m_strSamplePath = fileName.GetFullPath();

#endif //__WXMSW__

	// assign license text path + file name
	fileName = wxFileName(m_strBasePath, wxT(FCG_LICENSE_TEXT_FILE_NAME));
	m_strLicenseTextPath = fileName.GetFullPath();

	//
	// parse command line options to understand if FasterCap must be run:
	// - as Automated object (MS Windows only), or
	// - in GUI-less mode from the Shell, or
	// - in GUI mode, which is the default (no specific command parameter)
	//

	//g_bIsConsole = false;
	bOption = '\0';
	isAutomated = false;
	m_bIsHelp = false;
	for(i=0; i<argc; i++) {

		// if started as console application
		if(argv[i][0] == '-') {
			if(argv[i][1] == 'b') {
                // this is now already detected in Initialize(), but further options
                // are not processed there, so must be done here
				//g_bIsConsole = true;
				// store further option character
				bOption = argv[i][2];
			}
			else if(argv[i][1] == 'h') {
				m_bIsHelp = true;
			}
		}
		if(argv[i][0] == '-' || argv[i][0] == '/') {
			// if launched with automated switch, i.e. -Embedding or -Automated or /Embedding or /Automated
			if( strcmp(argv[i], "-Embedding") == 0 || strcmp(argv[i], "/Embedding") == 0 ||
                strcmp(argv[i], "-Automated") == 0 || strcmp(argv[i], "/Automated") == 0) {

				isAutomated = true;
			}
		}
	}

	// added manually outside wxSmith controlled section (otherwise wxSmith will overwrite)

	if(g_bIsConsole == false || m_bIsHelp == true) {

#ifdef __WXMSW__
		// MS Win specific: detach from the console in case the application has to be
		// launched as GUI, so the DOS shell can terminate.
		// Drawback: you see the DOS console flashing at startup
		// See below the 'else' branch for further explanations
		FreeConsole();
#endif

		// see below wxSmith generated code
		wxInitAllImageHandlers();

		//
		// init help system
		//

		// MS help system
#if wxUSE_MS_HTML_HELP
		m_pHelp = new wxCHMHelpController();

		helpRet = m_pHelp->Initialize(m_strFasterCapPath + wxT("\\FasterCapHelp"));

		if( !helpRet ) {
			wxMessageBox(wxT("Failed initializing help system (MS HTML Help system Initialize() failed)"));
		}

		// if using CHM help, we cannot use the '-h' option, since
		// m_pHelp->GetFrame() is not supportedm_strBasePath
		if(m_bIsHelp == true) {
			m_bIsHelp = false;
		}

#else
		m_pHelp = new wxHtmlHelpController(wxHF_DEFAULT_STYLE|wxHF_OPEN_FILES);

		// Start up the handler for Zip files.
		// since (a) AddHandler is a static method, and (b) the handlers are deleted in wxFileSystem's destructor
		// you can pass the parameter (wxFileSystemHandler *) in this way.
		// without this initialization (that works globally for the class since the method is static)
		// you cannot open zip files, for instance for the help system (where all the help is in one
		// compressed file)
		wxFileSystem::AddHandler(new wxArchiveFSHandler);

		// use the registry or files (Unix) to store configuration variables for the help (? no clear idea what)
		//m_clsHelp.UseConfig(wxConfig::Get());
		// path for storing temporary files - cached binary versions of index and contents files
		// not automatically deleted on exit. Used to improve help speed
		//m_clsHelp.SetTempDir(wxT("."));
		// add zip collection of help files
		// remark: following lines not needed, seems that without path
		// the function AddBook() defaults to executable path anyway
		//wxString exePath = wxFileName(wxStandardPaths::Get().GetExecutablePath()).GetPath();
		//helpRet = m_clsHelp.AddBook(wxFileName(exePath, wxT("FasterCapHelp.zip")));
		helpRet = m_pHelp->AddBook(wxFileName("FasterCapHelp.zip"));
		// other possibilty, if the help files are not zipped
		//helpRet = m_clsHelp.AddBook(wxFileName("FasterCapHelp/FasterCap.hhp"));
		if (! helpRet) {
			wxMessageBox(wxT("Failed initializing help system (AddBook() failed)"));
		}

		// if using CHM help, we cannot use the '-h' option, since
		// m_pHelp->GetFrame() is not supported; therefore the following section
		// is covered by #if wxUSE_MS_HTML_HELP as well

		if(m_bIsHelp == true) {
			// Don't exit on frame deletion, since the help window is programmed
			// to cause the app to exit even if it is still open. We need to have the app
			// close by other means.
			// Note by Enrico: seems not true - if set to false, FasterCap never exits,
			// even calling ExitMainLoop()
			//SetExitOnFrameDelete(false);
			DisplayHelp();
			SetTopWindow(m_pHelp->GetFrame());
		}
		else {
#endif

// must comment 'outside' because this section is managed bu wxSmith, and we need to customize it
// (if done directly in the below section, it would be overwritten by wxSmith)
// Unfortunately, this confuses the Source Code Formatter (AStyle) of Code::Blocks,
// but it is a necessary evil
		/*
		        //(*AppInitialize
		        bool wxsOK = true;
		        wxInitAllImageHandlers();
		        if ( wxsOK )
		        {
		        	FasterCapFrame* Frame = new FasterCapFrame(0);
		        	Frame->Show();
		        	SetTopWindow(Frame);
		        }
		        //*)
		*/
		m_pFasterCapFrame = new FasterCapFrame(0, isAutomated);
		// if automated, create the frame but do not show it
		if(isAutomated == false) {
			m_pFasterCapFrame->Show();
		}
		SetTopWindow(m_pFasterCapFrame);
#if !wxUSE_MS_HTML_HELP
	}
#endif
}
else {
	// remark: to enable a console application, i.e. that does not detach from
	// the DOS shell once launched from a DOS shell, you need to change the Code::Blocks option
	// in Project->Properties->Bulid Targets->Type to 'Console Application'.
	// In this case however it would never detach from the shell,
	// so you would have an open DOS shell also for GUI applications
	//
	// You may also try the viceversa, i.e. start as GUI and then
	// attach again to the console, but in this case the
	// console file handlers for stdin / out / err must be reopened, and
	// in any case the prompt has already returned, writing a new prompt line
	//if(AttachConsole(ATTACH_PARENT_PROCESS) == 0 ) {
	//    wxMessageBox("The process was started with the Console option set, but it cannot detect any launching console\n(maybe you started the application from an icon shortcut passing the wrong command line parameters?)");
	//}

    console.main(argc, argv, bOption);
}

// return true to contine execution; console vs. GUI is resolved inside OnRun()
return true;
}

int FasterCapApp::OnRun()
{
	if(g_bIsConsole == true) {
	    // must explicitly call OnExit(), since the main event loop is never entered. Is this true? Seems not true.
	    // OnExit() does some clean-up.
        //OnExit();
		// if this was called as GUI-less shell application, just return the exit code
		return m_iRetStatus;
	}
	else {
		// else call the base member
		return wxApp::OnRun();
	}
}

/*
void FasterCapApp::OnIdle(wxIdleEvent& event)
{
    // only if open as help only
    if(m_bIsHelp == true) {
        // in case of help only, the top window has been assigned
        // as the help window; so if it does not exist any more,
        // the user has closed the help window

        if (!GetTopWindow()) {
            ExitMainLoop();
        }
        // don't care about the specific event
        event.Skip();
        // wait (to lower the load on the system, otherwise busy throwing idle events)
        ::wxMilliSleep(10);
        // and request to post another idle event, so we monitor
        // the status contiunously
        event.RequestMore();
    }
}
*/

int FasterCapApp::OnExit()
{
	if(m_pHelp != NULL) {
		delete m_pHelp;
	}

	return 0;
}

void FasterCapApp::DisplayHelp()
{
// MS help system
#if wxUSE_MS_HTML_HELP
	m_pHelp->DisplaySection(wxT("WelcometoFasterCap.htm"));
#else
	// open help to main page
	m_pHelp->Display(wxT("WelcometoFasterCap.htm"));
#endif
}

bool FasterCapApp::IsFirstUse()
{
	if(m_lIsFirstUse == 1) {
		return true;
	}
	else {
		return false;
	}
}

void FasterCapApp::SetNotFirstUse()
{
	m_lIsFirstUse = 0;
}

wxString FasterCapApp::GetBasePath()
{
	return m_strBasePath;
}

wxString FasterCapApp::GetSamplePath()
{
	return m_strSamplePath;
}

wxString FasterCapApp::GetLicenseTextPath()
{
	return m_strLicenseTextPath;
}

///////////////////////////////////////////////////
// MixAppTraits class
///////////////////////////////////////////////////

MixAppTraits::MixAppTraits(bool is_console) : wxGUIAppTraits(), m_bIsConsole(is_console)
{
}

wxLog *MixAppTraits::CreateLogTarget() {
    if (m_bIsConsole == false) {
        return wxGUIAppTraits::CreateLogTarget();
    }
    else {
        return new wxLogStderr;
    }
}
