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

#include <wx/msgdlg.h>

//(*InternalHeaders(FasterCapFrame)
#include <wx/artprov.h>
#include <wx/bitmap.h>
#include <wx/font.h>
#include <wx/intl.h>
#include <wx/image.h>
#include <wx/string.h>
//*)

// to use standard paths in a portable way through wxStandardPaths
#include <wx/stdpaths.h>
// virtual file system handlers definitions for wxFileSystemHandler
#include <wx/fs_arc.h>
// include the application icon
#include "res/FasterCap_16x16.xpm"

// for COM Automation
#ifdef __WXMSW__
#   include <windows.h>
#   include <objbase.h>
#   include <initguid.h>
#   include <wx/msw/ole/oleutils.h>
#endif //__WXMSW__

#include "FasterCapMain.h"
#include "FasterCapGlobal.h"
#include "Solver/SolveCapacitance.h"
#include "FasterCapConsole.h"

// about dialog class
#include "AboutBox.h"
// license dialog class
#include "LicenseBox.h"

#ifdef DEBUG_TEST_POT
#  include "test.h"
#endif //DEBUG_TEST_POT

// wx.pdf version 2.8.12, page 2100
// for custom WINDOW ID generation
#define FCM_OUTPUT_TEXT_ID      wxID_HIGHEST + 1
#define FCM_THREAD_END_ID       wxID_HIGHEST + 2
#define FCM_THREAD_TERM_ID       wxID_HIGHEST + 3

//
// global variables
//
// variables to control the working thread
//
// controls the exit from the infinite while() loop of the worker thread
volatile bool g_bExitFCThread = false;
// global pointer to main frame (for log functions)
FasterCapFrame *g_pFrame = NULL;


#ifdef __WXMSW__

#define FASTERCAPMAIM_AUTOWRAP_BUFFER   256

/////////////////////////////////////////////////////////////////////////////
// FasterCap COM Automation
/////////////////////////////////////////////////////////////////////////////

// GUID defines are here because if used in a header file, cause a 'multiple inclusion error'
// of the linker, since DEFINE_GUID is a variable definition, and the variable is not
// declared extern, so it is defined multiple times

// FasterCap type library's GUID
// {F99B4BC5-434E-4d11-8656-AE0AC816E586}
DEFINE_GUID(CLSID_TypeLib, 0xf99b4bc5, 0x434e, 0x4d11, 0x86, 0x56, 0xae, 0xa, 0xc8, 0x16, 0xe5, 0x86);

// FasterCap object's GUID
// {4ACA654C-C56C-47d5-86F4-7B721F06E056}
DEFINE_GUID(CLSID_IFasterCap, 0x4aca654c, 0xc56c, 0x47d5, 0x86, 0xf4, 0x7b, 0x72, 0x1f, 0x6, 0xe0, 0x56);

// FasterCap VTable's GUID
// {80FDCF45-7F79-4409-8195-4FF9BB4BD5A3}
DEFINE_GUID(IID_IFasterCap, 0x80fdcf45, 0x7f79, 0x4409, 0x81, 0x95, 0x4f, 0xf9, 0xbb, 0x4b, 0xd5, 0xa3);


// A count of how many objects we created (when clients call the IClassFactory object's CreateInstance()),
// which have not yet been Release()'d by the client
static LONG g_LOutstandingObjects;

// pointer to the type library's TYPEINFO
static ITypeInfo *IFasterCapTypeInfo;

// Counter of the number of apps that have locked our server
// via our IClassFactory object's LockServer()
static LONG g_LLockCount;

// our class factory
IFasterCapClassFactory IFCClassFactoryObj;

// register number of the COM interface
DWORD g_DWRegisterNum;

#endif //__WXMSW__

//helper functions
enum wxbuildinfoformat {
	short_f, long_f
};

wxString wxbuildinfo(wxbuildinfoformat format)
{
	wxString wxbuild(wxVERSION_STRING);

	if (format == long_f ) {
#if defined(__WXMSW__)
		wxbuild << _T("-Windows");
#elif defined(__UNIX__)
		wxbuild << _T("-Linux");
#endif

#if wxUSE_UNICODE
		wxbuild << _T("-Unicode build");
#else
		wxbuild << _T("-ANSI build");
#endif // wxUSE_UNICODE
	}

	return wxbuild;
}

//(*IdInit(FasterCapFrame)
const long FasterCapFrame::ID_TEXTCTRL1 = wxNewId();
const long FasterCapFrame::ID_PANEL1 = wxNewId();
const long FasterCapFrame::idMenuQuit = wxNewId();
const long FasterCapFrame::ID_MENUITEM1 = wxNewId();
const long FasterCapFrame::ID_MENU_RUN = wxNewId();
const long FasterCapFrame::ID_MENU_STOP = wxNewId();
const long FasterCapFrame::ID_MENU_HELP_TOPICS = wxNewId();
const long FasterCapFrame::idMenuAbout = wxNewId();
const long FasterCapFrame::ID_STATUSBAR1 = wxNewId();
const long FasterCapFrame::ID_TOOLBARITEM_RUN = wxNewId();
const long FasterCapFrame::ID_TOOLBAR1 = wxNewId();
//*)

// wx.pdf version 2.8.12, page 2103
// In wxWidgets 2.9.x seems to be replaced by wxDEFINE_EVENT macro
DEFINE_EVENT_TYPE(UWM_OUTPUT_TEXT)
DEFINE_EVENT_TYPE(UWM_THREAD_END)
DEFINE_EVENT_TYPE(UWM_THREAD_TERM)

BEGIN_EVENT_TABLE(FasterCapFrame,wxFrame)
	//(*EventTable(FasterCapFrame)
	//*)

	// Process a command, supplying the window
	// identifier, command event identifier, and
	// member function. Expects a member function with a
	// wxCommandEvent argument.
	// UWM_OUTPUT_TEXT is the custom event type passed with wxCommandEvent,
	// FCM_OUTPUT_TEXT_ID is the custom command ID (i.e. if UWM_OUTPUT_TEXT were identifying
	// the wxFrame menu, FCM_OUTPUT_TEXT_ID should be one of the items in the menu)
	EVT_COMMAND(FCM_OUTPUT_TEXT_ID, UWM_OUTPUT_TEXT, FasterCapFrame::OnOutputText)
	EVT_COMMAND(FCM_THREAD_END_ID, UWM_THREAD_END, FasterCapFrame::OnThreadEnd)
	EVT_COMMAND(FCM_THREAD_TERM_ID, UWM_THREAD_TERM, FasterCapFrame::OnThreadTerminated)

END_EVENT_TABLE()

FasterCapFrame::FasterCapFrame(wxWindow* parent, bool isAutomated, wxWindowID id)
{
	wxString exePath;

	m_bIsAutomated = isAutomated;

	//(*Initialize(FasterCapFrame)
	wxMenuItem* MenuItem5;
	wxMenuItem* MenuItem2;
	wxMenu* Menu3;
	wxMenuItem* MenuItem1;
	wxMenuItem* MenuItem4;
	wxMenu* Menu1;
	wxMenuItem* MenuItem3;
	wxBoxSizer* BoxSizer1;
	wxMenuBar* MenuBar1;
	wxFlexGridSizer* FlexGridSizer1;
	wxMenu* Menu2;
	wxMenu* Menu5;
	wxMenu* Menu4;

	Create(parent, id, _("FasterCap"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE, _T("id"));
	BoxSizer1 = new wxBoxSizer(wxHORIZONTAL);
	Panel1 = new wxPanel(this, ID_PANEL1, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL, _T("ID_PANEL1"));
	FlexGridSizer1 = new wxFlexGridSizer(0, 1, 0, 0);
	FlexGridSizer1->AddGrowableCol(0);
	FlexGridSizer1->AddGrowableRow(0);
	TextOutputWindow = new wxTextCtrl(Panel1, ID_TEXTCTRL1, wxEmptyString, wxDefaultPosition, wxSize(1000,400), wxTE_AUTO_SCROLL|wxTE_PROCESS_ENTER|wxTE_PROCESS_TAB|wxTE_MULTILINE|wxTE_RICH|wxTE_RICH2, wxDefaultValidator, _T("ID_TEXTCTRL1"));
	wxFont TextOutputWindowFont(10,wxTELETYPE,wxFONTSTYLE_NORMAL,wxNORMAL,false,wxEmptyString,wxFONTENCODING_DEFAULT);
	TextOutputWindow->SetFont(TextOutputWindowFont);
	FlexGridSizer1->Add(TextOutputWindow, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	Panel1->SetSizer(FlexGridSizer1);
	FlexGridSizer1->Fit(Panel1);
	FlexGridSizer1->SetSizeHints(Panel1);
	BoxSizer1->Add(Panel1, 1, wxALL|wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 0);
	SetSizer(BoxSizer1);
	MenuBar1 = new wxMenuBar();
	Menu1 = new wxMenu();
	MenuItem1 = new wxMenuItem(Menu1, idMenuQuit, _("Quit\tAlt-F4"), _("Quit the application"), wxITEM_NORMAL);
	Menu1->Append(MenuItem1);
	MenuBar1->Append(Menu1, _("&File"));
	Menu3 = new wxMenu();
	MenuItem6 = new wxMenuItem(Menu3, ID_MENUITEM1, _("Clear &All"), wxEmptyString, wxITEM_NORMAL);
	Menu3->Append(MenuItem6);
	MenuBar1->Append(Menu3, _("&Edit"));
	Menu4 = new wxMenu();
	MenuBar1->Append(Menu4, _("&View"));
	Menu5 = new wxMenu();
	MenuItem3 = new wxMenuItem(Menu5, ID_MENU_RUN, _("Run"), _("Open and run a FasterCap file"), wxITEM_NORMAL);
	Menu5->Append(MenuItem3);
	MenuItem4 = new wxMenuItem(Menu5, ID_MENU_STOP, _("Stop"), _("Stop FasterCap execution"), wxITEM_NORMAL);
	Menu5->Append(MenuItem4);
	MenuBar1->Append(Menu5, _("FasterCap"));
	Menu2 = new wxMenu();
	MenuItem5 = new wxMenuItem(Menu2, ID_MENU_HELP_TOPICS, _("Help Topics\tF1"), _("Display help"), wxITEM_NORMAL);
	Menu2->Append(MenuItem5);
	Menu2->AppendSeparator();
	MenuItem2 = new wxMenuItem(Menu2, idMenuAbout, _("About FasterCap...\tF2"), _("Show program information, version and copyright"), wxITEM_NORMAL);
	Menu2->Append(MenuItem2);
	MenuBar1->Append(Menu2, _("Help"));
	SetMenuBar(MenuBar1);
	StatusBar1 = new wxStatusBar(this, ID_STATUSBAR1, 0, _T("ID_STATUSBAR1"));
	int __wxStatusBarWidths_1[1] = { -1 };
	int __wxStatusBarStyles_1[1] = { wxSB_NORMAL };
	StatusBar1->SetFieldsCount(1,__wxStatusBarWidths_1);
	StatusBar1->SetStatusStyles(1,__wxStatusBarStyles_1);
	SetStatusBar(StatusBar1);
	ToolBar1 = new wxToolBar(this, ID_TOOLBAR1, wxDefaultPosition, wxDefaultSize, wxTB_HORIZONTAL|wxNO_BORDER, _T("ID_TOOLBAR1"));
	ToolBarItem1 = ToolBar1->AddTool(ID_TOOLBARITEM_RUN, _("Open"), wxArtProvider::GetBitmap(wxART_MAKE_ART_ID_FROM_STR(_T("wxART_FILE_OPEN")),wxART_TOOLBAR), wxNullBitmap, wxITEM_NORMAL, wxEmptyString, wxEmptyString);
	ToolBar1->Realize();
	SetToolBar(ToolBar1);
	BoxSizer1->Fit(this);
	BoxSizer1->SetSizeHints(this);

	Connect(idMenuQuit,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&FasterCapFrame::OnQuit);
	Connect(ID_MENUITEM1,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&FasterCapFrame::OnMenuEditClearAll);
	Connect(ID_MENU_RUN,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&FasterCapFrame::OnMenuRun);
	Connect(ID_MENU_STOP,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&FasterCapFrame::OnMenuStop);
	Connect(ID_MENU_HELP_TOPICS,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&FasterCapFrame::OnHelp);
	Connect(idMenuAbout,wxEVT_COMMAND_MENU_SELECTED,(wxObjectEventFunction)&FasterCapFrame::OnAbout);
	Connect(ID_TOOLBARITEM_RUN,wxEVT_COMMAND_TOOL_CLICKED,(wxObjectEventFunction)&FasterCapFrame::OnMenuRun);
	Connect(wxID_ANY,wxEVT_CLOSE_WINDOW,(wxObjectEventFunction)&FasterCapFrame::OnClose);
	//*)

	// set application icon
	SetIcon(wxIcon(FasterCap_16x16_xpm));

	//
	// create the 'Run' dialog
	//

	m_clsRunDlg = new RunDialog(this);

	//
	// start worker thread
	//

	// copy the main frame pointer to the global pointer
	g_pFrame = this;

	// signal globally that the thread can continue (no user request to terminate)
	g_bFCContinue = true;
	// and that the thread loop must go on indefinitely (no application close request)
	g_bExitFCThread = false;

	// create the condition. Must do it on the heap because there
	// is no standard constructor wxCondition::wxCondition()
	// (must then remember to delete in class destructor)
	// In wxWidgets, the wxCondition has an associated mutex
	// to control who is using the condtion, to avoid signalling
	// it while nobody is waiting. However if nobody is waiting,
	// and therefore the mutex is locked, and you try to signal,
	// the process trying to signal will block until the mutex
	// is released (because it is waiting for the mutex to be free)
	m_pRunEventCondition = new wxCondition(m_clsIsRunningMutex);
	// lock the mutex
	// The mutex will then be used to know if the worker thread is running, since
	// once the wxCondition is signalled, Wait() will lock the mutex and wake up
	// the thread, until another call to Wait() will happen.
	// So we can test with TryLock();
	m_clsIsRunningMutex.Lock();

	m_clsFCThread = new RunFCThread(this, m_pRunEventCondition, &m_clsIsRunningMutex);

	if ( m_clsFCThread->Create() != wxTHREAD_NO_ERROR ) {
		wxMessageBox(wxT("Can't create worker thread!"));
	}
	else if ( m_clsFCThread->Run() != wxTHREAD_NO_ERROR ) {
		wxMessageBox(wxT("Can't start worker thread!"));
	}

	// release the mutex, so the thread can lock it and start waiting
	m_clsIsRunningMutex.Unlock();

#ifdef __WXMSW__

	HRESULT	hr;
	wxWCharBuffer wcbuf;

#  ifdef DEBUG_REGISTER_TYPELIB

	int res;
	ITypeLib *ptlib;
	//HKEY hkey;
	wchar_t *pExePathWC;

	// register the type library
	// now this can give an error for an unpriviledged user; MS provides RegisterTypeLibraryForUser()
	// that anyway is not supported here. The w/a is to do it manually, i.e. map the HKEY_CLASSES_ROOT
	// registry subtree to the HKEY_CURRENT_USER registry subtree, as also suggested in MS KB935200

	// Remark: RegOverridePredefKey() not supported in current MinGW implementation, so must move
	// this section under the installer
	// Inno Setup has a command to do that (both in Pascal scripting, function RegisterTypeLibrary,
	// as well under [files] section, flag 'regtypelib'

	//    RegOpenKeyEx(HKEY_CURRENT_USER, "Software\Classes", 0, KEY_WRITE, &hkey);
	//    RegOverridePredefKey(HKEY_CLASSES_ROOT, hkey);

	//do original work
	exePath = wxFileName(wxStandardPaths::Get().GetExecutablePath()).GetPath();
	exePath += wxT("\\IFasterCap.tlb");

	// wc_str() creates a pointer to a temporary object that is already invalid at the next line,
	// so we must use a wxWCharBuffer to store the wide characters
	wcbuf = exePath.wc_str(wxConvLocal);
	pExePathWC = wcbuf.data();
	res = LoadTypeLib(&pExePathWC[0], &ptlib);
	if (!res) {
		res = RegisterTypeLib(ptlib, &pExePathWC[0], 0);
		if(res != S_OK) {
			wxMessageBox(wxT("Internal Error: RegisterTypeLib() returns error!"));
		}
		ptlib->Release();
	}
	else {
		wxMessageBox(wxT("Internal Error: cannot register the COM Automation type library!"));
	}

	//stop hacking
	//    RegOverridePredefKey(HKEY_CLASSES_ROOT, NULL);
	//    RegCloseKey(hkey);

#  endif //DEBUG_REGISTER_TYPELIB

	if(m_bIsAutomated == true) {
		//
		// COM automation
		//

		// Init the pointer to the TypeInfo - no TypeInfo yet loaded
		IFasterCapTypeInfo = NULL;
		// Init the global lock count
		g_LLockCount = 0;

		// Initialize COM
		hr = CoInitialize(NULL);
		if(hr != S_OK && hr != S_FALSE) {
			wxMessageBox(wxT("Internal Error: CoInitialize() failed!"));
			switch(hr) {
			case E_INVALIDARG:
				wxMessageBox(wxT("Internal Error: CoInitialize() returned E_INVALIDARG"));
				break;
			case E_OUTOFMEMORY:
				wxMessageBox(wxT("Internal Error: CoInitialize() returned E_OUTOFMEMORY"));
				break;
			case E_UNEXPECTED:
				wxMessageBox(wxT("Internal Error: CoInitialize() returned E_UNEXPECTED"));
				break;
			case RPC_E_CHANGED_MODE:
				wxMessageBox(wxT("Internal Error: CoInitialize() returned RPC_E_CHANGED_MODE"));
				break;
			default:
				wxMessageBox(wxT("Internal Error: CoInitialize() returned an unknown error"));
			}
		}

		// Because CoRegisterClassObject() will IClassFactory AddRef(), while Release()
		// is called only by CoRevokeClassObject(), let's initialize 'g_LOutstandingObjects'
		// to -1 so that we don't really count this first AddRef(), and it represents only those
		// AddRef's caused by an application using this EXE
		g_LOutstandingObjects = (LONG)-1;

		// Add this EXE to COM's Running Task table. Save the token that CoRegisterClassObject() returns so we
		// can later remove this EXE from COM's task table
//		hr = CoRegisterClassObject(CLSID_IFasterCap, (IUnknown *)&IFCClassFactoryObj, CLSCTX_LOCAL_SERVER, REGCLS_SINGLEUSE, &g_DWRegisterNum);
		hr = CoRegisterClassObject(CLSID_IFasterCap, (IUnknown *)&IFCClassFactoryObj, CLSCTX_SERVER, REGCLS_MULTIPLEUSE, &g_DWRegisterNum);
		if (hr)	{
			wxMessageBox(wxT("Internal Error: cannot add this exe to COM Automation running task list!"));
		}
	}

	// init IDispatch pointers for callbacks
	m_pDispEndCallback = NULL;
	m_pDispLogCallback = NULL;

#endif //__WXMSW__

}


FasterCapFrame::~FasterCapFrame()
{
	//(*Destroy(FasterCapFrame)
	//*)

	// destroy the 'Run' dialog
	m_clsRunDlg->Destroy();

	delete m_pRunEventCondition;

#ifdef __WXMSW__

	HRESULT hr;

	if(m_bIsAutomated == true) {
		//
		// COM Automation
		//

		// release automation callbacks
		if(m_pDispEndCallback != NULL) {
			// tell Automation to remove the reference to this object,
			// so memory can be released
			m_pDispEndCallback->Release();
		}
		if(m_pDispLogCallback != NULL) {
			// tell Automation to remove the reference to this object,
			// so memory can be released
			m_pDispLogCallback->Release();
		}

		// Release any ITypeInfo we got
		if (IFasterCapTypeInfo) {
			IFasterCapTypeInfo->Release();
		}

		// Remove this EXE from COM's Running task table
		hr = CoRevokeClassObject(g_DWRegisterNum);
		if(hr != S_OK) {
			wxMessageBox(wxT("Internal Error: CoRevokeClassObject() failed"));
		}

		// Free up COM
		CoUninitialize();
	}
#endif //__WXMSW__
}


#ifdef __WXMSW__
// as per faqmsw.htm file in wxWidget distribution,
// point "How do I handle Windows messages in my wxWidgets program?",
// to process Windows messages you must override MSWWindowProc()
// We'll do it only for WXMSW
// Remark: if you need to send a message instead, you can use the following format:
// return_value = SendMessage(HWND_BROADCAST, my_custom_discover_message, (WPARAM) GetHandle(), 0);
WXLRESULT FasterCapFrame::MSWWindowProc(WXUINT message, WXWPARAM wParam, WXLPARAM lParam)
{
	WXLRESULT rc;
	COPYDATASTRUCT *pCopyDataStruct;
	const char *rxCopyData;

	if (message == WM_COPYDATA) {
		// handle the message
		pCopyDataStruct = (PCOPYDATASTRUCT)(lParam);
		rxCopyData = (LPCSTR) (pCopyDataStruct->lpData);
		// the first 10 chars are the command that has been sent; the following are the arguments
		if (strncmp( rxCopyData, "path name ", 10) == 0) {
			// the argument is only the path and file name of the
			// FasterCap input file. The options have to be chosen locally.
			LaunchRunDialog(rxCopyData + 10);
		}
		// return status
		rc = 0;
	}
	else {
		// otherwise, pass the message to a default handler for processing
		rc = wxFrame::MSWWindowProc(message, wParam, lParam);
	}

	return rc;
}
#endif

// this function does not need to be multi-thread safe, because
// it is always called by the main process, either directly
// or in case of a message event from the worker thread
// (and the event is in itself safe: queued in the main process
// until it resumes from pre-emption)
void FasterCapFrame::OutputText(wxString text, int color)
{
	// from wxTextCtrl online docs, "wxTextCtrl Styles" section

	if(color == FCW_BLACK) {
		TextOutputWindow->SetDefaultStyle(wxTextAttr(*wxBLACK));
	}
	else if (color == FCW_RED) {
		TextOutputWindow->SetDefaultStyle(wxTextAttr(*wxRED));
	}
	else {
		// default
		TextOutputWindow->SetDefaultStyle(wxTextAttr(*wxBLACK));
	}

	TextOutputWindow->AppendText(text);
}

void FasterCapFrame::OnTextOutputWindowText(wxCommandEvent& )
{
}

// this function should be called either by the menu 'Run' command
// or by FastModel (or other controlling applications)
void FasterCapFrame::OnMenuRun(wxCommandEvent& )
{
	// default constructor creates an empty string
	string emptyString;

#ifdef DEBUG_TEST_POT
	test_pot2D();
#else
	// test if first launch, and if so
	if((Globals::GetApp())->IsFirstUse() == true) {
		emptyString = (Globals::GetApp())->GetSamplePath();
		(Globals::GetApp())->SetNotFirstUse();
	}
	LaunchRunDialog(emptyString);
#endif //DEBUG_TEST_POT
}

void FasterCapFrame::LaunchRunDialog(string fileIn)
{
	// test the mutex to understand if the worker thread is already busy.
	// if we can lock, then the thread is idle.
	// This also prevents Windows OLE messages (which are pumped
	// together with other window messages while waiting
	// for user input, see DoModal() ) to start FasterCap
	// thread in the meanwhile
	if(m_clsIsRunningMutex.TryLock() == wxMUTEX_NO_ERROR) {
		// if a file path / name was specified, set it in the run dialog
		if(fileIn.empty() == false) {
			m_clsRunDlg->SetGlobalVarsFileIn(fileIn);
		}
		// set the samples installation directory
		// as the first directory to look at
//        m_clsRunDlg->m_strSamplePath = m_strSamplePath;

		if( m_clsRunDlg->ShowModal() == wxID_OK ) {
			LaunchFasterCap(m_clsRunDlg->GetGlobalVars());
		}
		else {
			// if the user did not press OK, release the lock
			m_clsIsRunningMutex.Unlock();
		}
	}
	else {
		OutputText(wxT("\nError: 'Run' request received by User, but FasterCap is already running\n"), FCW_RED);
	}
}

// FasterCap launch through Automation
void FasterCapFrame::LaunchAutomation(CAutoRefGlobalVars globalVars)
{
	// test the mutex to understand if the worker thread is already busy.
	// if we can lock, then the thread is idle.
	// This also prevents Windows OLE messages (which are pumped
	// together with other window messages while waiting
	// for user input, see DoModal() ) to start FasterCap
	// thread in the meanwhile
	if(m_clsIsRunningMutex.TryLock() == wxMUTEX_NO_ERROR) {
		LaunchFasterCap(globalVars);
	}
	else {
		OutputText(wxT("\nError: 'Run' request received via Automation, but FasterCap is already running\n"), FCW_RED);
	}
}

// Remark: no synchronization object is needed, since it's always
// the same process which launchs FasterCap, both in case of user's
// request or OLE request through Windows messages, so access is
// never concurrent
// However, the assumption is that the mutex object 'm_clsIsRunningMutex'
// has ALREADY been locked by the calling process, that uses TryLock() also
// to test if FasterCap is already running
void FasterCapFrame::LaunchFasterCap(CAutoRefGlobalVars globalVars)
{
	wxCondError retCond;
	wxMutexError retMutex;

	// update the frame title
	SetLabel(globalVars.m_sFileIn + wxT(" - FasterCap"));

	// set the global vars to be used by the thread
	m_clsGlobalVars = globalVars;

	// signal globally that the thread can continue (no user request to terminate)
	g_bFCContinue = true;

	// signal the worker thread to start execution
	retCond = m_pRunEventCondition->Signal();

	if(retCond != wxCOND_NO_ERROR) {
		OutputText(wxT("\nInternal Error: Failed to signal the condition for unlocking worker thread (wxCondition::Signal() failed)\n"), FCW_RED);
	}

	// unlock the mutex
	retMutex = m_clsIsRunningMutex.Unlock();

	if(retMutex != wxMUTEX_NO_ERROR) {
		OutputText(wxT("\nInternal Error: Failed to unlock the mutex in LaunchFasterCap() (wxMutex::Unlock() failed)\n"), FCW_RED);
	}
}

void FasterCapFrame::OnMenuStop(wxCommandEvent& )
{
	StopFasterCap();
}

void FasterCapFrame::StopFasterCap()
{
	wxMutexError retMutex;

	// test the mutex to understand if the worker thread is already busy.
	// if we can lock, then the thread is idle
	if(m_clsIsRunningMutex.TryLock() == wxMUTEX_NO_ERROR) {
		// write error message: the worker thread is idle, nothing to stop
		OutputText(wxT("\nWarning: user break request received, but FasterCap is already idle\n"), FCW_RED);

		// unlock the mutex
		retMutex = m_clsIsRunningMutex.Unlock();
		if(retMutex != wxMUTEX_NO_ERROR) {
			wxMessageBox(wxT("Internal Error: failed to unlock the mutex (wxMutex::Unlock() failed)"));
		}
	}
	else {
		// signal to the working thread to stop execution
		g_bFCContinue = false;

		// write message we are going to stop
		OutputText(wxT("\nWarning: user break request received, stopping execution, please wait..\n"), FCW_RED);
	}
}

bool FasterCapFrame::IsFasterCapRunning()
{
	wxMutexError retMutex;
	bool ret;

	// test the mutex to understand if the worker thread is already busy.
	// if we can lock, then the thread is idle
	if(m_clsIsRunningMutex.TryLock() == wxMUTEX_NO_ERROR) {
		ret = false;
		// unlock the mutex
		retMutex = m_clsIsRunningMutex.Unlock();
		if(retMutex != wxMUTEX_NO_ERROR) {
			wxMessageBox(wxT("Internal Error: failed to unlock the mutex in IsFasterCapRunning() (wxMutex::Unlock() failed)"));
		}
	}
	else {
		ret = true;
	}

	return ret;
}

void FasterCapFrame::OnMenuEditClearAll(wxCommandEvent& )
{
	TextOutputWindow->Clear();
}

void FasterCapFrame::OnOutputText(wxCommandEvent& event)
{
	wxString text;
	int color;

	text = event.GetString();
	color = event.GetInt();

	OutputText(text, color);

#ifdef __WXMSW__
	// call the document callback function for Automation, if configured

	// Remark: automation cannot make calls during input-synchronous calls,
	// like SendMessage(); so post a message to the main process to ask
	// for Automation log. The log message is stored during the processing
	// of the previous SendMessage()
	// In wxWidgets version, wxPostMessage() should achieve this goal
	// instead of the MS Win non-portable PostMessage() used in the
	// MFC / VS FasterCap non-portable old 3.x version

	VARIANT parm[2];
	AutoBSTR textAuto;

	if(m_pDispLogCallback != NULL) {

		// call callback routine, passing the oldest element of the
		// log list as argument
		textAuto = text;
		parm[0].vt = VT_BSTR;
		parm[0].bstrVal = BSTR(textAuto);
		parm[1].vt = VT_I4;
		parm[1].lVal = (long)color;

		AutoWrap(DISPATCH_METHOD, NULL, m_pDispLogCallback, BSTR(m_clsLogCallbackName), 2, parm[1], parm[0]);

		// remark: not needed for VT_I4, and AutoBTSR is able to auto-deallocate itself upon destruction
		//VariantClear(&parm[0]);
		//VariantClear(&parm[1]);
	}
#endif //__WXMSW__
}

void FasterCapFrame::OnThreadEnd(wxCommandEvent& )
{
#ifdef __WXMSW__
	// call the document callback function for Automation, if configured
	if(m_pDispEndCallback != NULL) {
		AutoWrap(DISPATCH_METHOD, NULL, m_pDispEndCallback, BSTR(m_clsEndCallbackName), 0);
	}
#endif //__WXMSW__
}

void FasterCapFrame::OnThreadTerminated(wxCommandEvent& )
{
	// and really destroy the FasterCapFrame
	Destroy();
}

void FasterCapFrame::OnHelp(wxCommandEvent& )
{
	(Globals::GetApp())->DisplayHelp();
}

void FasterCapFrame::OnAbout(wxCommandEvent& )
{
	AboutBox aboutDlg(this);
	aboutDlg.ShowModal();
}

void FasterCapFrame::OnQuit(wxCommandEvent& )
{
	Quit();
}

void FasterCapFrame::Quit()
{
	// The Close() function simply generates a wxCloseEvent whose handler tries to close the window.
	// It doesn't close the window itself, however.
	// If the passed parameter is 'false', the window's close handler is able to veto the destruction
	// of this window, if 'true' it should not.
	Close(false);
}

void FasterCapFrame::OnClose(wxCloseEvent& event)
{
	if ( event.CanVeto() ) {
		// test the mutex to understand if the worker thread is already busy.
		// if we cannot lock, then the thread is running
		if(m_clsIsRunningMutex.TryLock() != wxMUTEX_NO_ERROR) {
			event.Veto();
			OutputText(wxT("\nError on user 'Close' request: cannot close window while FasterCap is running\n"), FCW_RED);
			return;
		}
		else {
			// we locked the mutex, so terminate the thread:
			// 1) change the global var saying we can exit the thread
			// 2) signal the wxCondition to unblock the thread
			g_bExitFCThread = true;
			m_pRunEventCondition->Signal();
			// unlock the mutex
			m_clsIsRunningMutex.Unlock();
		}
	}
	else {
		// if we cannot veto, then we signal to the thread to terminate execution
		// (will probably do it when it is finished running)
		// (remark: cannot use m_clsFCThread->Delete() because it waits for the
		// thread to terminate)
		g_bExitFCThread = true;
	}

	// close the frame
	//

	// Destroy() the frame. We may also call event.Skip(),
	// since the default event handler does call Destroy(), too
	Destroy();
	//event.Skip();
}

// ----------------------------------------------------------------------------
// worker thread
// ----------------------------------------------------------------------------

RunFCThread::RunFCThread(FasterCapFrame *pframe, wxCondition *runEventCondition, wxMutex *isRunning)
	: wxThread()
{
	m_pFrame = pframe;
	m_pRunEventCondition = runEventCondition;
	m_pIsRunningMutex = isRunning;
}

void RunFCThread::OnExit()
{
}

void *RunFCThread::Entry()
{
	CSolveCap solveMain;
	StlStringList::iterator stlListIter;

	// lock the mutex (it will be automatically unlocked by wxCondition::Wait()
	// when the worker thread is ready waiting for the start trigger)
	m_pIsRunningMutex->Lock();

	// loop until closed down by the controlling application
	while( g_bExitFCThread == false ) {
		// wait for the condition event to become signalled from the controlling thread
		m_pRunEventCondition->Wait();

		// check if we were asked to exit
		if ( g_bExitFCThread == true ) {
			break;
		}

		(Globals::GetApp())->m_iRetStatus = FC_NORMAL_END;

		(Globals::GetApp())->m_iRetStatus = solveMain.Run(m_pFrame->m_clsGlobalVars);

		// if return code is an error, print it
		solveMain.PrintRetError((Globals::GetApp())->m_iRetStatus);

		// close all open files (in case)
		// Not used for Unix: as per GNU standard C library: "this function should be used only in special situations (...)
		// It is also problematic since the standard streams will also be closed."
		// fcloseall();

		//
		// copy solve time and memory
		//

		m_pFrame->m_fSolveTime = g_fSolveTime;
		m_pFrame->m_lSolveMemory = g_clsMemUsageCopy.GetTotal();

		//
		// copy number of panels and links
		//

		m_pFrame->m_lPanelsNum = g_lPanelsNum;
		m_pFrame->m_lLinksNum = g_lLinksNum;

#ifdef __WXMSW__

		AutoVariant data;
		DWORD numElements[2];
		long index[2];

		//
		// copy capacitance and conductance matrices into safe-array
		//

		// Remark: the safe array contains VARIANTs, which are of type VT_R8 (double) in their turn;
		// safe array could contain VT_R8 directly, but since VBScript is not able
		// to deal with arrays containing anything but VARIANTs, we need to wrap them into variants.
		// See MS KB174576, HOWTO: Return Arrays from Server-Side Objects in ASP

		// destroy previous matrix
		m_pFrame->m_clsCapMatrix.Clear();
		m_pFrame->m_clsCondMatrix.Clear();
		// init capacitance matrix dimensions
		numElements[0] = g_clsCapMatrixRe.num_rows();
		numElements[1] = g_clsCapMatrixRe.num_cols();
		// create the safe-array...
		m_pFrame->m_clsCapMatrix.Create(VT_VARIANT, 2, numElements);
		m_pFrame->m_clsCondMatrix.Create(VT_VARIANT, 2, numElements);
		// copy values

		for(index[0]=0; index[0]<(long)g_clsCapMatrixRe.num_rows(); index[0]++) {
			for(index[1]=0; index[1]<(long)g_clsCapMatrixRe.num_cols(); index[1]++) {
				data = g_clsCapMatrixRe[(CLin_subscript)index[0]][(CLin_subscript)index[1]];
				m_pFrame->m_clsCapMatrix.PutElement(index, &data);
			}
		}
		for(index[0]=0; index[0]<(long)g_clsCapMatrixRe.num_rows(); index[0]++) {
			for(index[1]=0; index[1]<(long)g_clsCapMatrixRe.num_cols(); index[1]++) {
				data = g_clsCapMatrixIm[(CLin_subscript)index[0]][(CLin_subscript)index[1]];
				m_pFrame->m_clsCondMatrix.PutElement(index, &data);
			}
		}

		//
		// copy conductor names into safe-array
		//

		// Remark: safe array contains VARIANTs which are of type VT_BSTR;
		// See MS KB174576, HOWTO: Return Arrays from Server-Side Objects in ASP

		// destroy previous OLE matrix
		m_pFrame->m_clsCondNames.Clear();
		// init conductor names array dimensions
		numElements[0] = g_stlCondNames.size();
		// create the safe-array...
		m_pFrame->m_clsCondNames.Create(VT_VARIANT, 1, numElements);

		// copy values
		for(index[0]=0, stlListIter = g_stlCondNames.begin(); index[0]<(long)numElements[0] && stlListIter != g_stlCondNames.end(); index[0]++, stlListIter++) {
			data = stlListIter->c_str();
			m_pFrame->m_clsCondNames.PutElement(index, &data);
		}
#endif //__WXMSW__

		// inform main process that thread has stopped (is going to stop)
		// create the output text command event. This is used for callbacks when
		// controlled from another application (e.g. OLE in MS Windows)
		wxCommandEvent event(UWM_THREAD_END, FCM_THREAD_END_ID);
		// send in a thread-safe way (when m_pFrame is a GUI object)
		wxPostEvent(g_pFrame, event);
	}

	// inform main process that thread has terminated
	wxCommandEvent event(UWM_THREAD_TERM, FCM_THREAD_TERM_ID);
	// send in a thread-safe way (when m_pFrame is a GUI object)
	wxPostEvent(g_pFrame, event);

	return NULL;
}

// ----------------------------------------------------------------------------
// global functions accessed by the routines called by the worker thread, to print
// messages to the GUI window
// ----------------------------------------------------------------------------

int LogMsg(const char *fmt,...)
{
	wxString msg;
	int ret;
	va_list arg_ptr;
	char buf[FCM_LOG_BUF_SIZE];

    // The reason for using vsnprintf() instead of wxString::Format() or
    // wxString::PrintfV() is that, since we passed variable args,
    // then the proper casting is 'gone' beecause it's too late to do anything about
    // the arguments when they're already in va_list form
    // so wxString::PrintfV() expects strings in Unicode in a Unicode build,
    // or in UTF-8 in a non-unicode build, i.e. build-depending.
    // Since our LogMsg / ErrMsg is used with UTF-8 parameters in many cases
    // (i.e. passing a buffer of char), this does not work (interprets
    // the content of the buffer as Unicode). So we use vsnprintf()
    // that is always UTF-8 (viceversa, vswprintf() is always Unicode).
    // Note: the wx vararg functions (in this case wxString::Format)
    // accept any kind of strings. But this works only because they're not vararg
    // functions any longer, in fact, but rather pseudo-variadic templates.
    // (as per wxWidgets explanation in a forum thread)

	va_start(arg_ptr, fmt);
    ret = vsnprintf(buf, FCM_LOG_BUF_SIZE, fmt, arg_ptr);
	va_end(arg_ptr);

    msg = buf;

	// if PrintfV() returns error
	if(ret < 0) {
		msg = wxT("Error: LogMsg() is not able to print the message\n");
	}

	// if the main frame exists and this is not a console
	if(g_pFrame != NULL && g_bIsConsole == false) {
		// create the output text command event
		wxCommandEvent event(UWM_OUTPUT_TEXT, FCM_OUTPUT_TEXT_ID);
		// set color
		event.SetInt(FCW_BLACK);
		// set message
		event.SetString(msg);
		// send in a thread-safe way (when m_pFrame is a GUI object)
		wxPostEvent(g_pFrame, event);
	}
	else if(g_bIsConsole == true) {
		std::cout << msg;
	}
	else {
		wxMessageBox(wxT("Error: LogMsg() has an invalid reference to the main frame"));
		ret = -1;
	}

	return ret;
}

int ErrMsg(const char *fmt,...)
{
	wxString msg;
	int ret;
	va_list arg_ptr;
	char buf[FCM_LOG_BUF_SIZE];

	va_start(arg_ptr, fmt);
    ret = vsnprintf(buf, FCM_LOG_BUF_SIZE, fmt, arg_ptr);
	va_end(arg_ptr);

    msg = buf;

	// if PrintfV() returns error
	if(ret < 0) {
		msg = wxT("Error: ErrMsg() is not able to print the message\n");
	}

	// if the main frame exists and this is not a console
	if(g_pFrame != NULL && g_bIsConsole == false) {
		// create the output text command event
		wxCommandEvent event(UWM_OUTPUT_TEXT, FCM_OUTPUT_TEXT_ID);
		// set color
		event.SetInt(FCW_RED);
		// set message
		event.SetString(msg);
		// send in a thread-safe way (when m_pFrame is a GUI object)
		wxPostEvent(g_pFrame, event);
	}
	else if(g_bIsConsole == true) {
		std::cout << msg;
	}
	else {
		wxMessageBox(wxT("Error: ErrMsg() has an invalid reference to the main frame"));
		ret = -1;
	}

	return ret;
}

#ifdef __WXMSW__
/////////////////////////////////////////////////////////////////////////////
// FasterCap COM Automation
/////////////////////////////////////////////////////////////////////////////


//
// Implementation of the IFasterCap interface, based on IDispatch
// We implement a dual interface, based on a type library, so all members
// characteristic of the IDispatch function rely upon having a type library
// registered in the Registry. For example, GetIDsOfNames() does not have the command ID information
// corresponding to the command string embedded in the body of the function, but retrieves
// it from the type library information in the Registry.
// However, registring the type library has to be done in the installer, since
// this information is in the protected section of the Registry.
// To create the type library, use 'midl IFasterCap.idl' from a DOS shell;
// this produces the file 'IFasterCap.tlb' that must be used by the routine
// that registers the type library, and left in that location thereafter,
// since the path to the file is then stored in the Registry. Recommendation
// is to store it together with the .exe (there is also the option to embed
// the .tlb in the resources, and then reference only the .exe, but in this
// way it is clearer).
// To use 64 bits version of FasterCap, you need a typelib with 64 bits stubs.
// midl.exe compiler by defaults generates both 64bits and 32bits stubs, but
// you need to use a recent version of the midl compiler (after 6.0).
// Testing now the version 7.x shipped with Win7 SDK 7.1.
// To register the .tlb output file manually, you need the command
// regtlibv12.exe which is found, for 64bit, under Windows\Microsoft.NET\Framework64\xxx
// (where xxx is the latest version of the subfolders). Note that you need to run
// the command prompt as administrator (launch from Start wih right mouse button,
// selecting 'run as administrator').
// Note .tlb files canNOT be registered with regsvr32.exe, any version.
// Once the dual interface is available, also LibreOffice can automate the exe,
// while for MFC project Automation via LibreOffice did not work (IDispatch interface appears
// not fully implemented; not clear why in Excel this could work)
// Relevant information can be found in Jeff Glatt "COM in plain C" part 1 & 2, and
// "COM in C++" source example files by valeriyabobko, both available
// on www.codeblocks.com
//


// Implementation of IUnknown virtual functions of IFasterCap
// Every COM object's interface must have the implementation of the virtual functions of IUnknown:
// QueryInterface(), AddRef(), and Release().
//

// IFasterCap QueryInterface()
HRESULT STDMETHODCALLTYPE IFasterCap::QueryInterface(REFIID vTableGuid, void **ppv)
{
	// Check if the GUID matches IFasterCap VTable's GUID, or IUnknown's, or IDispatch's
	if (!IsEqualIID(vTableGuid, IID_IUnknown) && !IsEqualIID(vTableGuid, IID_IDispatch) && !IsEqualIID(vTableGuid, IID_IFasterCap) ) {
		// GUID not recognized, return error and set the handle to NULL
		*ppv = NULL;
		return(E_NOINTERFACE);
	}

	// return this object in the caller's handle
	*ppv = this;

	// Increment the count of callers who have an outstanding pointer to this object
	AddRef();

	return(NOERROR);
}

// IFasterCap AddRef()
ULONG STDMETHODCALLTYPE IFasterCap::AddRef()
{
	// Increment IFasterCap's reference count, and return the updated value.
	count++;

	return(count);
}

// IFasterCap Release()
ULONG STDMETHODCALLTYPE IFasterCap::Release()
{
	// decrement IFasterCap reference count
	count--;
	// if 0, then we can safely free this IFasterCap
	if (count == 0) {
//		if (string) {
//			SysFreeString(string);
//		}
		delete this;

		// decrement the count of outstanding objects
		InterlockedDecrement(&g_LOutstandingObjects);

		return(0);
	}
	return(count);
}

//
// Implementation of IDispatch virtual functions of IFasterCap
// Every COM object's interface must have the implementation of the virtual functions of IDispatch:
// GetTypeInfoCount(), GetTypeInfo(), GetIDsOfNames(), and Invoke().
//

// Helper function for the IDispatch functions below
HRESULT IFasterCap::LoadMyTypeInfo()
{
	register HRESULT	hr;
	LPTYPELIB			pTypeLib;

	// Load our type library and get a pointer to its TYPELIB
	// Remark: this calls implicitly pTypeLib->AddRef(pTypeLib)
	hr = LoadRegTypeLib(CLSID_TypeLib, 1, 0, 0, &pTypeLib);
	if( !hr ) {
		// if we could load our type library, get Microsoft's generic ITypeInfo,
		// giving it our loaded type library, and telling the funciton that
		// this is for our IFasterCap VTable, through the VTable GUID
		pTypeLib->GetTypeInfoOfGuid(IID_IFasterCap, &IFasterCapTypeInfo);
		if (!hr ) {
			// we no longer need the ptr to the TYPELIB, so call Release()
			pTypeLib->Release();
			// Since we must return the ITypeInfo pointer, let's increment
			// the reference count. Up to the caller to Release() when finished
			IFasterCapTypeInfo->AddRef();
		}
	}

	return(hr);
}


// IFasterCap GetTypeInfoCount()
HRESULT STDMETHODCALLTYPE IFasterCap::GetTypeInfoCount(UINT *pCount)
{
	*pCount = 1;
	return(S_OK);
}

// IFasterCap GetTypeInfo()
HRESULT STDMETHODCALLTYPE IFasterCap::GetTypeInfo(UINT itinfo, LCID , ITypeInfo **pTypeInfo)
{
	register HRESULT  hr;

	// Assume an error
	*pTypeInfo = 0;

	if (itinfo) {
		hr = ResultFromScode(DISP_E_BADINDEX);
	}
	else if (IFasterCapTypeInfo) {
		// if IFasterCapTypeInfo is already created, just increment its reference counter
		// Remark: we ignore the LCID since we support just one (default) language.
		IFasterCapTypeInfo->AddRef();
		hr = 0;
	}
	else {
		// load our type library
		hr = LoadMyTypeInfo();
	}

	// if correct return code (hr != 0), then return in 'pTypeInfo' the pointer
	if (!hr) {
		*pTypeInfo = IFasterCapTypeInfo;
	}

	return(hr);
}

// IFasterCap GetIDsOfNames()
HRESULT STDMETHODCALLTYPE IFasterCap::GetIDsOfNames(REFIID , LPOLESTR *rgszNames, UINT cNames, LCID , DISPID *rgdispid)
{
	register HRESULT  hr;

	if (!IFasterCapTypeInfo) {
		hr = LoadMyTypeInfo();
		if (hr) {
			return(hr);
		}
	}

	// Let OLE32.DLL's DispGetIDsOfNames() do all the real work of using our type
	// library to look up the DISPID of the requested function in our object
	return(DispGetIDsOfNames(IFasterCapTypeInfo, rgszNames, cNames, rgdispid));
}

// IFasterCap Invoke()
HRESULT STDMETHODCALLTYPE IFasterCap::Invoke(DISPID dispid, REFIID riid, LCID , WORD wFlags, DISPPARAMS *params, VARIANT *result, EXCEPINFO *pexcepinfo, UINT *puArgErr)
{
	register HRESULT  hr;

	// the "default" interface is the only one implemented
	if (!IsEqualIID(riid, IID_NULL))
		return(DISP_E_UNKNOWNINTERFACE);

	// if the TYPEINFO is not available, let's load it
	if (!IFasterCapTypeInfo) {
		hr = LoadMyTypeInfo();
		if (hr) {
			return(hr);
		}
	}

	// Let OLE32.DLL's DispInvoke() do all the real work of calling the appropriate
	// function in our object, and converting the passed arguments to the correct format
	return(DispInvoke(this, IFasterCapTypeInfo, dispid, wFlags, params, result, pexcepinfo, puArgErr));
}

//
// IFasterCap methods implementation
//

HRESULT STDMETHODCALLTYPE IFasterCap::IsRunning(VARIANT_BOOL *ret)
{
	bool isrunning;

	isrunning = ((Globals::GetApp())->GetFasterCapFrame())->IsFasterCapRunning();
	if(isrunning == true) {
		*ret = VARIANT_TRUE;
	}
	else {
		*ret = VARIANT_FALSE;
	}

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::ShowWindow()
{
	((Globals::GetApp())->GetFasterCapFrame())->Show();

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::Quit(VARIANT_BOOL *ret)
{
	bool isrunning;

	isrunning = ((Globals::GetApp())->GetFasterCapFrame())->IsFasterCapRunning();
	if(isrunning == true) {
		*ret = VARIANT_TRUE;
	}
	else {
		((Globals::GetApp())->GetFasterCapFrame())->Close();
		*ret = VARIANT_FALSE;
	}

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetCapacitance(VARIANT *ret)
{
	VariantCopy(ret, (LPVARIANT) (((Globals::GetApp())->GetFasterCapFrame())->m_clsCapMatrix) );

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetConductance(VARIANT *ret)
{
	VariantCopy(ret, (LPVARIANT) (((Globals::GetApp())->GetFasterCapFrame())->m_clsCondMatrix) );

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetCondNames(VARIANT *ret)
{
	VariantCopy(ret, (LPVARIANT) (((Globals::GetApp())->GetFasterCapFrame())->m_clsCondNames) );

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::Run(BSTR commandLine, VARIANT_BOOL *ret)
{
	bool status, parserRet;
	CAutoRefGlobalVars globalVars;
	FasterCapConsole console;
	wxString cmdLine, errMsg;

	*ret = VARIANT_FALSE;

	status = ((Globals::GetApp())->GetFasterCapFrame())->IsFasterCapRunning();

	// convert from BSTR to wxString using the helper function in wx/msw/ole/oleutils.cpp & oleutils.h
	// (not documented in the official docs)
	cmdLine = wxConvertStringFromOle(commandLine);

	if( status == true) {
		ErrMsg("Error: Trying to start FasterCap via Automation when FasterCap is already running\n");
		ErrMsg("       Command line is: %s\n", (const char*)cmdLine);
	}
	else {
		// parse command line
		//

		parserRet = console.ParseCmdLine((const char*)cmdLine, globalVars, errMsg);

		// if parser error
		if(parserRet == true) {
			// debug
			ErrMsg("Error: parsing the input command line provided via automation failed with the following error:\n");
			// print parser error
			ErrMsg((const char*)errMsg);

			(Globals::GetApp())->m_iRetStatus = FC_COMMAND_LINE_ERROR;
		}
		else {
			// set the input file in the run dialog for future possible manual runs
			((Globals::GetApp())->GetFasterCapFrame())->m_clsRunDlg->SetGlobalVarsFileIn(globalVars.m_sFileIn);
			// and launch FasterCap
			((Globals::GetApp())->GetFasterCapFrame())->LaunchAutomation(globalVars);

			*ret = VARIANT_TRUE;
		}
	}

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::Stop()
{
	((Globals::GetApp())->GetFasterCapFrame())->StopFasterCap();

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetSolveTime(float *ret)
{
	*ret = ((Globals::GetApp())->GetFasterCapFrame())->m_fSolveTime;

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetSolveMemory(long *ret)
{
	*ret = ((Globals::GetApp())->GetFasterCapFrame())->m_lSolveMemory;

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetPanelsNum(long *ret)
{
	*ret = ((Globals::GetApp())->GetFasterCapFrame())->m_lPanelsNum;

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetLinksNum(long *ret)
{
	*ret = ((Globals::GetApp())->GetFasterCapFrame())->m_lLinksNum;

	return(NOERROR);
}

HRESULT STDMETHODCALLTYPE IFasterCap::GetReturnStatus(short *ret)
{
	*ret = (Globals::GetApp())->m_iRetStatus;

	return(NOERROR);
}

// set an automation callback upon FasterCap working thread end
HRESULT STDMETHODCALLTYPE IFasterCap::SetEndCallback(LPDISPATCH callback, BSTR cbName,  VARIANT_BOOL *ret)
{
	FasterCapFrame *pFrame;

	pFrame = (Globals::GetApp())->GetFasterCapFrame();

	if(callback != NULL) {
		pFrame->m_pDispEndCallback = callback;
		pFrame->m_clsEndCallbackName = cbName;
		// tell Automation that there is a new reference to this object
		// so that it does not prematurely release it;
		// otherwise, the reference to the 'callback' object would go
		// out of scope once the SetEndCallback() method is complete,
		// and Automation therefore would dereference and releases the memory.
		pFrame->m_pDispEndCallback->AddRef();
		*ret = VARIANT_TRUE;
	}
	else {
		// no callback anymore
		pFrame->m_pDispEndCallback->Release();
		pFrame->m_pDispEndCallback = NULL;
		*ret = VARIANT_FALSE;
	}

	return(NOERROR);
}

// set an automation callback upon FasterCap working thread log messages
HRESULT STDMETHODCALLTYPE IFasterCap::SetLogCallback(LPDISPATCH callback, BSTR cbName,  VARIANT_BOOL *ret)
{
	FasterCapFrame *pFrame;

	pFrame = (Globals::GetApp())->GetFasterCapFrame();

	if(callback != NULL) {
		pFrame->m_pDispLogCallback = callback;
		pFrame->m_clsLogCallbackName = cbName;
		// tell Automation that there is a new reference to this object
		// so that it does not prematurely release it;
		// otherwise, the reference to the 'callback' object would go
		// out of scope once the SetEndCallback() method is complete,
		// and Automation therefore would dereference and releases the memory.
		pFrame->m_pDispLogCallback->AddRef();
		*ret = VARIANT_TRUE;
	}
	else {
		// no callback anymore
		pFrame->m_pDispLogCallback->Release();
		pFrame->m_pDispLogCallback = NULL;
		*ret = VARIANT_FALSE;
	}

	return(NOERROR);
}


//
// Implementation of IClassFactory virtual functions.
//

// IFasterCapClassFactory AddRef()
ULONG STDMETHODCALLTYPE IFasterCapClassFactory::AddRef()
{
	// Increment the count of pointers
	InterlockedIncrement(&g_LOutstandingObjects);

	// Since we never actually allocate/free an IClassFactory,
	// there is no need to keep reference count
	return(1);
}

// IFasterCapClassFactory QueryInterface()
HRESULT STDMETHODCALLTYPE IFasterCapClassFactory::QueryInterface(REFIID factoryGuid, void **ppv)
{
	// Check that the caller wants either an IUnknown or an IClassFactory.
	if (IsEqualIID(factoryGuid, IID_IUnknown) || IsEqualIID(factoryGuid, IID_IClassFactory)) {
		AddRef();

		// Return a pointer to this IFasterCapClassFactory
		*ppv = this;

		return(NOERROR);
	}

	// GUID is unknown
	*ppv = 0;
	return(E_NOINTERFACE);
}

// IFasterCapClassFactory Release()
ULONG STDMETHODCALLTYPE IFasterCapClassFactory::Release()
{
	ULONG count;

	// One less object that an app has not yet Release()'ed
	count = InterlockedDecrement(&g_LOutstandingObjects);

	// If we can unload this EXE now, then post a WM_QUIT message so WinMain
	// will drop out of the loop and this EXE terminates
	//if (!DllCanUnloadNow()) PostMessage(0, WM_QUIT, 0, 0);
	// TBC warning: if the COM object is released, we should quit FasterCap!

	return(count);
}

// IFasterCapClassFactory CreateInstance()
// Called by someone who has a pointer to this implementation
// of the IClassFactory object and now wants to create and retrieve a pointer to IFasterCap
HRESULT STDMETHODCALLTYPE IFasterCapClassFactory::CreateInstance(IUnknown *punkOuter, REFIID vTableGuid, void **objHandle)
{
	HRESULT  hr;
	register IFasterCap *FasterCapObj;

	// init the handler to zero, in case of error
	*objHandle = 0;

	// signal we don't support aggregation
	if (punkOuter) {
		hr = CLASS_E_NOAGGREGATION;
	}
	else {
		// allocate a IFasterCap object
		FasterCapObj = new IFasterCap();
		if (!FasterCapObj) {
			hr = E_OUTOFMEMORY;
		}
		else {
			// Increment the reference count so we can call Release() below and
			// it will deallocate only if there is an error with QueryInterface()
			FasterCapObj->AddRef();

			// Return a pointer to the IFasterCap object we created
			// To do that, let's use IFasterCap QueryInterface() function:
			// it checks the GUID that was passed by the caller and increments the
			// reference count
			hr = FasterCapObj->QueryInterface(vTableGuid, objHandle);

			// Decrement reference count. NOTE: If there was an error in QueryInterface()
			// then Release() will be decrementing the count back to 0 and will free the
			// IFasterCap object for us.
			FasterCapObj->Release();

			// If success, increment the static object count, that tracks if the
			// application should stay loaded, of if it can close down
			if (!hr) {
				InterlockedIncrement(&g_LOutstandingObjects);
			}
		}
	}

	return(hr);
}

// IFasterCapClassFactory LockServer()
HRESULT STDMETHODCALLTYPE IFasterCapClassFactory::LockServer(BOOL flock)
{
	if (flock) {
		InterlockedIncrement(&g_LLockCount);
	}
	else {
		InterlockedDecrement(&g_LLockCount);
	}

	return(NOERROR);
}

// AutoWrap helper function to implement a client Automation interface
//
// It is used here to call callback functions by our Automation server
//
// For information about callbacks, see online help under:
//
// - HOWTO: Automate Excel From VC++ Without Using MFC, KB216686
// - b2c.exe VisualBasic to VC++ Automation client converter by Microsoft (web site)
// - Building OLE Automation Service Components in Visual Basic and Visual C++
//   (Hotel Manager Example) and related example project in VBasic
// - Locating Resources To Study Automation, KB152023
// - Office Automation Using VC++, KB196776
// - Platform SDK -> COM and ActiveX -> Automation -> Overview of Automation ->
//   How Do Clients and Objects Interact?
// - Platform SDK -> COM and ActiveX -> Automation -> Accessing ActiveX Objects ->
//   Creating Applications and Tools That Access Objects ->
//   Accessing Members Through IDispatch / Accessing Members Through VTBLs
// - Inside OLE 2dn Edition
HRESULT AutoWrap(int autoType, VARIANT *pvResult, IDispatch *pDisp, LPOLESTR ptName, int cArgs...)
{
	// Begin variable-argument list...
	va_list marker;
	va_start(marker, cArgs);

	if(!pDisp) {
		ErrMsg("Error: NULL IDispatch passed to AutoWrap()");
		return(E_POINTER);
	}

	// Variables used...
	DISPPARAMS dp = { NULL, NULL, 0, 0 };
	DISPID dispidNamed = DISPID_PROPERTYPUT;
	DISPID dispID;
	HRESULT hr;
	char szName[FASTERCAPMAIM_AUTOWRAP_BUFFER];

	// Convert down to ANSI
	WideCharToMultiByte(CP_ACP, 0, ptName, -1, szName, FASTERCAPMAIM_AUTOWRAP_BUFFER, NULL, NULL);

	// Get DISPID for name passed...
	hr = pDisp->GetIDsOfNames(IID_NULL, &ptName, 1, LOCALE_USER_DEFAULT, &dispID);
	if(FAILED(hr)) {
		ErrMsg("Error: IDispatch::GetIDsOfNames(\"%s\") failed with error 0x%08lx", szName, hr);
		ErrMsg("       Most probably there is no automation interface available in the client application");
		return hr;
	}

	// Allocate memory for arguments...
	VARIANT *pArgs = new VARIANT[cArgs+1];
	// Extract arguments...
	for(int i=0; i<cArgs; i++) {
		pArgs[i] = va_arg(marker, VARIANT);
	}

	// Build DISPPARAMS
	dp.cArgs = cArgs;
	dp.rgvarg = pArgs;

	// Handle special-case for property-puts!
	if(autoType & DISPATCH_PROPERTYPUT) {
		dp.cNamedArgs = 1;
		dp.rgdispidNamedArgs = &dispidNamed;
	}

	// Make the call!
	hr = pDisp->Invoke(dispID, IID_NULL, LOCALE_SYSTEM_DEFAULT, autoType, &dp, pvResult, NULL, NULL);
	if(FAILED(hr)) {
		ErrMsg("Error: IDispatch::Invoke(\"%s\"=%08lx) failed w/err 0x%08lx", szName, dispID, hr);
		ErrMsg("       Most probably there is no automation interface available in the client application");
		return hr;
	}
	// End variable-argument section...
	va_end(marker);

	delete [] pArgs;

	return hr;
}


#endif //__WXMSW__


