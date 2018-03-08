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


#ifndef FASTERCAPMAIN_H
#define FASTERCAPMAIN_H

#include "Solver/SolverGlobal.h"

//(*Headers(FasterCapFrame)
#include <wx/sizer.h>
#include <wx/menu.h>
#include <wx/textctrl.h>
#include <wx/toolbar.h>
#include <wx/panel.h>
#include <wx/frame.h>
#include <wx/statusbr.h>
//*)

// for COM Automation
#ifdef __WXMSW__
// local automation helper classes
#   include "AutomationHelper.h"
#endif //__WXMSW__

// run dialog class
#include "RunDialog.h"

// prototypes
class RunFCThread;


class FasterCapFrame: public wxFrame
{
    public:

        FasterCapFrame(wxWindow* parent, bool isAutomated, wxWindowID id = -1);
        virtual ~FasterCapFrame();
        void OutputText(wxString text, int color);
        void LaunchAutomation(CAutoRefGlobalVars globalVars);
        bool IsFasterCapRunning();
        void Quit();
        void StopFasterCap();

        CAutoRefGlobalVars m_clsGlobalVars;
		float m_fSolveTime;
		long m_lSolveMemory;
		long m_lPanelsNum;
		long m_lLinksNum;
        RunDialog *m_clsRunDlg;
#ifdef __WXMSW__
        AutoSafeArray m_clsCapMatrix;
        AutoSafeArray m_clsCondMatrix;
        AutoSafeArray m_clsCondNames;
        IDispatch *m_pDispEndCallback;
        AutoBSTR m_clsEndCallbackName;
        IDispatch *m_pDispLogCallback;
        AutoBSTR m_clsLogCallbackName;
#endif //__WXMSW__

    private:

        //(*Handlers(FasterCapFrame)
        void OnQuit(wxCommandEvent& event);
        void OnAbout(wxCommandEvent& event);
        void OnTextOutputWindowText(wxCommandEvent& event);
        void OnMenuRun(wxCommandEvent& event);
        void OnHelp(wxCommandEvent& event);
        void OnClose(wxCloseEvent& event);
        void OnMenuStop(wxCommandEvent& event);
        void OnMenuEditClearAll(wxCommandEvent& event);
        //*)
        void OnOutputText(wxCommandEvent& event);
        void OnThreadEnd(wxCommandEvent& event);
        void OnThreadTerminated(wxCommandEvent& event);

        void LaunchRunDialog(std::string fileIn);
        void LaunchFasterCap(CAutoRefGlobalVars globalVars);

        // overrides
#ifdef __WXMSW__
        WXLRESULT MSWWindowProc(WXUINT message, WXWPARAM wParam, WXLPARAM lParam);
#endif

        //(*Identifiers(FasterCapFrame)
        static const long ID_TEXTCTRL1;
        static const long ID_PANEL1;
        static const long idMenuQuit;
        static const long ID_MENUITEM1;
        static const long ID_MENU_RUN;
        static const long ID_MENU_STOP;
        static const long ID_MENU_HELP_TOPICS;
        static const long idMenuAbout;
        static const long ID_STATUSBAR1;
        static const long ID_TOOLBARITEM_RUN;
        static const long ID_TOOLBAR1;
        //*)

        //(*Declarations(FasterCapFrame)
        wxToolBar* ToolBar1;
        wxTextCtrl* TextOutputWindow;
        wxPanel* Panel1;
        wxToolBarToolBase* ToolBarItem1;
        wxStatusBar* StatusBar1;
        wxMenuItem* MenuItem6;
        //*)

        RunFCThread *m_clsFCThread;
        wxMutex m_clsIsRunningMutex;
        wxCondition *m_pRunEventCondition;
        bool m_bIsAutomated;

        DECLARE_EVENT_TABLE()
};

// ----------------------------------------------------------------------------
// worker thread
// ----------------------------------------------------------------------------

class RunFCThread : public wxThread
{
public:
    RunFCThread(FasterCapFrame *pframe, wxCondition *runEventCondition, wxMutex *isRunning);

    // thread execution starts here
    virtual void *Entry();

    // called when the thread exits - whether it terminates normally or is
    // stopped with Delete() (but not when it is Kill()ed!)
    virtual void OnExit();

public:
    FasterCapFrame *m_pFrame;

protected:
    wxCondition *m_pRunEventCondition;
    wxMutex *m_pIsRunningMutex;
};

#ifdef __WXMSW__
/////////////////////////////////////////////////////////////////////////////
// FasterCap COM Automation
/////////////////////////////////////////////////////////////////////////////

// FasterCap's VTable
#undef  INTERFACE
#define INTERFACE IFasterCap
//DECLARE_INTERFACE_ (INTERFACE, IDispatch)
struct FAR IFasterCap : public IDispatch
{
	// IUnknown functions
	STDMETHOD  (QueryInterface)		(THIS_ REFIID, void **);
	STDMETHOD_ (ULONG, AddRef)		(THIS);
	STDMETHOD_ (ULONG, Release)		(THIS);
	// IDispatch functions
	// Remark: original usage of STDMETHOD_(UINT, xxx), i.e. with a return value, gives an error, modified
	STDMETHOD (GetTypeInfoCount)(THIS_ UINT *);
	STDMETHOD (GetTypeInfo)		(THIS_ UINT, LCID, ITypeInfo **);
	STDMETHOD (GetIDsOfNames)	(THIS_ REFIID, LPOLESTR *, UINT, LCID, DISPID *);
	STDMETHOD (Invoke)			(THIS_ DISPID, REFIID, LCID, WORD, DISPPARAMS *, VARIANT *, EXCEPINFO *, UINT *);
	// Extra functions
	STDMETHOD  (IsRunning)      (THIS_ VARIANT_BOOL *);
	STDMETHOD  (ShowWindow)     (THIS);
	STDMETHOD  (Quit)           (THIS_ VARIANT_BOOL *);
	STDMETHOD  (GetCapacitance) (THIS_ VARIANT *);
	STDMETHOD  (GetConductance) (THIS_ VARIANT *);
	STDMETHOD  (GetCondNames)   (THIS_ VARIANT *);
	STDMETHOD  (Run)            (THIS_ BSTR, VARIANT_BOOL *);
    STDMETHOD  (Stop)           (THIS);
	STDMETHOD  (GetSolveTime)   (THIS_ float *);
	STDMETHOD  (GetSolveMemory) (THIS_ long *);
	STDMETHOD  (GetPanelsNum)   (THIS_ long *);
	STDMETHOD  (GetLinksNum)    (THIS_ long *);
	STDMETHOD  (GetReturnStatus)(THIS_ short *);
    STDMETHOD  (SetEndCallback) (THIS_ LPDISPATCH, BSTR, VARIANT_BOOL *);
    STDMETHOD  (SetLogCallback) (THIS_ LPDISPATCH, BSTR, VARIANT_BOOL *);

protected:
    HRESULT LoadMyTypeInfo();

    DWORD   count;
    BSTR    string;
};


// Issue with x64 crash of COM server at activation.
// By MS #defines (see WinNT.h), the DECLARE_INTERFACE_ is expanded into ' struct __declspec(novtable) '
// this is creating a crash when CoRegisterClassObject() on 'IFasterCapClassFactory' is called.
// Not expanding in this way, but as 'class' or 'struct FAR' is avoiding the crash.
// Note that __declspec(novtable) is, as per MS documentation:
// "This form of __declspec can be applied to any class declaration, but should only be applied
// to pure interface classes, that is, classes that will never be instantiated on their own."
// so it seems it makes sense it crashes. Not clear why the macro is expanded in this way,
// when MS BaseTyps.h says in the comments to the DECLARE_INTERFACE_ macro that:
// "      Example C++ expansion:
//  struct FAR IClassFactory : public IUnknown
//        NOTE: Our documentation says '#define interface class' but we use
//              'struct' instead of 'class' to keep a lot of 'public:' lines
//              out of the interfaces.  The 'FAR' forces the 'this' pointers to
//              be far, which is what we need."

// FasterCap's VTable
#undef  INTERFACE
#define INTERFACE IFasterCapClassFactory
//DECLARE_INTERFACE_(INTERFACE,IClassFactory)
struct FAR IFasterCapClassFactory : public IClassFactory
{
	// IUnknown functions
	STDMETHOD  (QueryInterface)		(THIS_ REFIID, void **);
	STDMETHOD_ (ULONG, AddRef)		(THIS);
	STDMETHOD_ (ULONG, Release)		(THIS);
	// IDispatch functions
	// Remark: original usage of STDMETHOD_(UINT, xxx), i.e. with a return value, gives an error, modified
	STDMETHOD (CreateInstance) (THIS_ LPUNKNOWN, const IID&, void**);
	STDMETHOD (LockServer) (THIS_ BOOL);
};

// AutoWrap helper function to implement a client Automation interface
HRESULT AutoWrap(int autoType, VARIANT *pvResult, IDispatch *pDisp, LPOLESTR ptName, int cArgs...);

#endif //__WXMSW__

// this is typically in a header: it just declares UWM_OUTPUT_TEXT event type
// in the header of the source file
// You can ignore the value parameter of the DECLARE_EVENT_TYPE macro since it
// used only for backwards compatibility
// wx.pdf version 2.8.12, page 2103
// In wxWidgets 2.9.x seems to be replaced by wxDECLARE_EVENT macro
DECLARE_EVENT_TYPE(UWM_OUTPUT_TEXT, -1)
DECLARE_EVENT_TYPE(UWM_THREAD_END, -1)
DECLARE_EVENT_TYPE(UWM_THREAD_TERM, -1)

#endif // FASTERCAPMAIN_H
