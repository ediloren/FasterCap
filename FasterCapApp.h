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


#ifndef FASTERCAPAPP_H
#define FASTERCAPAPP_H

#include <wx/app.h>
#include "FasterCapMain.h"

// generic help system
#include <wx/html/helpctrl.h>
// MS help system
#if wxUSE_MS_HTML_HELP
#   include <wx/msw/helpchm.h>
#endif
#include <wx/apptrait.h>

// This class is used in conjunction with FasterCapApp to enable
// a pure console mode.
// The new AppTraits is needed to be sure that the log target is the console
// and not the wxLogGui
class MixAppTraits : public wxGUIAppTraits
{
    public:
        MixAppTraits(bool is_console);
        wxLog *CreateLogTarget();

    protected:
        bool m_bIsConsole;
};

// FasterCapApp enables also a pure console mode that can be used
// GUI-less without even
// initializing any GUI related stuff, so it can be used in
// a remote console without exported DISPLAY in Linux
// This is not an issue in Win as in this case the console
// is always opened at start-up and the windows created only
// on demand.
// Similar to what shown at
// http://compgroups.net/comp.soft-sys.wxwindows/hybrid-gui-console-app/2686190
class FasterCapApp : public wxApp
{
    public:
        inline FasterCapFrame *GetFasterCapFrame() {
            return m_pFasterCapFrame;
        }
        void DisplayHelp();
        bool IsFirstUse();
        void SetNotFirstUse();
        wxString GetBasePath();
        wxString GetSamplePath();
        wxString GetLicenseTextPath();

        int m_iRetStatus;

protected:
        // overrides
        bool Initialize(int& argc, wxChar **argv);
        void CleanUp();
        MixAppTraits *CreateTraits();
        // virtuals
        virtual bool OnInit();
        virtual int OnRun();
//        virtual void OnIdle(wxIdleEvent& event);
        virtual int OnExit();

        FasterCapFrame *m_pFasterCapFrame;
        long m_lIsFirstUse;
        wxString m_strBasePath;
        wxString m_strFasterCapPath;
        wxString m_strSamplePath;
        wxString m_strLicenseTextPath;
#if wxUSE_MS_HTML_HELP
        // MS help support
        wxCHMHelpController *m_pHelp;
#else
        // generic help support
        wxHtmlHelpController *m_pHelp;
#endif
        bool m_bIsHelp;

//        DECLARE_EVENT_TABLE()
};
#endif // FASTERCAPAPP_H
