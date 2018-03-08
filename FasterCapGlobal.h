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


// CPP header file containing function prototypes linking
// FasterCap with the GUI I/O window.
// Actual function bodies, when needed, are in FasterCapMain.cpp

#ifndef FASTERCAPGLOBAL_H
#define FASTERCAPGLOBAL_H

#include <stdio.h>
#include <stdarg.h>

#include <wx/string.h>

// debug defines
//#define DEBUG_REGISTER_TYPELIB
//#define DEBUG_TEST_POT
//#define DEBUG_LOG_ENABLE

//#define DEBUG_DUMP_BASIC
//#define DEBUG_DUMP_POT
//#define DEBUG_DUMP_UNCOMP_POT
//#define DEBUG_DUMP_OTHER

// define headless
//#define FCG_HEADLESS

// version and copyright header messages
#define FCG_LIC_STD_HEADER1     "FasterCap License"
#define FCG_HEADER_VERSION      "FasterCap version 6.0.7"
#define FCG_HEADER_COPYRIGHT    "Copyright 2019 FastFieldSolvers S.R.L."
#define FCG_HEADER_WEBSITE      "http://www.fastfieldsolvers.com, All Rights reserved"
// license text file name
#define FCG_LICENSE_TEXT_FILE_NAME      "LICENCE.txt"

// must include global definitions for CMemoryUsage class definition
// and FC return code defines, etc.
#include "Solver/SolverGlobal.h"

#ifndef FCG_HEADLESS
// must include for App pointer global visibility
#include "FasterCapApp.h"
#endif // !FCG_HEADLESS

// Needed for CLin_Matrix
#include "LinAlgebra/Mtx.h"

#define MAX_TITLE_LENGHT 64

// color IDs for OutputText() in FasterCapMain.cpp
#define FCW_BLACK   1
#define FCW_RED     2

// temp buffer for LogMsg() and ErrMsg()
#define FCM_LOG_BUF_SIZE    1024


#ifndef FCG_HEADLESS
// global, exported variables & functions

class Globals {

public:
    inline static FasterCapApp *GetApp() {
        return m_pApp;
    }
    inline static void SetApp(FasterCapApp *fca) {
        m_pApp = fca;
    }


protected:
    static FasterCapApp *m_pApp;
};

// this function is used for system messages (console stdout or
// GUI MessageBox according to 'FasterCapApp' being a pure console
// or a GUI)
void SysMsg(wxString message, wxString caption="Message", long style = wxOK|wxCENTRE);

#endif //FCG_HEADLESS


// interrupt FasterCap worker thread
extern volatile bool g_bFCContinue;

// later
//extern char g_sTitle[64];
extern CLin_Matrix g_clsCapMatrixRe, g_clsCapMatrixIm;
extern StlStringList g_stlCondNames;
extern CMemoryUsage g_clsMemUsageCopy, g_clsMemUsage;
extern float g_fSolveTime;
extern long g_lPanelsNum;
extern long g_lLinksNum;
extern bool g_bIsConsole;

// this function collects all log messages
int LogMsg(const char *fmt,...);
// this function collects all error messages
int ErrMsg(const char *fmt,...);
// this function collects all debug messages (for console only)
int DebugMsg(const char *fmt,...);



#endif //FASTERCAPGLOBAL_H

