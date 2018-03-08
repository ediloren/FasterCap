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


#include "FasterCapGlobal.h"

#ifndef FCG_HEADLESS
    #include <wx/msgdlg.h>
#endif // !FCG_HEADLESS

// signals FasterCap to abort simulation and start waiting again for commands
volatile bool g_bFCContinue = true;

#ifdef DEBUG_LOG_ENABLE

int DebugMsg(const char *fmt,...)
{
    int ret = 0;

	wxString msg;
	va_list arg_ptr;
	wxDateTime currtime;

	va_start(arg_ptr, fmt);

	ret = msg.PrintfV(fmt, arg_ptr);

	va_end(arg_ptr);

	// if PrintfV() returns error
	if(ret < 0) {
		msg = wxT("Error: DebugMsg() is not able to print the message\n");
	}

	// add timestamp
	currtime.SetToCurrent();
	msg = currtime.FormatISODate() + wxT("@") + currtime.FormatISOTime() + wxT(" ") + msg;

	std::cout << msg;
	std::cout.flush();

	return ret;
}

#else //not DEBUG_LOG_ENABLE

int DebugMsg(const char *,...)
{

	return 0;
}

#endif //DEBUG_LOG_ENABLE


#ifndef FCG_HEADLESS

// Static data members must be initialized at file scope, even if private.
FasterCapApp *Globals::m_pApp = NULL;

// remark: cannot be used before initializing 'g_bIsConsole'
// also, in console mode 'style' is ignored
void SysMsg(wxString message, wxString caption, long style)
{
    if(g_bIsConsole == true) {
        std::cout << caption << " : " << message << std::endl;
	}
	else {
		wxMessageBox(message, caption, style);
    }
}

// in case of GUI, the LogMsg(), ErrMsg() functions etc. are implemented in FasterCapApp.cpp
 
#else // if FCG_HEADLESS

// ----------------------------------------------------------------------------
// global functions for printing messages to the console window
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

	std::cout << msg;

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

	std::cout << msg;

	return ret;
}


#endif // !FCG_HEADLESS

