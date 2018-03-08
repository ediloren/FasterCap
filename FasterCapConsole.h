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


#ifndef FASTERCAPCONSOLE_H
#define FASTERCAPCONSOLE_H

#include "FasterCapGlobal.h"
#include "Solver/SolverGlobal.h"

class FasterCapConsole
{
    public:
#ifdef FCG_HEADLESS
        int main(int& argc, char **argv, char bOption='\0');
#else
        int main(int& argc, wxChar **argv, char bOption='\0');
#endif // FCG_HEADLESS        int main(int& argc, wxChar **argv, char bOption = '\0');
        bool ParseCmdLine(const char *commandStr, CAutoRefGlobalVars &globalVars, wxString &errMsg);

    protected:
        int getSubstring(const char *buffer, char *substr, int *skip);

    private:
};

#endif // FASTERCAPCONSOLE_H
