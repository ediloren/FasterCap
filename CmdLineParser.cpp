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


#include "CmdLineParser.h"

#define CMDLINEPARSER_MAX_INPUT_STR_LEN     4094

// returns 'true' if errors in the parsing
bool CmdLineParser::ParseCmdLine(const char *commandStr, CAutoRefGlobalVars &globalVars, wxString &errMsg)
{
	char argStr[CMDLINEPARSER_MAX_INPUT_STR_LEN], *cmdStr;
	int res, skip;
	bool cmderr;

    // clean error message
    errMsg.Clear();

	// copy pointer, not to change 'commandStr'
	// (required when the function is called by automation passing a 'const')
	cmdStr = (char *)commandStr;

	// defaults

	globalVars.Reset();

	// parse input arguments
	//

	res = getSubstring(cmdStr, argStr, &skip);
	cmdStr += skip;
	if( res == 0 || res == EOF)  {
		errMsg = wxString::Format(wxT("FasterCap launched without enough arguments! You must specify at least the input file name\n"));
		errMsg += wxString::Format(wxT("Command line is: %s\n"), cmdStr);
		return true;
	}

	cmderr = false;
	while(res != 0 && res != EOF && cmderr != true) {
		if(argStr[0] == '-') {

			// '-e' is tolerance for refinement
			if(argStr[1] == 'e') {
				if(sscanf(&(argStr[2]), "%lf", &(globalVars.m_dEps)) != 1) {
					cmderr = true;
					//errMsg = wxString::Format(wxT("%s: bad mutual potential epsilon value '%s'\n"), commandStr, &argStr[2]);
					errMsg = wxString::Format(wxT("%s: unsupported parameter '%s'\n"), commandStr, argStr);
				}
			}

			// '-s' is maximum side left after discretization
			else if(argStr[1] == 's') {
				if(sscanf(&(argStr[2]), "%lf", &(globalVars.m_dMaxDiscSide)) != 1) {
					cmderr = true;
					//errMsg = wxString::Format(wxT("%s: bad side discretization value '%s'\n"), commandStr, &argStr[2]);
					errMsg = wxString::Format(wxT("%s: unsupported parameter '%s'\n"), commandStr, argStr);
				}
			}

			// '-m' is mesh relative refinement level
			// '-mc' is the curvature coefficient for gradient meshing near sharp edges
			else if(argStr[1] == 'm') {
				if(argStr[2] == 'c') {
					if(sscanf(&(argStr[3]), "%lf", &(globalVars.m_dMeshCurvCoeff)) != 1) {
						cmderr = true;
						errMsg = wxString::Format(wxT("%s: bad mesh curvature coefficient '%s'\n"), commandStr, &argStr[2]);
					}
				}
				else if(sscanf(&(argStr[2]), "%lf", &(globalVars.m_dMeshEps)) != 1) {
					cmderr = true;
					errMsg = wxString::Format(wxT("%s: bad mesh relative refinement value '%s'\n"), commandStr, &argStr[3]);
				}
			}

			// '-d' is direct potential interaction coefficient to mesh refinement ratio
			else if(argStr[1] == 'd') {
				if(sscanf(&(argStr[2]), "%lf", &(globalVars.m_dEpsRatio)) != 1) {
					cmderr = true;
					errMsg = wxString::Format(wxT("%s: bad direct potential interaction coefficient to mesh refinement ratio value '%s'\n"), commandStr, &argStr[2]);
				}
			}

			// '-f' is Out-Of-Core free memory to link memory condition (0 = don't go OOC)
			else if(argStr[1] == 'f') {
				if(sscanf(&(argStr[2]), "%lf", &(globalVars.m_dOutOfCoreRatio)) != 1) {
					cmderr = true;
					errMsg = wxString::Format(wxT("%s: Out-Of-Core free memory to link memory condition '%s'\n"), commandStr, &argStr[2]);
				}
			}

			// '-p' is type of preconditioner
			else if(argStr[1] == 'p') {
				if(argStr[2] == 'N' || argStr[2] == 'n') {
					globalVars.m_ucPrecondType = AUTOREFINE_PRECOND_NONE;
				}
				else if (argStr[2] == 'J' || argStr[2] == 'j') {
					globalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_JACOBI;
				}
				else if (argStr[2] == 'B' || argStr[2] == 'b') {
					globalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_BLOCK;
					if(sscanf(&(argStr[3]), "%u", &(globalVars.m_uiBlockPreSize)) != 1) {
						cmderr = true;
						errMsg = wxString::Format(wxT("%s: bad preconditioner dimension '%s'\n"), commandStr, &argStr[3]);
					}
				}
				else if (argStr[2] == 'S' || argStr[2] == 's') {
					globalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_SUPER;
					if(sscanf(&(argStr[3]), "%u", &(globalVars.m_uiSuperPreDim)) != 1) {
						cmderr = true;
						errMsg = wxString::Format(wxT("%s: bad preconditioner dimension '%s'\n"), commandStr, &argStr[3]);
					}
				}
				else if (argStr[2] == 'H' || argStr[2] == 'h') {
					globalVars.m_ucPrecondType |= AUTOREFINE_PRECOND_HIER;
					if (argStr[3] == 's') {
						if(sscanf(&(argStr[4]), "%lf", &(globalVars.m_dMaxHierPreDiscSide)) != 1) {
							cmderr = true;
							errMsg = wxString::Format(wxT("%s: bad preconditioner side discretization value '%s'\n"), commandStr, &argStr[4]);
						}
					}
					else if (argStr[3] == 'e') {
						if(sscanf(&(argStr[4]), "%lf", &(globalVars.m_dHierPreEps)) != 1) {
							cmderr = true;
							errMsg = wxString::Format(wxT("%s: bad preconditioner mutual potential epsilon value  '%s'\n"), commandStr, &argStr[4]);
						}
					}
					else if(argStr[1] == 't') {
						if(sscanf(&(argStr[4]), "%lf", &(globalVars.m_dHierPreGmresTol)) != 1) {
							cmderr = true;
							errMsg = wxString::Format(wxT("%s: bad preconditioner GMRES tolerance value '%s'\n"), commandStr, &argStr[4]);
						}
					}
				}
				else {
					cmderr = true;
					errMsg = wxString::Format(wxT("%s: bad preconditioner type '%s'\n"), commandStr, &argStr[2]);
				}
			}

			// '-kc' is keep charge information after termination
			else if(argStr[1] == 'k') {
				if(argStr[2] == 'c') {
					globalVars.m_bKeepCharge = true;
				}
			}

			// '-s' is refine mesh using calculated charges
			else if(argStr[1] == 's') {
				globalVars.m_bRefineCharge = true;
			}

			// '-km' is keep mesh used in previous run
			else if(argStr[1] == 'k') {
				if(argStr[2] == 'm') {
					globalVars.m_bKeepMesh = true;
				}
			}

			// '-r' is dump residual at every GMRES step
			else if(argStr[1] == 'r') {
				globalVars.m_bDumpResidual = true;
			}

			// '-i' is dump detailed time and memory information
			else if(argStr[1] == 'i') {
				globalVars.m_bDumpTimeMem = true;
			}

			// '-v' is provide verbose output
			else if(argStr[1] == 'v') {
				globalVars.m_bVerboseOutput = true;
			}

			// '-t' is GMRES tolerance to stop iteration
			else if(argStr[1] == 't') {
				if(sscanf(&(argStr[2]), "%lf", &(globalVars.m_dGmresTol)) != 1) {
					cmderr = true;
					errMsg = wxString::Format(wxT("%s: bad GMRES iteration tolerance '%s'\n"), commandStr, &argStr[2]);
				}
			}

			// '-o' is output refined geometry to file (in FastCap2 format)
			else if(argStr[1] == 'o') {
				globalVars.m_bOutputGeo = true;
			}

			// '-g' is galerkin scheme (as opposite to default collocation)
			else if(argStr[1] == 'g') {
				globalVars.m_cScheme = AUTOREFINE_GALERKIN;
			}

			// '-a' is automatic refinement and selection of parameters
			else if(argStr[1] == 'a') {
				if(argStr[2] == 'p' || argStr[2] == 'P') {
					globalVars.m_bAutoPrecond = true;
				}
				else {
					globalVars.m_bAuto = true;
					if(sscanf(&(argStr[2]), "%lf", &(globalVars.m_dAutoMaxErr)) != 1) {
						cmderr = true;
						errMsg = wxString::Format(wxT("%s: bad max error value for auto option'%s'\n"), commandStr, &argStr[2]);
					}
				}
			}

			// '-b' is the console option; not to be passed to FasterCap engine
			else if(argStr[1] != 'b') {
				errMsg = wxString::Format(wxT("%s: illegal option -- '-%s'\n"), commandStr, &(argStr[1]));
				cmderr = true;
			}
		}
		// isn't an option, must be the input file
		else {
			globalVars.m_sFileIn = argStr;
		}

		// get next argument
		res = getSubstring(cmdStr, argStr, &skip);
		cmdStr += skip;
	}

	return cmderr;
}

// Special version of sscanf which handles also file names with spaces inside,
// provided they are surrounded by '"'.
// Note that this routine is specialized to retrieve only strings
int CmdLineParser::getSubstring(const char *buffer, char *substr, int *skip)
{
	int res, openPos, deltaClosePos, startPos;
	char *openPosPtr, *closePosPtr, tmpStr[256];

	substr[0] = '\0';

	// read a piece of string
	res = sscanf(buffer, "%s%n", tmpStr, skip);
	// if not finished
	if( res != EOF ) {
		// find if it does contain any '"'
		openPosPtr = (char *)strchr(buffer, '"');
        // if found and if within tmpStr (before any space)
		openPos = (int)(openPosPtr - buffer);
		if( openPosPtr != NULL && openPos < (int)strlen(tmpStr) ) {
			// search for the closing '"'
			closePosPtr = (char *)strchr(buffer + openPos + 1, '"');
			// if no closing '"', assume the end is the end of the
			// string in the buffer, regardless of any spaces or tabs
			if(closePosPtr == NULL) {
				// build 'substr' string by collating the different pieces,
				// skipping the '"'
				strncpy(tmpStr, buffer, openPos);
				tmpStr[openPos] = '\0';
				strcat(tmpStr, buffer + openPos + 1);
				// remove any trailing space left
				startPos = strspn(tmpStr, " \t\n");
				strcpy(substr, tmpStr + startPos);
				*skip = strlen(buffer);
			}
			else {
				// build 'substr' string by collating the different pieces,
				// skipping the '"'s
				deltaClosePos = (int)(closePosPtr - openPosPtr);
				strncpy(tmpStr, buffer, openPos);
				tmpStr[openPos] = '\0';
				strncat(tmpStr, buffer + openPos + 1, deltaClosePos - 1);
				// remove any trailing space left
				startPos = strspn(tmpStr, " \t\n");
				strcpy(substr, tmpStr + startPos);
				*skip = openPos + deltaClosePos + 1;
			}
		}
		// if no '"' at all, simply copy scanned string
		else {
			strcpy(substr, tmpStr);
		}
	}

	return res;
}
