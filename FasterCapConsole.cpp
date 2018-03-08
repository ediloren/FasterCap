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


#include "FasterCapConsole.h"

#include "Solver/SolveCapacitance.h"

#define CMDLINEPARSER_MAX_INPUT_STR_LEN     4094

#ifdef FCG_HEADLESS

int main(int argc, char **argv)
{
    int ret;

    FasterCapConsole console;

    ret = console.main(argc, argv);

    return ret;
}

int FasterCapConsole::main(int& argc, char **argv, char bOption)
#else
int FasterCapConsole::main(int& argc, wxChar **argv, char bOption)
#endif // FCG_HEADLESS
{
    bool ret;
	CSolveCap solveMain;
	int i, retStatus;
	wxString inputArgs, errMsg;
	CAutoRefGlobalVars globalVars, defGlobalVars;

	// copy input arguments in a string for proper parsing (e.g. taking care of file names with spaces)
	inputArgs.Clear();
	for(i=1; i<argc; i++) {
		inputArgs += argv[i];
		// this puts an extra trailing space in the end but that's ok anyway for the parser
		inputArgs += wxT(" ");
	}

	// dump command args for debug
	//printf("%s\n", (const char*)inputArgs);

	// parse command line
	ret = ParseCmdLine((const char*)inputArgs, globalVars, errMsg);

	// if no parser error
	if(ret == false && bOption != '?' && bOption != 'v') {
		// actually run FasterCap
		retStatus = FC_NORMAL_END;
		retStatus = solveMain.Run(globalVars);
	}
	else {
		if(bOption == '?' || bOption == 'v') {
			retStatus = FC_NORMAL_END;
			errMsg = wxT("");
		}
		else {
			retStatus = FC_COMMAND_LINE_ERROR;
		}

		LogMsg(FCG_HEADER_VERSION);
		LogMsg("\n");
		LogMsg(FCG_HEADER_COPYRIGHT);
		LogMsg(" ");
		LogMsg(FCG_HEADER_WEBSITE);
		LogMsg("\n");

		if(bOption != 'v') {
			// print error
			ErrMsg((const char*)errMsg);
			LogMsg("Usage: %s <input file> [-a<relative error>] [-ap]\n", (const char*)argv[0]);
			LogMsg("                 [-m<mesh>] [-mc<mesh curvature] [-t<tolerance>]\n");
			LogMsg("                 [-d<interaction coeff>] [-f<outofcore>] [-g]\n");
			LogMsg("                 [-pj] [-ps<dimension>] [-o] [-r] [-c] [-i] [-v]\n");
			LogMsg("                 [-b|-b?|-bv]\n");
			LogMsg("DEFAULT VALUES:\n");
			LogMsg("  -a:  Automatically calculate settings, stop when\n");
			LogMsg("       relative error is lower than <relative error>, e.g. 0.01\n");
			LogMsg("  -ap: Automatic preconditioner usage\n");
			LogMsg("  -m:  Mesh relative refinement value = %g\n", defGlobalVars.m_dMeshEps);
			LogMsg("  -mc: Mesh curvature coefficient = %g\n", defGlobalVars.m_dMeshCurvCoeff);
			LogMsg("  -t:  GMRES iteration tolerance = %g\n", defGlobalVars.m_dGmresTol);
			LogMsg("  -d:  Direct potential interaction coefficient to mesh refinement ratio = %g\n", defGlobalVars.m_dEpsRatio);
			LogMsg("  -f:  Out-Of-Core free memory to link memory condition = %g\n", defGlobalVars.m_dOutOfCoreRatio);
			LogMsg("  -g:  Use Galerkin scheme\n");
			LogMsg("  -pj: Use Jacobi Preconditioner\n");
			LogMsg("  -ps: Use two-levels preconditioner with dimension = %d\n", defGlobalVars.m_uiSuperPreDim);
			LogMsg("OPTIONS:\n");
			LogMsg("  -o:  Output refined geometry in FastCap2 format\n");
			LogMsg("  -oi: Dump input geometry in FasterCap format and stop\n");
			LogMsg("  -e:  Output capacitance matrix to file\n");
			LogMsg("  -r:  Dump Gmres residual at each iteration\n");
			LogMsg("  -c:  Dump charge densities in output file\n");
			LogMsg("  -i:  Dump detailed time and memory information\n");
			LogMsg("  -v:  Verbose output\n");
			LogMsg("  -b:  Launch as console/shell application without GUI\n");
			LogMsg("  -b?: Print console usage (this text)\n");
			LogMsg("  -bv: Print only the version\n");
#ifndef __WXMSW__
			LogMsg("  -h: Open only the help system\n");
#endif
		}
	}

	// if return code is an error, print it
	solveMain.PrintRetError(retStatus);

	return retStatus;
}


// returns 'true' if errors in the parsing
bool FasterCapConsole::ParseCmdLine(const char *commandStr, CAutoRefGlobalVars &globalVars, wxString &errMsg)
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

			// '-s' is maximum side left after discretization
			if(argStr[1] == 's') {
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

			// '-c' is dump charge densities in output file
			else if(argStr[1] == 'c') {
                globalVars.m_bOutputCharge = true;
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

			// '-o' is output refined geometry to file (in FastCap2 format),
			// '-oi' is dump input geometry, no refinement, in FastCap2 format and stop
			//       (so actual value of 'm_bOutputGeo' does not matter)
			else if(argStr[1] == 'o') {
				if(argStr[2] == 'i') {
					globalVars.m_bDumpInputGeo = true;
				}
				else {
                    globalVars.m_bOutputGeo = true;
				}
			}

			// '-e' is output capacitance matrix to file
			else if(argStr[1] == 'e') {
				globalVars.m_bOutputCapMtx = true;
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
int FasterCapConsole::getSubstring(const char *buffer, char *substr, int *skip)
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

