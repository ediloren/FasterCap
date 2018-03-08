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


// Solver Global definitions ("local" to the solver sources)

#if !defined(AFX_GLOBAL_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_GLOBAL_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_

#include <string>
#include <list>

// Needed for CLin_Matrix
#include "LinAlgebra/Mtx.h"

// error codes returned by FasterCap. Max value limited to int,
// however to support Linux return codes in a char, limited to 8 bits
// Error codes are mapped when possible to sysexits.h Linux exit codes
#define FC_NORMAL_END			0
// catchall for generic error
#define FC_GENERIC_ERROR		1
// command line usage error
#define FC_COMMAND_LINE_ERROR	64
#define FC_CANNOT_OPEN_FILE		66
#define FC_OUT_OF_MEMORY		71
#define FC_FILE_ERROR			74
#define FC_CANNOT_GO_OOC		97
#define FC_EXCEPTION_ERROR      98
#define FC_USER_BREAK			125


#define G_KILOBYTE	1024
#define G_MEGABYTE	1048576

// type of solver defines, for 'g_ucSolverType'
#define SOLVERGLOBAL_3DSOLVER   1
#define SOLVERGLOBAL_2DSOLVER   2

// constants used in CAutoRefGlobalVars and PotEstimateOpt() preconditioner type
#define AUTOREFINE_PRECOND_NONE				(unsigned char)0
// diagonal preconditioner
#define AUTOREFINE_PRECOND_JACOBI			(unsigned char)1
// block bottom-up preconditioner
#define AUTOREFINE_PRECOND_BLOCK			(unsigned char)2
// compressed matrix top-down preconditioner
#define AUTOREFINE_PRECOND_SUPER			(unsigned char)4
// compressed matrix top-down preconditioner
#define AUTOREFINE_PRECOND_HIER				(unsigned char)8
// value used in discretization pass for signaling to PotEstimateOpt() that should suppress the user warnings
#define AUTOREFINE_DISCRETIZE				(unsigned char)16

#define AUTOREFINE_COLLOCATION				1
#define AUTOREFINE_GALERKIN					2

// states for m_ucHasCmplxPerm
#define AUTOREFINE_REAL_PERM				0
#define AUTOREFINE_CPLX_PERM				1

// maximum size (number of panels) used for preconditioner
//
// super preconditioner
// remark: SOLVE_MAX_SUPER_PRECOND_NUM due to 'unsigned int' limit
// of 'm_puiSupPrecondIndex' dimension in SolveCapacitance.cpp, cannot exceed 65535
#define SOLVE_MAX_SUPER_PRECOND_NUM		4096
#define SOLVE_MIN_SUPER_PRECOND_NUM		16

using namespace std;

// trick to pass comma in the macro argument list
#define DEF_COMMA ,

// defines for memory allocation
//

#define SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM) try { VAR = new TYPE[LEN]; } \
                                catch (bad_alloc&) { return FC_OUT_OF_MEMORY; }    \
                                if( VAR == NULL) { return FC_OUT_OF_MEMORY; }      \
                                MEM += LEN * sizeof(TYPE);

#define SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN) try { VAR = new TYPE[LEN]; } \
                                catch (bad_alloc&) { return FC_OUT_OF_MEMORY; }    \
                                if( VAR == NULL) { return FC_OUT_OF_MEMORY; }

#define SAFENEW_MATRIX_RET(TYPE, VAR, LEN, MEM) try { VAR = new TYPE[LEN][LEN]; } \
                                catch (bad_alloc&) { return FC_OUT_OF_MEMORY; }    \
                                if( VAR == NULL) { return FC_OUT_OF_MEMORY; }      \
                                MEM += LEN * LEN * sizeof(TYPE);

#define SAFENEW_MATRIXNM_RET(TYPE, VAR, LENN, LENM, MEM) try { VAR = new TYPE[LENN][LENM]; } \
                                catch (bad_alloc&) { return FC_OUT_OF_MEMORY; }    \
                                if( VAR == NULL) { return FC_OUT_OF_MEMORY; }      \
                                MEM += LENN * LENM * sizeof(TYPE);

#define SAFENEW_RET(TYPE, VAR, MEM) try { VAR = new TYPE; } \
                                catch (bad_alloc&) { return FC_OUT_OF_MEMORY; }    \
                                if( VAR == NULL) { return FC_OUT_OF_MEMORY; }      \
                                MEM += sizeof(TYPE);

#define SAFENEW_CTOR_RET(TYPE, PARAM, VAR, MEM) try { VAR = new TYPE(PARAM); } \
                                catch (bad_alloc&) { return FC_OUT_OF_MEMORY; }    \
                                if( VAR == NULL) { return FC_OUT_OF_MEMORY; }      \
                                MEM += sizeof(TYPE);

#define SAFENEW_RET_NULL(TYPE, VAR, MEM) try { VAR = new TYPE; } \
                                catch (bad_alloc&) { return NULL; }    \
                                if( VAR == NULL) { return NULL; }      \
                                MEM += sizeof(TYPE);

#ifdef MS_VS
// when using MS VisualC++
#  define isnan _isnan
#  define isfinite _finite
#endif

class CMemoryUsage;

// global, type of solver invoked
extern unsigned char g_ucSolverType;
typedef std::list<std::string> StlStringList;
extern StlStringList g_stlCondNames;
extern CLin_Matrix g_clsCapMatrixRe, g_clsCapMatrixIm;
extern CMemoryUsage g_clsMemUsageCopy, g_clsMemUsage;

class CMemoryUsage
{
public:
	// constructor
	CMemoryUsage();
	void Clear();
	unsigned long GetTotal();
	unsigned long GetTotalKB();

	unsigned long m_ulPanelsMem;
	unsigned long m_ulLinksMem;
	unsigned long m_ulCondMem;
	unsigned long m_ulGmresMem;
	unsigned long m_ulPrecondMem;
	unsigned long m_ulChargesMem;
	unsigned long m_ulHierMem;
};

class CAutoRefGlobalVars
{
public:
	// constructor
	CAutoRefGlobalVars();
	void Reset();

	unsigned char m_ucHasCmplxPerm;
	unsigned char m_ucPrecondType;
	unsigned int m_uiSuperPreDim, m_uiBlockPreSize;
	std::string m_sFileIn;
	double m_dMaxDiscSide, m_dEps, m_dMeshEps, m_dEpsRatio, m_dMeshCurvCoeff, m_dGmresTol;
	double m_dMaxHierPreDiscSide, m_dHierPreEps, m_dHierPreGmresTol;
	double m_dAutoMaxErr, m_dOutOfCoreRatio;
	bool m_bDumpResidual, m_bVerboseOutput, m_bOutputGeo, m_bDumpInputGeo, m_bAuto, m_bAutoPrecond, m_bDumpTimeMem;
	bool m_bKeepCharge, m_bRefineCharge, m_bKeepMesh, m_bOutputCharge, m_bOutputCapMtx;
	char m_cScheme;

	// variables not linked to user options, but to global statuses
	bool m_bWarnGivenPre, m_bWarnGivenThin, m_bWarnGivenSkew, m_bWarnGivenSelfPot;
	bool m_bWarnGivenNaN;
};

#endif //AFX_GLOBAL_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_
