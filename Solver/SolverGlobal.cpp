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


// SolverGlobal.cpp : globals
//
// Enrico Di Lorenzo, 2012/10/03

#include "SolverGlobal.h"

// global, type of solver invoked
unsigned char g_ucSolverType = SOLVERGLOBAL_3DSOLVER;

CAutoRefGlobalVars::CAutoRefGlobalVars()
{
	Reset();
}

void CAutoRefGlobalVars::Reset()
{
	m_sFileIn = "\0";
	// if not specified, precond type is JACOBI
	m_ucPrecondType = AUTOREFINE_PRECOND_JACOBI;
	m_dMaxHierPreDiscSide = 128;
	m_uiSuperPreDim = 128;
	m_uiBlockPreSize = 32;
	m_bDumpResidual = false;
	m_bDumpTimeMem = false;
	m_bVerboseOutput = false;
	m_dGmresTol = 0.01;
	m_dHierPreGmresTol = 0.5;
	m_bOutputGeo = false;
	m_cScheme = AUTOREFINE_COLLOCATION;
	m_bAuto = false;
	m_bAutoPrecond = false;
	m_dAutoMaxErr = 0.01;
	m_dMeshEps = 0.1;
	m_dEpsRatio = 1.0;
	m_dMeshCurvCoeff = 3.0;
	m_dOutOfCoreRatio = 5.0;
	m_bOutputCharge = false;
	// not used any more in 'Run' dialog, but calculated offline and used globally
	m_dMaxDiscSide = 0.1;
	m_dEps = 0.3 * m_dMaxDiscSide;
	// not used any more in 'Run' dialog, but must be kept because
	// trigger operations in the code (watch out: may be used as unsupported options
	// via automation!)
	m_bKeepCharge = false;
	m_bRefineCharge = false;
	m_bKeepMesh = false;
    // not used in the 'Run' dialog, but globally used
	m_ucHasCmplxPerm = AUTOREFINE_REAL_PERM;
	m_bWarnGivenPre = false;
	m_bWarnGivenThin = false;
	m_bWarnGivenSkew = false;
	m_bWarnGivenSelfPot = false;
	m_bWarnGivenNaN = false;
}

CMemoryUsage::CMemoryUsage()
{
	Clear();
}

void CMemoryUsage::Clear()
{
	m_ulPanelsMem = 0;
	m_ulLinksMem = 0;
	m_ulCondMem = 0;
	m_ulGmresMem = 0;
	m_ulPrecondMem = 0;
	m_ulChargesMem = 0;
	m_ulHierMem = 0;
}

unsigned long CMemoryUsage::GetTotal()
{
	return (m_ulPanelsMem +	m_ulChargesMem + m_ulLinksMem + m_ulCondMem + m_ulGmresMem + m_ulPrecondMem + m_ulHierMem);
}

unsigned long CMemoryUsage::GetTotalKB()
{
	return GetTotal() / G_KILOBYTE;
}

