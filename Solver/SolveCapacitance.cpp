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


// CSolveCap.cpp : capacitance solve class
//
// Enrico Di Lorenzo, 2005/10/11

#include "../stdafx.h"

#include <fstream>
#include <time.h>
// for openmp
#include "omp.h"
// for _isnan() and _finite()
#include <float.h>

#include "SolveCapacitance.h"

// link with FasterCap main frame
#include "../FasterCapGlobal.h"

// max number of automatic iterations, after which we declare convergence failed
#define SOLVE_MAX_AUTO_ITERATIONS  128
// minimum increment factor in terms of # of panels / links for each iteration of the 'auto' option
#define SOLVE_AUTO_INCREMENT_FACTOR 1.1
// this is 2 / SIN(60), used to get square of equilateral triangle side, given the area
// (i.e. l^2 = A * 2 / sin(60)
#define SOLVE_SIDE2_FROM_AREA	2.3094


// constructor
CSolveCap::CSolveCap()
{

	// all structures pointed to with these pointers are deallocated
	// by DeallocateMemory()
	m_pclsSupPrecondMtx = NULL;
	m_pclsSupPotMtx = NULL;
	m_pclsTmpCapMtx = NULL;

	m_pCondCharges = NULL;

	m_pclsGmres_q = NULL;
	m_pclsGmres_h = NULL;
	m_pclsGmres_zf = NULL;
	m_pclsGmres1_q = NULL;
	m_pclsGmres1_h = NULL;

	m_clsSupPrecondElements = NULL;
	m_puiSupPrecondIndex = NULL;
	m_pfSupPrecondAreae = NULL;

	m_clsBlockPrecondElements = NULL;
	m_pdBlockPrecond = NULL;
	m_pucBlockPrecondDim = NULL;
}

CSolveCap::~CSolveCap()
{
	// free memory
	DeallocateMemory(AUTOREFINE_DEALLMEM_ALL, m_clsGlobalVars);

}

int CSolveCap::Run(CAutoRefGlobalVars globalVars)
{
	string path, filein;
	double start, finish;
	int ret;

	ret = FC_NORMAL_END;

	LogMsg("Running ");
	LogMsg(FCG_HEADER_VERSION);
	LogMsg("\n");
	LogMsg(FCG_HEADER_COPYRIGHT);
	LogMsg(" ");
	LogMsg(FCG_HEADER_WEBSITE);
	LogMsg("\n");

	LogMsg("Starting capacitance extraction with the following parameters:\n");
	LogMsg("Input file: %s%s\n", path.c_str(), (char*)globalVars.m_sFileIn.c_str());

	// defaults

	// calculate pEps from mesh eps and eps ratio
	globalVars.m_dEps = globalVars.m_dEpsRatio * globalVars.m_dMeshEps;

	if(globalVars.m_dEps == 0.0) {
		ErrMsg("Warning: no tolerance specified for refinement, using standard value 0.01\n");
		globalVars.m_dEps = 0.01;
	}
	if(globalVars.m_dMeshCurvCoeff < 1.0) {
		ErrMsg("Warning: curvature coefficient %f is below the minimum value of 1.0, setting to 1.0\n", globalVars.m_dMeshCurvCoeff);
		globalVars.m_dMeshCurvCoeff = 1.0;
	}
	if(globalVars.m_dOutOfCoreRatio < 0.0) {
		ErrMsg("Warning: Out-of-core free memory to link memory ratio %f is below the minimum value of 0.0, setting to 0.0\n", globalVars.m_dOutOfCoreRatio);
		globalVars.m_dOutOfCoreRatio = 0.0;
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0 && (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_JACOBI) != 0) {
		globalVars.m_ucPrecondType = AUTOREFINE_PRECOND_SUPER;
		ErrMsg("Warning: both two-levels preconditioner and Jacobi preconditioner specified, keeping only the two-levels preconditioner\n");
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0 && globalVars.m_dMaxHierPreDiscSide == 0.0) {
		globalVars.m_dMaxHierPreDiscSide = globalVars.m_dMaxDiscSide;
		ErrMsg("Warning: no preconditioner side discretization value specified, using %f\n", globalVars.m_dMaxHierPreDiscSide);
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0 && globalVars.m_dHierPreEps == 0.0) {
		globalVars.m_dHierPreEps = globalVars.m_dEps * 10.0;
		ErrMsg("Warning: no mutual potential epsilon value specified for hierarchical preconditioner, using %f\n", globalVars.m_dHierPreEps);
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0 && globalVars.m_bRefineCharge == true) {
		globalVars.m_bRefineCharge = false;
		ErrMsg("Warning: 'refine mesh using calculated charges (-s)' option is not compatible with hierarchical preconditioner, resetting option to false\n");
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0 ) {
		if(globalVars.m_uiSuperPreDim > SOLVE_MAX_SUPER_PRECOND_NUM) {
			globalVars.m_uiSuperPreDim = SOLVE_MAX_SUPER_PRECOND_NUM;
			ErrMsg("Warning: maximum two-levels preconditioner dimension is %d, using %d\n", SOLVE_MAX_SUPER_PRECOND_NUM, SOLVE_MAX_SUPER_PRECOND_NUM);
		}
		else if(globalVars.m_uiSuperPreDim < SOLVE_MIN_SUPER_PRECOND_NUM) {
			globalVars.m_uiSuperPreDim = SOLVE_MIN_SUPER_PRECOND_NUM;
			ErrMsg("Warning: minimum two-levels preconditioner dimension is %d, using %d\n", SOLVE_MIN_SUPER_PRECOND_NUM, SOLVE_MIN_SUPER_PRECOND_NUM);
		}
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0 ) {
		if(globalVars.m_uiBlockPreSize > SOLVE_MAX_BLOCK_PRECOND_NUM) {
			globalVars.m_uiBlockPreSize = SOLVE_MAX_BLOCK_PRECOND_NUM;
			ErrMsg("Warning: maximum block preconditioner dimension is %d, using %d\n", SOLVE_MAX_BLOCK_PRECOND_NUM, SOLVE_MAX_BLOCK_PRECOND_NUM);
		}
		else if(globalVars.m_uiBlockPreSize < SOLVE_MIN_BLOCK_PRECOND_NUM) {
			globalVars.m_uiBlockPreSize = SOLVE_MIN_BLOCK_PRECOND_NUM;
			ErrMsg("Warning: minimum block preconditioner dimension is %d, using %d\n", SOLVE_MIN_BLOCK_PRECOND_NUM, SOLVE_MIN_BLOCK_PRECOND_NUM);
		}
	}
	if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
		if( globalVars.m_cScheme == AUTOREFINE_GALERKIN ) {
			globalVars.m_cScheme = AUTOREFINE_COLLOCATION;
			ErrMsg("Warning: Galerkin scheme yet not supported in 2D, moving to Collocation\n");
		}
	}

	// start timer to time overall solve steps
	// remark: using omp_get_wtime() to measure wall time and not processor time,
	// that would provide a wrong information for parallel execution
	start = omp_get_wtime();

	try {
		ret = SolveMain(globalVars, g_clsCapMatrixRe, g_clsCapMatrixIm);
	}
	catch( bad_alloc& ) {
		// all memory exceptions should have been catched during allocations;
		// however, some may have been skipped (for instance stl ones when using MFC, old remark)

		ret = FC_OUT_OF_MEMORY;
	}
// in debug mode, we WANT the exceptions to be thrown
#ifndef _DEBUG
	catch(...) {
		// catch everything, just in case

		ErrMsg("ERROR: unhandled exception, exiting\n");

		ret = FC_EXCEPTION_ERROR;
	}
#endif

// see also wxWidgets "C++ exceptions" chapter in the docs. You have some exception handling capabilities:
// "To summarize, when you use exceptions in your code, you may handle them in the following places, in order of priority:
//  1.  In a try/catch block inside an event handler.
//  2.  In wxApp::OnExceptionInMainLoop().
//  3.  In wxApp::OnUnhandledException()."


	if(ret !=  FC_NORMAL_END) {
		// clean up memory
		//
		g_clsMemUsageCopy = g_clsMemUsage;
		DeallocateMemory(AUTOREFINE_DEALLMEM_ALL, globalVars);
		// was allocated in SolveMain()
		m_clsMulthier.DeallocateMemory();
	}

	finish = omp_get_wtime();

	g_fSolveTime = (float)(finish - start);

	LogMsg("\nTotal allocated memory: %d kilobytes\n", g_clsMemUsageCopy.GetTotalKB());
	LogMsg("Total time: ");
	PrintTime(g_fSolveTime);

	LogMsg("\n");

#ifdef DEBUG_DUMP_OTHER

	LogMsg("Link matrix:\n");
	for(i=0; i<16; i++) {
		for(j=0; j<16; j++) {
			LogMsg("%d ", m_clsSolve.m_clsMulthier.m_iaLinksBtwLevels[i][j]);
		}
		LogMsg("\n");
	}

	LogMsg("\n");

	LogMsg("Potest matrix:\n");
	for(i=0; i<16; i++) {
		for(j=0; j<16; j++) {
			LogMsg("%d ", m_clsSolve.m_clsMulthier.m_iaPotestBtwLevels[i][j]);
		}
		LogMsg("\n");
	}

	LogMsg("\n");

	LogMsg("Links histogram:\n");
	for(i=0; i<31; i++) {
		LogMsg("%d panels with %d links", (int)m_clsSolve.m_clsMulthier.m_laLinkNumStatistics[i], i);
		LogMsg("\n");
	}
	LogMsg("%d panels with more than 30 links\n", (int)m_clsSolve.m_clsMulthier.m_laLinkNumStatistics[i], i);

	LogMsg("\n");

	LogMsg("Average # of links per level:\n");
	for(i=0; i<31; i++) {
		LogMsg("%d panels at level %d with average %d links per panel", (int)m_clsSolve.m_clsMulthier.m_laPanelsPerLevel[i], i, (int) (((float)m_clsSolve.m_clsMulthier.m_laLinksPerLevel[i]) / m_clsSolve.m_clsMulthier.m_laPanelsPerLevel[i]) );
		LogMsg("\n");
	}

#endif

	return ret;
}


int CSolveCap::SolveMain(CAutoRefGlobalVars globalVars, CLin_Matrix &capRe, CLin_Matrix &capIm)
{
	CLin_Matrix	cRe[SOLVE_MAX_AUTO_ITERATIONS], cIm[SOLVE_MAX_AUTO_ITERATIONS];
	double diff;
	long ret;
	int i, j;
	unsigned long oldpanelsnum, oldlinksnum, newpanelsnum, newlinksnum, capnum, numrows;
	clock_t start, finish;
	float solveTime;
	bool not_refined_enough;
	StlAutoCondDeque::iterator itc1;
    
    // if just dumping the input file
    if(globalVars.m_bDumpInputGeo == true) {
        
		// clean memory (if charges / panels were still allocated after last run,
		// and are now not used for discretization)
		m_clsMulthier.Clean(AUTOREFINE_DEALLMEM_AT_START, globalVars);

        ret = m_clsMulthier.ReadFastCapFile(&globalVars);
        // if any error in ReadFastCapFile() (e.g. out of memory, user break, etc.)
        if(ret != FC_NORMAL_END) {
            return ret;
        }
        
        m_clsMulthier.OutputFastCapFile(globalVars.m_sFileIn, "dump");
    }
    // actually solving
    else {
        // init hierarchical multiplication structures (mainly allocating vectors)
        ret = m_clsMulthier.AllocateMemory();
        if(ret != FC_NORMAL_END) {
            return ret;
        }

        if(globalVars.m_bAuto == true) {

            LogMsg("Auto calculation with max error: %g\n", globalVars.m_dAutoMaxErr);
            LogMsg("Remark: Auto option overrides all other Manual settings\n");

            // start timer to time solve step
            start = clock();

            // clarify we start from scratch
            globalVars.m_bRefineCharge = false;
            globalVars.m_bKeepCharge = false;
            globalVars.m_bKeepMesh = false;
            // clean memory (if charges / panels were still allocated after last run,
            // and are now not used for discretization)
            m_clsMulthier.Clean(AUTOREFINE_DEALLMEM_AT_START, globalVars);

            // read input file and build super hierarchy;
            // will also modify 'globalVars' to report if there is complex permittivity
            ret = InputFile(&globalVars);
            if(ret != FC_NORMAL_END) {
                return ret;
            }

            globalVars.m_dGmresTol = globalVars.m_dAutoMaxErr / 2.0;
            // set eps to a value far beyond any threshold for refinement, so at first pass
            // there is no refinement at all
            globalVars.m_dMeshEps = 1E32;
            globalVars.m_dEps = globalVars.m_dEpsRatio * globalVars.m_dMeshEps;

            LogMsg("\n");
            OutputSolveParams(globalVars);
            LogMsg("Number of input panels to solver engine: %lu\n", m_clsMulthier.GetInputPanelNum());

            LogMsg("\nIteration number #0 ***************************\n\n");

            // set global vars for all subsequent calls to SolveCapacitance routines
            m_clsGlobalVars = globalVars;

            LogMsg("***************************************\n");
            LogMsg("Refining the geometry.. \n");
            ret = RefineGeoAndCountLinks();
            if(ret !=  FC_NORMAL_END) {
                return ret;
            }
            LogMsg("Refinement completed\n");
            OutputMeshParams(globalVars);

            ret = SolveComputeLinks();
            if(ret !=  FC_NORMAL_END) {
                return ret;
            }

            AutoSetPrecondType(m_clsMulthier.GetLinksNum(AUTOREFINE_HIER_PRE_0_LEVEL), m_clsMulthier.m_lCondNum, globalVars);
            OutputSolvePrecondType(globalVars);

            ret = SolveForCapacitance(&cRe[0], &cIm[0]);
            if(ret !=  FC_NORMAL_END) {
                return ret;
            }

            OutputCapMtx(cRe[0], cIm[0], globalVars);

            OutputSolveStats(globalVars);

            finish = clock();
            solveTime = (float)(finish - start) / CLOCKS_PER_SEC;

            LogMsg("Iteration time: ");
            PrintTime(solveTime);
            LogMsg("Iteration allocated memory: %d kilobytes\n", g_clsMemUsage.GetTotalKB());

            oldpanelsnum = m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL);
            oldlinksnum = m_clsMulthier.GetLinksNum(AUTOREFINE_HIER_PRE_0_LEVEL);

            // free memory, to be ready for next iteration
            DeallocateMemory(AUTOREFINE_DEALLMEM_AT_END, globalVars);


            for(i=1; i<SOLVE_MAX_AUTO_ITERATIONS; i++) {

                // start timer to time solve step
                start = clock();

                LogMsg("\nIteration number #%d ***************************\n\n", i);

                LogMsg("***************************************\n");
                LogMsg("Increasing the geometric refinement.. \n");
                not_refined_enough = true;
                j = 0;
                do {
                    // make eps converge
                    globalVars.m_dMeshEps = m_clsMulthier.m_dMaxMeshEps / sqrt(2.0);
                    globalVars.m_dEps = globalVars.m_dEpsRatio * globalVars.m_dMeshEps;

                    // clean memory (if charges / panels were still allocated after last run,
                    // and are now not used for discretization)
                    m_clsMulthier.Clean(AUTOREFINE_DEALLMEM_AT_START, globalVars);

                    // read input file and build super hierarchy;
                    // will also modify 'globalVars' to report if there is complex permittivity
                    ret = InputFile(&globalVars);
                    if(ret != FC_NORMAL_END) {
                        return ret;
                    }

                    // set global vars for all subsequent calls to SolveCapacitance routines
                    m_clsGlobalVars = globalVars;

                    ret = RefineGeoAndCountLinks();
                    if(ret !=  FC_NORMAL_END) {
                        return ret;
                    }

                    newpanelsnum = m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL);
                    newlinksnum = m_clsMulthier.GetLinksNum(AUTOREFINE_HIER_PRE_0_LEVEL);

                    // if links or panels number stayed the same, increasing the '-m' parameter had no or small effect, so go on
                    if(newlinksnum > oldlinksnum * SOLVE_AUTO_INCREMENT_FACTOR && newpanelsnum > oldpanelsnum * SOLVE_AUTO_INCREMENT_FACTOR) {
                        not_refined_enough = false;
                    }
                    else {
                        LogMsg("Delta in refined panel and link count w.r.t. previous iteration less than %d%%\n", (int)((SOLVE_AUTO_INCREMENT_FACTOR-1.0) * 100.0));
                        LogMsg("(Panels # %d, Links # %d)\n", newpanelsnum, newlinksnum);
                        if(globalVars.m_bVerboseOutput == true) {
                            LogMsg("(Mesh relative refinement value (-m): %g)\n", globalVars.m_dMeshEps);
                        }
                        LogMsg("Automatically increasing the refinement parameters\n");
                    }

                    j++;

                }
                while(not_refined_enough == true && j<SOLVE_MAX_AUTO_ITERATIONS);
                LogMsg("Refinement completed\n");
                OutputMeshParams(globalVars);

                ret = SolveComputeLinks();
                if(ret !=  FC_NORMAL_END) {
                    return ret;
                }

                AutoSetPrecondType(m_clsMulthier.GetLinksNum(AUTOREFINE_HIER_PRE_0_LEVEL), m_clsMulthier.m_lCondNum, globalVars);
                OutputSolvePrecondType(globalVars);

                ret = SolveForCapacitance(&cRe[i], &cIm[i]);
                if(ret !=  FC_NORMAL_END) {
                    return ret;
                }

                if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM) {
                    diff = FNorm(cRe[i] - cRe[i-1]) / FNorm(cRe[i]);
                }
                else {
                    diff = FNorm(cRe[i] - cRe[i-1], cIm[i] - cIm[i-1]) / FNorm(cRe[i], cIm[i]);
                }

                OutputCapMtx(cRe[i], cIm[i], globalVars);

                LogMsg("Weighted Frobenius norm of the difference between capacitance (auto option): %g\n", diff);
                OutputSolveStats(globalVars);

                finish = clock();
                solveTime = (float)(finish - start) / CLOCKS_PER_SEC;

                LogMsg("Iteration time: ");
                PrintTime(solveTime);
                LogMsg("Iteration allocated memory: %d kilobytes\n", g_clsMemUsage.GetTotalKB());

                if(diff <= globalVars.m_dAutoMaxErr) {
                    capRe = cRe[i];
                    capIm = cIm[i];
                    break;
                }

                // save number of panels and of links of last iteration
                // (only in this case: if the new number of links / panels was not enough, let's refine,
                // but compare with the oldest number of links / panels)
                oldpanelsnum = newpanelsnum;
                oldlinksnum = newlinksnum;


                // free memory, to be ready for next iteration
                DeallocateMemory(AUTOREFINE_DEALLMEM_AT_END, globalVars);
            }
        }
        else {
            // no automatic calculation of parameters
            //

            // clean memory (if charges / panels were still allocated after last run,
            // and are now not used for discretization)
            m_clsMulthier.Clean(AUTOREFINE_DEALLMEM_AT_START, globalVars);

            if(globalVars.m_bRefineCharge == false && globalVars.m_bKeepMesh == false) {
                // read input file and build super hierarchy;
                // will also modify 'globalVars' to report if there is complex permittivity
                ret = InputFile(&globalVars);
                if(ret != FC_NORMAL_END) {
                    return ret;
                }
            }

            LogMsg("\n");
            OutputSolveParams(globalVars);
            OutputMeshParams(globalVars);
            OutputSolvePrecondType(globalVars);
            LogMsg("Number of input panels to solver engine: %lu\n", m_clsMulthier.GetInputPanelNum());
            LogMsg("\n");

            // set global vars for all subsequent calls to SolveCapacitance routines
            m_clsGlobalVars = globalVars;

            LogMsg("***************************************\n");
            LogMsg("Refining the geometry.. \n");
            ret = RefineGeoAndCountLinks();
            if(ret !=  FC_NORMAL_END) {
                return ret;
            }
            LogMsg("Refinement completed\n");

            ret = SolveComputeLinks();
            if(ret !=  FC_NORMAL_END) {
                return ret;
            }

    // debug
    //m_clsMulthier.DebugDumpInteractions();


            ret = SolveForCapacitance(&capRe, &capIm);
            if(ret !=  FC_NORMAL_END) {
                return ret;
            }

            OutputCapMtx(capRe, capIm, globalVars);

            OutputSolveStats(globalVars);
        }
        
        // output the capacitance matrix to file, if requested, in CSV format
        if(globalVars.m_bOutputCapMtx == true) {
            OutputCapMtxToFile(capRe, capIm, globalVars);
        }

        if(globalVars.m_bOutputCharge == true) {
            LogMsg("\n");

            // get number of rows
            numrows = capRe.dim(1);

            // scan all conductors, but do not exceed the number of rows (condition 'capnum < numrows')
            // this is important in 2D case, where the matrix is not square, since the last reference GND
            // conductor is not considered for raising to 1V
            for(itc1 = m_clsMulthier.m_stlConductors.begin(), capnum=0;
                    itc1 != m_clsMulthier.m_stlConductors.end() && capnum < numrows; itc1++) {

                // output charges calculated for each conductor
                if( (*itc1)->m_bIsDiel == false) {

                    LogMsg("Outputting charge densities for conductor %s..\n", (*itc1)->m_sName);
                    // now scan tree for leaves and write out refined FastCap file
                    m_clsMulthier.OutputFastCapFile(globalVars.m_sFileIn, (*itc1)->m_sName, &(m_pCondCharges[capnum]));
                    LogMsg("Done\n");

                    capnum++;
                }
            }
        }

        if(globalVars.m_bOutputGeo == true) {
            LogMsg("\nOutputting refined geometry in FastCap compatible format..\n");
            // now scan tree for leaves and write out refined FastCap file
            m_clsMulthier.OutputFastCapFile(globalVars.m_sFileIn, "ref");
            LogMsg("Done\n");
        }

        CopyCondNames(g_stlCondNames);
    }
    
	g_clsMemUsageCopy = g_clsMemUsage;
	DeallocateMemory(AUTOREFINE_DEALLMEM_AT_END, globalVars);
	m_clsMulthier.DeallocateMemory();

	return ret;
}

void CSolveCap::AutoSetPrecondType(unsigned long numOfLinks, long condNum, CAutoRefGlobalVars &globalVars)
{
	long linksXcond;

	// if no auto setting of the preconditioner, keep the general setting
	if( globalVars.m_bAutoPrecond == false ) {
		return;
	}

	// scale down not to overcome the dimension of the 'long' var
	numOfLinks /= 1000;
	linksXcond = numOfLinks * condNum;

	if( linksXcond <= 500 ) {
		globalVars.m_ucPrecondType = AUTOREFINE_PRECOND_JACOBI;
	}
	else if( (linksXcond > 500) && (linksXcond <= 5000) ) {
		globalVars.m_ucPrecondType = AUTOREFINE_PRECOND_SUPER;
		globalVars.m_uiSuperPreDim = 128;
	}
	else if( (linksXcond > 5000) && (linksXcond <= 10000) ) {
		globalVars.m_ucPrecondType = AUTOREFINE_PRECOND_SUPER;
		globalVars.m_uiSuperPreDim = 512;
	}
	else if( (linksXcond > 10000) ) {
		globalVars.m_ucPrecondType = AUTOREFINE_PRECOND_SUPER;
		globalVars.m_uiSuperPreDim = 1024;
	}
}

void CSolveCap::OutputSolveParams(CAutoRefGlobalVars globalVars)
{
// debug
//	LogMsg("Eps (-e): %g, max discr side (-s): %g:\n", globalVars.m_dEps, globalVars.m_dMaxDiscSide);
//	LogMsg("Eps (-e): %g, mesh refinement (-m): %g, mesh curvature (-mc): %g\n", globalVars.m_dEps, globalVars.m_dMeshEps, globalVars.m_dMeshCurvCoeff);

    if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
		LogMsg("3D Solver Engine invoked\n");
	}
	else {
		LogMsg("2D Solver Engine invoked\n");
	}

	if(globalVars.m_bOutputGeo == true) {
		LogMsg("Output refined geometry in FastCap compatible format (-o)\n");
	}
	if(globalVars.m_bOutputCharge == true) {
		LogMsg("Output charge densities information (-c)\n");
	}
	if(globalVars.m_bOutputCapMtx == true) {
		LogMsg("Output capacitance matrix to file (-e)\n");
	}

	LogMsg("Solution scheme (-g): ");
	switch(globalVars.m_cScheme) {
	case AUTOREFINE_GALERKIN:
		LogMsg("Galerkin, ");
		break;
	default:
		LogMsg("Collocation, ");
	}

	LogMsg("GMRES tolerance (-t): %g\n", globalVars.m_dGmresTol);

	if(globalVars.m_bRefineCharge == true) {
		LogMsg("Refine mesh using charges (-s)\n");
	}
	if(globalVars.m_bKeepCharge == true) {
		LogMsg("Keep charge information (-kc)\n");
	}
	if(globalVars.m_bKeepMesh == true) {
		LogMsg("Keep start mesh information from previous run (-km)\n");
	}
	LogMsg("Out-of-core free memory to link memory condition (-f): %g\n", globalVars.m_dOutOfCoreRatio);
	LogMsg("Potential interaction coefficient to mesh refinement ratio (-d): %g\n", globalVars.m_dEpsRatio);
	LogMsg("Mesh curvature (-mc): %g\n", globalVars.m_dMeshCurvCoeff);
}

void CSolveCap::OutputMeshParams(CAutoRefGlobalVars globalVars)
{
	LogMsg("Mesh refinement (-m): %g\n", globalVars.m_dMeshEps);
}

void CSolveCap::OutputSolvePrecondType(CAutoRefGlobalVars globalVars)
{
	LogMsg("Precond Type(s) (-p): ");
	if( globalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE ) {
		LogMsg("None");
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_JACOBI) != 0 ) {
		LogMsg("Jacobi ");
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0 ) {
		LogMsg("Block, block preconditioner dimension (-pb): %d ", globalVars.m_uiBlockPreSize);
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0 ) {
		LogMsg("Two-levels, two-levels preconditioner dimension (-ps): %d ", globalVars.m_uiSuperPreDim);
	}
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0 ) {
		LogMsg("Hierarchical, Eps (-phe): %g, max discr side (-phs): %g, tol (-pht): %g ", globalVars.m_dHierPreEps, globalVars.m_dMaxHierPreDiscSide, globalVars.m_dHierPreGmresTol);
	}
	LogMsg("\n");
}

void CSolveCap::OutputSolveStats(CAutoRefGlobalVars globalVars)
{
	g_lPanelsNum = m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL);
	// links # do not consider auto potential calculations. Must add to have the overall number of links.
	g_lLinksNum = m_clsMulthier.GetTotalLinksNum(AUTOREFINE_HIER_PRE_0_LEVEL);

	LogMsg("\nSolve statistics:\n");
	LogMsg("Number of input panels: %lu of which %lu conductors and %lu dielectric\n", m_clsMulthier.m_ulInputRawCondPanelNum+m_clsMulthier.m_ulInputRawDielPanelNum, m_clsMulthier.m_ulInputRawCondPanelNum, m_clsMulthier.m_ulInputRawDielPanelNum);
	LogMsg("Number of input panels to solver engine: %lu\n", m_clsMulthier.m_ulInputPanelNum);
	LogMsg("Number of panels after refinement: %lu\n", m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL));
	LogMsg("Number of potential estimates: %lu\n", m_clsMulthier.m_ulNumofpotest);
	LogMsg("Number of links: %lu (uncompressed %lu, compression ratio is %.1f%%)\n", g_lLinksNum, g_lPanelsNum * g_lPanelsNum, 100.0f - ((float)(g_lLinksNum)) / ((float)(g_lPanelsNum * g_lPanelsNum)) * 100.0);;
	if( (globalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0 ) {
		LogMsg("Number of precond panels: %lu\n", m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_1_LEVEL));
		LogMsg("Number of precond links: %lu\n", m_clsMulthier.GetLinksNum(AUTOREFINE_HIER_PRE_1_LEVEL));
	}
	LogMsg("Max recursion level: %d\n", m_clsMulthier.m_iMaxLevel);
	LogMsg("Max Mesh relative refinement value: %g\n", m_clsMulthier.m_dMaxMeshEps);


	if(m_clsGlobalVars.m_bDumpTimeMem == true) {

		LogMsg("Time for reading input file: %fs\n", m_clsMulthier.m_fDurationReadFile);
		LogMsg("Time for building super hierarchy: %fs\n", m_clsMulthier.m_fDurationSuperH);
		LogMsg("Time for discretization: %fs\n", m_clsMulthier.m_fDurationDiscretize);
		LogMsg("Time for building potential matrix: %fs\n", m_clsMulthier.m_fDurationRefine);
		LogMsg("Time for precond calculation: %fs\n", m_fDurationPrecond);
		LogMsg("Time for gmres solving: %fs\n", m_fDurationSolve);

		LogMsg("Memory allocated for panel hierarchy: %lu bytes\n", g_clsMemUsage.m_ulPanelsMem);
		LogMsg("Memory allocated for links structure: %lu bytes\n", g_clsMemUsage.m_ulLinksMem);
		LogMsg("Memory allocated for conductor list: %lu bytes\n", g_clsMemUsage.m_ulCondMem);
		LogMsg("Memory allocated for Gmres: %lu bytes\n", g_clsMemUsage.m_ulGmresMem);
		LogMsg("Memory allocated for preconditioner: %lu bytes\n", g_clsMemUsage.m_ulPrecondMem);
		LogMsg("Memory allocated for charge vectors: %lu bytes\n", g_clsMemUsage.m_ulChargesMem);
		LogMsg("Memory allocated for hierarchical multiplication: %lu bytes\n", g_clsMemUsage.m_ulHierMem);
	}
}

// output capacitance matrix
int CSolveCap::OutputCapMtx(const CLin_Matrix &matrixRe, const CLin_Matrix &matrixIm, CAutoRefGlobalVars globalVars)
{
	CLin_subscript M, N, i, j;
	long ret;
	StlStringList stringList;
	StlStringList::iterator its;
	double sum;

	LogMsg("Capacitance matrix is:\n");

	CopyCondNames(stringList);

	M = matrixRe.num_rows();
	N = matrixRe.num_cols();

	// verify that matrix dimension is the same for the real and imaginary parts
	ASSERT(M == matrixIm.num_rows());
	ASSERT(N == matrixIm.num_cols());

	// if 2D solver, remove last row and column, unless verbose
	if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
		// verify that matrix dimension is compatible with the size of the conductor names array
		ASSERT(M == stringList.size()-1);
		if(globalVars.m_bVerboseOutput == false) {
			if(N>1) {
				N--;
			}
			else {
				ErrMsg("Error: only one input conductor defined for a 2D problem. Result is meaningless.");
			}
		}
	}
	else {
		// verify that matrix dimension is compatible with the size of the conductor names array
		ASSERT(M == stringList.size());
	}

	ret = LogMsg("Dimension %d x %d\n", M, N);

	for (i=0, its = stringList.begin(); i<M && its != stringList.end(); i++, its++) {
		// conductor name
		LogMsg("%s  ", its->c_str());
		for (j=0; j<N; j++) {

			if( globalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
				if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
					ret += LogMsg("%g ", matrixRe[i][j]);
				}
				else {
					ret += LogMsg("%g ", matrixRe[i][j] * TWO_PI_TIMES_E0);
				}
			}
			else {
				if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
					ret += LogMsg("%g%+gj ", matrixRe[i][j], matrixIm[i][j]);
				}
				else {
					ret += LogMsg("%g%+gj ", matrixRe[i][j] * TWO_PI_TIMES_E0, matrixIm[i][j] * TWO_PI_TIMES_E0);
				}
			}
		}
		ret += LogMsg("\n");
	}

	// issue warnings, if this is the case
	//

	// check for non-negative off-diagonals
	for (i=0; i<M; i++) {
		for (j=0; j<N; j++) {
			if( i != j && matrixRe[i][j] >= 0.0 ) {
				ErrMsg("Warning: capacitance matrix has a non-negative off-diagonal element at row %d col %d\n", i+1, j+1);
			}
		}
	}

	// check diagonal dominance
	// (remark: the assumption is that for a real capacitance matrix, the imaginary part is all zeros)
	for (i=0; i<M; i++) {
		sum = 0.0;
		for (j=0; j<N; j++) {
			if( i != j ) {
				sum += sqrt(matrixRe[i][j]*matrixRe[i][j]+matrixIm[i][j]*matrixIm[i][j]);
			}
		}
		if( sqrt(matrixRe[i][i]*matrixRe[i][i]+matrixIm[i][i]*matrixIm[i][i]) < sum ) {
			ErrMsg("Warning: capacitance matrix is not diagonally dominant due to row %d\n", i+1);
		}
	}

	LogMsg("\n");
   
	return ret;
}

void CSolveCap::CopyCondNames(StlStringList &stringList)
{
	StlAutoCondDeque::iterator itc1;

	// delete string list
	stringList.clear();

	for(itc1 = m_clsMulthier.m_stlConductors.begin(); itc1 != m_clsMulthier.m_stlConductors.end(); itc1++) {
		// if conductor is real conductor and not dielectric
		if( (*itc1)->m_bIsDiel == false) {
			stringList.push_back((*itc1)->m_sName);
		}
	}
}

// output capacitance matrix to file
int CSolveCap::OutputCapMtxToFile(const CLin_Matrix &matrixRe, const CLin_Matrix &matrixIm, CAutoRefGlobalVars globalVars)
{
	CLin_subscript M, N, i, j;
	FILE *fout;
	std::string basefilename, fileoutname;
	std::string::size_type pos;


    M = matrixRe.num_rows();
    N = matrixRe.num_cols();

    // verify that matrix dimension is the same for the real and imaginary parts
    ASSERT(M == matrixIm.num_rows());
    ASSERT(N == matrixIm.num_cols());

    // if 2D solver, remove last row and column
    if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
        if(N>1) {
            N--;
        }
        else {
            ErrMsg("Error: only one input conductor defined for a 2D problem. Result is meaningless.");
        }
    }

    // build file out name from 'globalVars.m_sFileIn', stripping old extension if any
    // and adding ".csv"

    basefilename = globalVars.m_sFileIn;
    pos = basefilename.rfind(".");
    // if there was an extension, strip it, otherwise nothing
    if(pos != basefilename.npos ) {
        basefilename.resize(pos);
    }
    fileoutname = basefilename;
    fileoutname += ".csv";

    fout = fopen(fileoutname.c_str(), "w");

    if(fout == NULL) {
        ErrMsg("Warning: cannot open \"%s\" for writing, skipping\n", fileoutname.c_str());
        return FC_FILE_ERROR;
    }

    LogMsg("Generating capacitance matrix file \"%s\"\n", fileoutname.c_str());
   
    for (i=0; i<M; i++) {
         for (j=0; j<N; j++) {

            if( globalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
                if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
                    fprintf(fout, "%g", matrixRe[i][j]);
                }
                else {
                    fprintf(fout, "%g", matrixRe[i][j] * TWO_PI_TIMES_E0);
                }
            }
            else {
                if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
                    fprintf(fout, "%g%+gj", matrixRe[i][j], matrixIm[i][j]);
                }
                else {
                    fprintf(fout, "%g%+gj", matrixRe[i][j] * TWO_PI_TIMES_E0, matrixIm[i][j] * TWO_PI_TIMES_E0);
                }
            }
            
            if(j<N-1) {
                fprintf(fout, ", ");
            }
        }
        fprintf(fout, "\n");
    }
    
    fclose(fout);
            
	return FC_NORMAL_END;
}

void CSolveCap::PrintTime(float solveTime)
{
	long days, hours, minutes, seconds;

	// as per ANSI C spec, if both integer operands are positive or unsigned,
	// the result is truncated toward 0.
	seconds = long(solveTime);
	minutes = seconds / 60;
	hours = minutes / 60;
	days = hours / 24;
	// so now calculate the reminders
	seconds -= minutes * 60;
	minutes -= hours * 60;
	hours -= days * 24;

	LogMsg("%fs (%d days, %d hours, %d mins, %d s)\n", solveTime, days, hours, minutes, seconds);
}


void CSolveCap::PrintRetError(int retErr)
{

	if(retErr == FC_GENERIC_ERROR) {
		ErrMsg("\n\nError: Generic error, program execution stopped!\n");
	}
	else if(retErr == FC_COMMAND_LINE_ERROR) {
		ErrMsg("\n\nError: Command line error, solver not launched\n");
	}
	else if(retErr == FC_CANNOT_OPEN_FILE) {
		ErrMsg("\n\nError: Cannot open input file, program execution stopped!\n");
	}
	else if(retErr == FC_OUT_OF_MEMORY) {
		ErrMsg("\n\nError: Out of memory, program execution stopped!\n");
	}
	else if(retErr == FC_FILE_ERROR) {
		ErrMsg("\n\nError: File error, program execution stopped!\n");
	}
	else if(retErr == FC_CANNOT_GO_OOC) {
		ErrMsg("\n\nError: Cannot go Out-of-Core (possibly disk full)!\n");
	}
	else if(retErr == FC_EXCEPTION_ERROR) {
		ErrMsg("\n\nError: Unknown exception catched!\n");
	}
	else if(retErr == FC_USER_BREAK) {
		ErrMsg("\n\nWarning: Program execution stopped on user request!\n");
	}
	else if(retErr != FC_NORMAL_END) {
		ErrMsg("\n\nError: Unknown error code, program execution stopped!\n");
	}

	ErrMsg("\n");
}

int CSolveCap::InputFile(CAutoRefGlobalVars *globalVars)
{
	int ret;

	ret = m_clsMulthier.ReadFastCapFile(globalVars);
	// if any error in ReadFastCapFile() (e.g. out of memory, user break, etc.)
	if(ret != FC_NORMAL_END) {
		return ret;
	}

	// build panels superhierarchy
	ret = m_clsMulthier.BuildSuperHierarchy();
	// if any error in BuildSuperHierarchy() (e.g. out of memory, user break, etc.)
	if(ret != FC_NORMAL_END) {
		return ret;
	}

	return ret;
}

int CSolveCap::RefineGeoAndCountLinks()
{
	int ret;
	CAutoRefGlobalVars globalVars;

	m_clsMulthier.SetInteractionLevel(AUTOREFINE_HIER_PRE_0_LEVEL);

	if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0) {

		// first step of refinement to form the hierarchical super pre-conditioner
		//

		globalVars.m_dEps = m_clsGlobalVars.m_dHierPreEps;
		globalVars.m_dMaxDiscSide = m_clsGlobalVars.m_dMaxHierPreDiscSide;
		globalVars.m_bOutputGeo = false;
		globalVars.m_bOutputCapMtx = false;

		ret = m_clsMulthier.AutoRefinePanels(globalVars, AUTOREFINE_HIER_PRE_1_LEVEL);
		ret = m_clsMulthier.AutoRefineLinks(globalVars);

		if(ret !=  FC_NORMAL_END) {
			return ret;
		}

		// second step of refinement for actual solution
		//

		globalVars = m_clsGlobalVars;
		ret = m_clsMulthier.AutoRefinePanels(globalVars, AUTOREFINE_HIER_PRE_0_LEVEL);
	}
	else {
		// refinement for actual solution
		globalVars = m_clsGlobalVars;
		ret = m_clsMulthier.AutoRefinePanels(globalVars, AUTOREFINE_HIER_PRE_0_LEVEL);
	}

	return ret;
}

int CSolveCap::SolveComputeLinks()
{
	int ret;

	LogMsg("***************************************\n");
	LogMsg("Computing the links.. \n");
	LogMsg("Number of panels after refinement: %lu\n", m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL));
	LogMsg("Number of links to be computed: %lu\n", m_clsMulthier.GetLinksNum(AUTOREFINE_HIER_PRE_0_LEVEL));

	m_clsMulthier.SetInteractionLevel(AUTOREFINE_HIER_PRE_0_LEVEL);
	ret = m_clsMulthier.AutoRefineLinks(m_clsGlobalVars);

	LogMsg("Done computing links\n");
	LogMsg("***************************************\n");

	return ret;
}

int CSolveCap::SolveForCapacitance(CLin_Matrix *cRe, CLin_Matrix *cIm)
{
	int ret;

	// allocate memory
	ret = AllocateMemory();
	if(ret != FC_NORMAL_END) {
		return ret;
	}

	MakePreconditioner();

	ret = Solve(cRe, cIm);

	return ret;
}

// Make the preconditioner
void CSolveCap::MakePreconditioner()
{
	StlAutoCondDeque::iterator itc1, itc2;
	clock_t start, finish;
	unsigned long roots;

	// start timer to time preconditioner computation
	start = clock();

	if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0) {

		// compute super preconditioner variables
		m_uiSupPreNum = 0;

		roots = m_clsMulthier.m_lCondNum + m_clsMulthier.m_lDielNum;
		// if more roots than super precond size
		if( roots >= m_clsGlobalVars.m_uiSuperPreDim ) {
			// mark direct super preconditioner as not available - must rely only on what we have
			m_uiSupPreLevel = 0;
			LogMsg("Number of conductors and dielectric interfaces greater than %d; cannot use two-levels preconditioner\r\n", m_clsGlobalVars.m_uiSuperPreDim);
			m_clsGlobalVars.m_ucPrecondType &= ~(AUTOREFINE_PRECOND_SUPER);
		}
		else {
			// calculate how deep we can dive into each conductor's tree for calculating
			// the preconditioner, given the maximum dimension of preconditioner matrix
			// and the number of conductors
			m_uiSupPreLevel = (unsigned int) (log( (double)(m_clsGlobalVars.m_uiSuperPreDim / roots) ) / log(2.0));
		}
	}

	// recursively compute Super precond, Block precond (if needed)
	if(m_clsGlobalVars.m_ucPrecondType != AUTOREFINE_PRECOND_NONE) {
		m_ulPanelNum = 0;
		m_ulBlockPreBaseNum = 0;
		for(itc1 = m_clsMulthier.m_stlConductors.begin(); itc1 != m_clsMulthier.m_stlConductors.end(); itc1++) {

			m_iLevel = -1;
			// only if using Block precond
			if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0 ) {
				m_bIsComputingBlock = false;
			}
			m_pCurrCond = *itc1;
			RecurseComputePrecond((*itc1)->m_uTopElement.m_pTopElement);
		}
	}

	// only if using Super precond or Hier precond and Super precond can be used to accelerate iterations
	if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0) {

		// and finally compute the super preconditioner (inverting the matrix)
		ComputeSuperPrecond();
	}

	finish = clock();
	m_fDurationPrecond = (float)(finish - start) / CLOCKS_PER_SEC;



#ifdef DEBUG_DUMP_POT
	DebugDumpPotMtx();
#endif //DEBUG_DUMP_POT
#ifdef DEBUG_DUMP_UNCOMP_POT
	DebugDumpUncomprPotMtx();
#endif //DEBUG_DUMP_UNCOMP_POT
#ifdef DEBUG_DUMP_BASIC
#  ifdef DEBUG_DUMP_OTHER
	// debug only matrix dump
//	DebugDumpPrecondMtx();
////	m_clsMulthier.DebugOutputTree();
////	m_clsMulthier.DebugOutputTree2();
#  endif
#endif
}


// Remark: it is not necessary to allocate
// solution matrix 'cRe' and 'cIm' beforehand
int CSolveCap::Solve(CLin_Matrix *cRe, CLin_Matrix *cIm)
{
	StlAutoCondDeque::iterator itc1, itc2;
	unsigned long i, j, potindex, chindex, caprow, capcol, panelnum, potVectorDim;
	long firstCondElemIndex, rowNum;
	double capacitance[2], outpermRe, outpermIm, chargeRe, chargeIm;
	double start, finish;
	int ret;

	// init 'ret' status, in case anything fails
	ret = FC_GENERIC_ERROR;

	panelnum = m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL);

	// allocate vectors
	if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
		potVectorDim = panelnum;
	}
	else {
		potVectorDim = panelnum * 2;
	}
	CLin_Vector potential(potVectorDim, 0.0);
	CLin_Vector charge(potVectorDim);

	// if 2D solver, remove last row (no need, and saves time)
	rowNum = m_clsMulthier.m_lCondNum;
	if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
		if(m_clsMulthier.m_lCondNum>1) {
			rowNum--;
		}
		else {
			ErrMsg("Error: only one input conductor defined for a 2D problem. Result is meaningless.");
		}
	}
	// allocate solution matrix
	cRe->newsize(rowNum, m_clsMulthier.m_lCondNum);
	// init and fill up with zeros, in case the matrix is used (e.g. norm calculation etc.)
	cIm->newsize(rowNum, m_clsMulthier.m_lCondNum);
	*cIm = 0.0;


//-firststep-	CLin_Vector X0(m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL));

	// start timer to time gmres solve
	// remark: using omp_get_wtime() to measure wall time and not processor time,
	// that would provide a wrong information for parallel execution
	start = omp_get_wtime();

	// init variable to alternate precond schemes
	m_ucAlternatePrecond = 0;

	// init leaf panel counter for potential computations
	potindex = 0;

	// init index of first conductor element
	firstCondElemIndex = -1;

	// scan every conductor to find the charge when the current
	// conductor is raised to unit potential while the others
	// are kept at zero
	for(itc1 = m_clsMulthier.m_stlConductors.begin(), caprow=0;
	        itc1 != m_clsMulthier.m_stlConductors.end(); itc1++) {

		// if conductor is real conductor and not dielectric
		if( (*itc1)->m_bIsDiel == false) {

			if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {

				// raise conductor's potential to one
				// remark: also for complex potentials, voltages are 1.0+j0.0,
				// so no need to change the value of the Im part
				for(i=potindex; i < potindex + (*itc1)->m_ulLeafPanelNum; i++) {
					potential[i] = 1.0;
				}
			}
			else if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
				// if 2D, must correct for the 'k' integration constant;
				// we use the work-around from A. R. Djordjevic, R. F. Harrington, T. K. Sarkar,
				// "Evaluation of quasi-static matrix parameters for multiconductor transmission
				// lines using Galerkin's method", IEEE Transactions on Microwave Theory and Techniques,
				// Vol. 42, No. 7, 1994

				// So subtract last element of the potential vector from all the other elements.
				// However the last element is always zero unless we are raising the last conductor to 1.0 potential.
				// For complex potentials, voltages are 1.0+j0.0, so imaginary part does not make a difference.
				// Remark:  if dielectrics are presents, the 'work-around' must consider the last row
				// *belonging to a conductor*, not the last row in absolute (since this could belong to a dielectric)

				// mark the index of the first element belonging to the first conductor in the list
				// (by construction of the list, first there are all the dielectrics, then all the conductors)
				if(firstCondElemIndex == -1) {
					firstCondElemIndex = (long)potindex;
					m_clsMulthier.m_ulFirstCondElemIndex = potindex;
				}

				if(potindex + (*itc1)->m_ulLeafPanelNum == panelnum) {
					// no need to calculate the last row, that is simply the forced balance
					// (results always equal to the algebraic sum of the previous rows)
					continue;
					//for(i=(unsigned long)firstCondElemIndex; i < potindex; i++) {
					//	potential[i] = -1.0;
					//}
				}
				else {
					// raise conductor's potential to one
					// remark: also for complex potentials, voltages are 1.0+j0.0,
					// so no need to change the value of the Im part
					for(i=potindex; i < potindex + (*itc1)->m_ulLeafPanelNum; i++) {
						potential[i] = 1.0;
					}
				}
			}
			else {
				ASSERT(false);
			}

			// solve for charge

			// if block precond and a second precond method, we use alternate preconditioning,
			// so we need the flex gmres version
			if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0  &&
			        (m_clsGlobalVars.m_ucPrecondType & ~(AUTOREFINE_PRECOND_BLOCK)) != 0) {

				// this is not supported at the moment, also considering complex capacitance
				ErrMsg("Error: using unsupported preconditioner 'b'!");

				ret = gmresFlexPrecondSFastAll(&potential, &charge, m_clsGlobalVars.m_dGmresTol);

				if(ret !=  FC_NORMAL_END) {
					return ret;
				}
			}
			// if hierarchical precond, need flex gmres version (upper step, being approximated
			// through a second gmres iteration, gives different results at each main gmres iteration)
			else if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0) {

				// this is not supported at the moment, also considering complex capacitance
				ErrMsg("Error: using unsupported preconditioner 'h'!");

				ret = gmresFlexPrecondSFastAll(&potential, &charge, m_clsGlobalVars.m_dGmresTol);

				if(ret !=  FC_NORMAL_END) {
					return ret;
				}
			}
			else {
				ret = gmresPrecondSFastAll(&potential, &charge, m_clsGlobalVars.m_dGmresTol);

				if(ret !=  FC_NORMAL_END) {
					return ret;
				}
			}

// debug
//for(j=0; j < panelnum; j++) {
//	chargeRe = charge[j];
//	chargeIm = charge[j+panelnum];
//	LogMsg("Index %d, chargeRe %g, chargeIm %g\n", j, chargeRe, chargeIm);
//}

			//
			// sum up charge to find conductor's capacitance
			//

			// init leaf panel counter for potential computations
			chindex = 0;

			for(itc2 = m_clsMulthier.m_stlConductors.begin(), capcol=0;
			        itc2 != m_clsMulthier.m_stlConductors.end(); itc2++) {

				// sum up only if current conductor is real conductor and not dielectric
				if( (*itc2)->m_bIsDiel == false) {

					capacitance[0] = 0.0;
					capacitance[1] = 0.0;

					// sum this conductor's charges
					//
					// note that charge for this conductor is multiplied by outperm;
					// this is to consider the medium surrounding the conductor,
					// for detailed explanation, see Keith Nabors's PhD Thesis,
					// "Efficient three-dimensional capacitance calculation", p.76
					// as well as "The electrostatic field of conductive bodies in
					// multiple dielectric media", Rao & Sarkar; another good reference is
					// C. Wei, R. F. Harrington, J. R. Mautz, T. K. Sarkar
					// "Multiconductor transmission lines in multilayered dielectric media", IEEE Transactions on Microwave
					// Theory and Techniques, Vol. 32, No. 4, Apr 1985
					for(j=chindex; j < chindex + (*itc2)->m_ulLeafPanelNum; j++) {

						if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
							_ASSERT(!isnan(charge[j]));
							_ASSERT(isfinite(charge[j]));
							// get the outer permittivity for this panel
							outpermRe = (*itc2)->m_dSurfOutperm[m_clsMulthier.m_pucDielIndex[j]][0];
							// and adjust the charge
							capacitance[0] += outpermRe*charge[j];
						}
						else {
							// real part
							_ASSERT(!isnan(charge[j]));
							_ASSERT(isfinite(charge[j]));
							// imaginary part
							_ASSERT(!isnan(charge[j+panelnum]));
							_ASSERT(isfinite(charge[j+panelnum]));
							// get the complex outer permittivity for this panel
							outpermRe = (*itc2)->m_dSurfOutperm[m_clsMulthier.m_pucDielIndex[j]][0];
							outpermIm = (*itc2)->m_dSurfOutperm[m_clsMulthier.m_pucDielIndex[j]][1];

							// and adjust the charge
							// remark: a*b where a and b are complex is (aRe+j*aIm)*(bRe+j*bIm) = aRe*bRe-aIm*bIm +j*(aIm*bRe+aRe*bIm)
							chargeRe = charge[j];
							chargeIm = charge[j+panelnum];
							capacitance[0] += outpermRe*chargeRe - outpermIm*chargeIm;
							capacitance[1] += outpermIm*chargeRe + outpermRe*chargeIm;
//debug
//LogMsg("outpemRe %g, outpermIm %g, chargeRe %g, chargeIm %g; c[0] %g, c[1] %g\n", outpermRe, outpermIm, chargeRe, chargeIm, capacitance[0], capacitance[1]);
						}

					}
					// store capacitance value in matrix
					(*cRe)[caprow][capcol] = capacitance[0];
					(*cIm)[caprow][capcol] = capacitance[1];
					// increase capacitance matrix col index
					capcol++;
				}

				// increment chindex
				chindex += (*itc2)->m_ulLeafPanelNum;
			}

			// store charges for later use, if requested
			// Remark: charges from the solve on each conductor are summed up;
			// this is a simplification, in order not to have to save the charges
			// for each conductor and use it to discretize differently the geometry
			// for each solve pass (would not be good for speed improvement)
			if(m_clsGlobalVars.m_bKeepCharge == true) {
				// charges
				// remark: you must modify also the functions CAutoRefine::CopyCharges()
				// and CAutoRefine::OutputPanelTree(), dividing charge by the panel area,
				// to work with charge densities
				for(i=0; i < panelnum; i++) {
//					m_clsMulthier.m_fGlobalCharges[i] += fabs(charge[i]);
					m_clsMulthier.m_fGlobalCharges[i] += charge[i];
				}
			}

			if(m_clsGlobalVars.m_bOutputCharge == true) {
				// charges
				// remark: the using functions (e.g. CAutoRefine::CopyCharges()
				// and CAutoRefine::OutputPanelTree() ), will divide the charge by the panel area,
				// to work with charge densities
				//
				// update vector dimension
				m_pCondCharges[caprow].newsize(panelnum);
				// and copy charges into the vector
				// remark: this is only the real part, in case of complex charges
				for(i=0; i < panelnum; i++) {
					if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
						m_pCondCharges[caprow][i] = charge[i];
					}
					else {
						m_pCondCharges[caprow][i] = charge[i] * TWO_PI_TIMES_E0;
					}
				}
			}

			// clean up
			//

			// reset conductor's potential to zero
			for(i=potindex; i < potindex + (*itc1)->m_ulLeafPanelNum; i++) {
				potential[i] = 0.0;
			}

			// increase capacitance matrix row index
			caprow++;

		}

		// increase counter
		potindex += (*itc1)->m_ulLeafPanelNum;
	}

	finish = omp_get_wtime();
	m_fDurationSolve = (float)(finish - start);

	return ret;
}

// pre-allocate gmres vectors
int CSolveCap::AllocateMemory()
{
	unsigned long numElems_0, numPanels_0;
	bool ret;

	//
	// declare large arrays / matrices on the heap, not on the stack (otherwise will run out of stack
	// at run time)
	//

	numPanels_0 = m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL);
	if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
		numElems_0 = numPanels_0;
	}
	else {
		numElems_0 = numPanels_0 * 2;
	}

	// array of vectors to store charges, in case charge densities dump is required
	// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
	SAFENEW_ARRAY_NOMEM_RET(CLin_Vector, m_pCondCharges, m_clsMulthier.m_lCondNum)

	// preallocate vectors and matrices for gmres on bottom level

	//
	// allocate matrices
	// (note that we use 0-index arrays)
	//

	// q and h matrices
	if(m_pclsGmres_q == NULL) {
		// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
		SAFENEW_ARRAY_NOMEM_RET(CLin_Vector, m_pclsGmres_q, SOLVE_GMRES_ITER_MAX+1)
	}
	if(m_pclsGmres_q == NULL) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += (SOLVE_GMRES_ITER_MAX+1) * sizeof(CLin_Vector);
	// init first and second element
	ret = m_pclsGmres_q[0].newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
	ret = m_pclsGmres_q[1].newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);

	if(m_pclsGmres_h == NULL) {
		// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
		SAFENEW_ARRAY_NOMEM_RET(CLin_Vector, m_pclsGmres_h, SOLVE_GMRES_ITER_MAX)
	}
	if(m_pclsGmres_h == NULL) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += SOLVE_GMRES_ITER_MAX * sizeof(CLin_Vector);
	// init first element
	ret = m_pclsGmres_h[0].newsize(2);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += 2 * sizeof(double);

	// zf matrix (used in flex variant of gmres)
	if(m_pclsGmres_zf == NULL) {
		// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
		SAFENEW_ARRAY_NOMEM_RET(CLin_Vector, m_pclsGmres_zf, SOLVE_GMRES_ITER_MAX)
	}
	if(m_pclsGmres_zf == NULL) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += SOLVE_GMRES_ITER_MAX * sizeof(CLin_Vector);
	// init first element
	ret = m_pclsGmres_zf[0].newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);

	//
	// declare vectors (note that we use
	// 0-index arrays)
	//
	// g, v, w, z, y vectors
	ret = m_clsGmres_g.newsize(numElems_0+1);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += (numElems_0 + 1) * sizeof(double);
	ret = m_clsGmres_v.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
	ret = m_clsGmres_w.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
	ret = m_clsGmres_z.newsize(numElems_0+1);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += (numElems_0 + 1) * sizeof(double);
	ret = m_clsGmres_y.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
	// c,s vectors
	ret = m_clsGmres_c.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
	ret = m_clsGmres_s.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
	// r vector
	ret = m_clsGmres_r.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);

	// Pq vector (preconditioned vector)
	ret = m_clsGmres_Pq.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
	// x0 vector (initial solution vector)
	ret = m_clsGmres_x0.newsize(numElems_0);
	if(ret == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);

	// record up to which iteration the arrays have been pre-allocated
	m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_0_LEVEL] = 0;

	// only in case the user is requesting a hierarchical preconditioner,
	// pre-allocate matrices for gmres on upper levels of hierarchy
	if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_HIER) != 0) {

		//
		// allocate matrices
		// (note that we use 0-index arrays)
		//

		// dimension is the same number of elements as in bottom level,
		// due to the nature of the hierarchical preconditioner (leaf
		// elements on the diagonal)

		// q and h matrices
		if(m_pclsGmres1_q == NULL) {
			// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
			SAFENEW_ARRAY_NOMEM_RET(CLin_Vector, m_pclsGmres1_q, SOLVE_GMRES_ITER_MAX+1)
		}
		if(m_pclsGmres1_q == NULL) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += (SOLVE_GMRES_ITER_MAX+1) * sizeof(CLin_Vector);
		// init first and second element
		ret = m_pclsGmres1_q[0].newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
		ret = m_pclsGmres1_q[1].newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);

		if(m_pclsGmres1_h == NULL) {
			// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
			SAFENEW_ARRAY_NOMEM_RET(CLin_Vector, m_pclsGmres1_h, SOLVE_GMRES_ITER_MAX)
		}
		if(m_pclsGmres1_h == NULL) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += SOLVE_GMRES_ITER_MAX * sizeof(CLin_Vector);
		// init first element
		ret = m_pclsGmres1_h[0].newsize(2);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += 2 * sizeof(double);

		//
		// declare vectors (note that we use
		// 0-index arrays)
		//
		// g, v, w, z, y vectors
		ret = m_clsGmres1_g.newsize(numElems_0+1);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += (numElems_0 + 1) * sizeof(double);
		ret = m_clsGmres1_v.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
		ret = m_clsGmres1_w.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
		ret = m_clsGmres1_z.newsize(numElems_0+1);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += (numElems_0 + 1) * sizeof(double);
		ret = m_clsGmres1_y.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
		// c,s vectors
		ret = m_clsGmres1_c.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
		ret = m_clsGmres1_s.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
		// r vector
		ret = m_clsGmres1_r.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);

		// Pq vector (preconditioned vector)
		ret = m_clsGmres1_Pq.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);
		// x0 vector (initial solution vector)
		ret = m_clsGmres1_x0.newsize(numElems_0);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += numElems_0 * sizeof(double);

		// record up to which iteration the arrays have been pre-allocated
		m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_1_LEVEL] = 0;
	}


	if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0) {

		// super preconditioner
		//
		// allocate only if not already allocated (in previous call)
		if(m_pclsSupPrecondMtx == NULL) {
			// SAFENEW_MATRIX_RET(TYPE, VAR, LEN, MEM)
			SAFENEW_MATRIX_RET(double, m_pclsSupPrecondMtx, SOLVE_MAX_SUPER_PRECOND_NUM, g_clsMemUsage.m_ulPrecondMem)
		}
		if(m_pclsSupPotMtx == NULL) {
			// SAFENEW_MATRIX_RET(TYPE, VAR, LEN, MEM)
			SAFENEW_MATRIX_RET(double, m_pclsSupPotMtx, SOLVE_MAX_SUPER_PRECOND_NUM, g_clsMemUsage.m_ulPrecondMem)
		}

		// allocate arrays for super preconditioner
		//
		// super preconditioner elements
		if(m_clsSupPrecondElements == NULL) {
			// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
			SAFENEW_ARRAY_RET(CSuperPrecondElement, m_clsSupPrecondElements, SOLVE_MAX_SUPER_PRECOND_NUM, g_clsMemUsage.m_ulPrecondMem)
		}

		// leaf pointers to super preconditioner high level leaves;
		// in case it was already declared, destroy and re-declare (can change length)
		if(m_puiSupPrecondIndex != NULL) {
			delete m_puiSupPrecondIndex;
		}
		// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
		SAFENEW_ARRAY_RET(unsigned int, m_puiSupPrecondIndex, numPanels_0, g_clsMemUsage.m_ulPrecondMem)

		// array of leaf panel areae, accessible via index, to correctly extend / reduce
		// the results of the super precond multiplication to the leaves (correct weighting,
		// in case of non-uniform panel discretizations)
		if(m_pfSupPrecondAreae != NULL) {
			delete m_pfSupPrecondAreae;
		}
		// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
		SAFENEW_ARRAY_RET(float, m_pfSupPrecondAreae, numPanels_0, g_clsMemUsage.m_ulPrecondMem)
	}

	// if block preconditioner
	if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0) {

		// allocate arrays for block preconditioner
		//

		// block preconditioner elements
		if(m_clsBlockPrecondElements == NULL) {
			// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
			SAFENEW_ARRAY_RET(CAutoPanel*, m_clsBlockPrecondElements, SOLVE_MAX_BLOCK_PRECOND_NUM, g_clsMemUsage.m_ulPrecondMem)
		}
		// block preconditioner
		if(m_pdBlockPrecond != NULL) {
			delete []m_pdBlockPrecond;
		}
		// SAFENEW_MATRIXNM_RET(TYPE, VAR, LEN, MEM)
		SAFENEW_MATRIXNM_RET(double, m_pdBlockPrecond, numElems_0, SOLVE_MAX_BLOCK_PRECOND_NUM, g_clsMemUsage.m_ulPrecondMem)

		// number of elements in each block (in case there are less than 'm_clsGlobalVars.m_uiBlockPreSize')
		if(m_pucBlockPrecondDim != NULL) {
			delete m_pucBlockPrecondDim;
		}
		// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
		SAFENEW_ARRAY_RET(unsigned char, m_pucBlockPrecondDim, numElems_0, g_clsMemUsage.m_ulPrecondMem);

		// borrow also the super potential matrix
		if(m_pclsSupPotMtx == NULL) {
			// SAFENEW_MATRIX_RET(TYPE, VAR, LEN, MEM)
			SAFENEW_MATRIX_RET(double, m_pclsSupPotMtx, SOLVE_MAX_SUPER_PRECOND_NUM, g_clsMemUsage.m_ulPrecondMem)
		}
		// and allocate temporary block to store the inverted 'm_pclsSupPotMtx' when calculating the precond
		if(m_pclsTmpCapMtx == NULL) {
			// SAFENEW_MATRIX_RET(TYPE, VAR, LEN, MEM)
			SAFENEW_MATRIX_RET(double, m_pclsTmpCapMtx, SOLVE_MAX_SUPER_PRECOND_NUM, g_clsMemUsage.m_ulPrecondMem)
		}
	}

	return FC_NORMAL_END;
}

// delete structures
void CSolveCap::DeallocateMemory(int command, CAutoRefGlobalVars globalVars)
{
	m_clsMulthier.Clean(command, globalVars);

	if(m_pCondCharges != NULL) {
		delete [] m_pCondCharges;
		m_pCondCharges = NULL;
	}

	if(m_pclsTmpCapMtx != NULL) {
		delete [] m_pclsTmpCapMtx;
		m_pclsTmpCapMtx = NULL;
	}

	if(m_pclsGmres_q != NULL) {
		delete []m_pclsGmres_q;
		m_pclsGmres_q = NULL;
	}
	if(m_pclsGmres_h != NULL) {
		delete []m_pclsGmres_h;
		m_pclsGmres_h = NULL;
	}
	if(m_pclsGmres_zf != NULL) {
		delete []m_pclsGmres_zf;
		m_pclsGmres_zf = NULL;
	}
	if(m_pclsGmres1_q != NULL) {
		delete []m_pclsGmres1_q;
		m_pclsGmres1_q = NULL;
	}
	if(m_pclsGmres1_h != NULL) {
		delete []m_pclsGmres1_h;
		m_pclsGmres1_h = NULL;
	}

	DeallocatePrecond();

	m_clsGmres_g.destroy();
	m_clsGmres_v.destroy();
	m_clsGmres_w.destroy();
	m_clsGmres_z.destroy();
	m_clsGmres_y.destroy();
	m_clsGmres_c.destroy();
	m_clsGmres_s.destroy();
	m_clsGmres_r.destroy();
	m_clsGmres_Pq.destroy();
	m_clsGmres_x0.destroy();

	m_clsGmres1_g.destroy();
	m_clsGmres1_v.destroy();
	m_clsGmres1_w.destroy();
	m_clsGmres1_z.destroy();
	m_clsGmres1_y.destroy();
	m_clsGmres1_c.destroy();
	m_clsGmres1_s.destroy();
	m_clsGmres1_r.destroy();
	m_clsGmres1_Pq.destroy();
	m_clsGmres1_x0.destroy();

	g_clsMemUsage.m_ulGmresMem = 0;
	g_clsMemUsage.m_ulPrecondMem = 0;
}

void CSolveCap::DeallocatePrecond()
{
	// used in super and hierarchical preconditioner

	if(m_pclsSupPrecondMtx != NULL) {
		delete [] m_pclsSupPrecondMtx;
		m_pclsSupPrecondMtx = NULL;
	}
	if(m_pclsSupPotMtx != NULL) {
		delete [] m_pclsSupPotMtx;
		m_pclsSupPotMtx = NULL;
	}
	if(m_clsSupPrecondElements != NULL) {
		delete m_clsSupPrecondElements;
		m_clsSupPrecondElements = NULL;
	}
	if(m_puiSupPrecondIndex != NULL) {
		delete m_puiSupPrecondIndex;
		m_puiSupPrecondIndex = NULL;
	}
	if(m_pfSupPrecondAreae != NULL) {
		delete m_pfSupPrecondAreae;
		m_pfSupPrecondAreae = NULL;
	}


	// used in block preconditioner

	if(m_clsBlockPrecondElements != NULL) {
		delete m_clsBlockPrecondElements;
		m_clsBlockPrecondElements = NULL;
	}
	if(m_pdBlockPrecond != NULL) {
		delete [] m_pdBlockPrecond;
		m_pdBlockPrecond = NULL;
	}
	if(m_pucBlockPrecondDim != NULL) {
		delete m_pucBlockPrecondDim;
		m_pucBlockPrecondDim = NULL;
	}

}

// recursively count panels and build super preconditioner
void CSolveCap::RecurseComputePrecond(CAutoElement* element)
{
	InteractionC3DList::iterator iti;
	bool calculateBlock;

	m_iLevel++;

	// only if using Block precond
	calculateBlock = false;
	if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0 ) {

		// if we have reached the right depth, or if from the start
		// there is not enough depth, reset element pointer array
		if(element->m_lNumOfChildren <= m_clsGlobalVars.m_uiBlockPreSize &&
		        m_bIsComputingBlock == false) {

			m_uiBlockPreNum = 0;
			m_bIsComputingBlock = true;
			calculateBlock = true;
		}
	}

	if(element->IsLeaf() == true) {

		// compute Jacobi preconditioner
//		if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_JACOBI) != 0 ) {
//
// old, when Jacobi precond was not calculated on the spot, and there was no dedicated array for self-potential
//			// self potential is always the first one, see RefineSelf()
//			m_clsJacPrecond[m_ulPanelNum] = 1 / m_clsMulthier.GetAutoPotCoeff(element);
//
// older
//			for(iti = element->m_stlInteractions[m_clsMulthier.GetInteractionLevel()].begin(); iti != element->m_stlInteractions[m_clsMulthier.GetInteractionLevel()].end(); iti++) {
//				if( (*iti).m_clsPanel == element ) {
//					m_clsJacPrecond[m_ulPanelNum] = 1 / (*iti).m_dPotCoeff;
//					break;
//				}
//			}
//		}

		// only if using Block precond, store current leaf element in the block precond elements array
		if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0 ) {
			m_clsBlockPrecondElements[m_uiBlockPreNum] = (CAutoPanel*)element;
			m_uiBlockPreNum++;
		}

		// only if using Super precond
		if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0) {

			// if in this branch there is no depth enough, use leaf for preconditioner
			if(m_iLevel <= (int)m_uiSupPreLevel) {
				// hi level preconditioner elements array
				m_clsSupPrecondElements[m_uiSupPreNum].m_pclsHiLevLeaf = element;
//				m_clsSupPrecondElements[m_uiSupPreNum].m_lLoLevLeavesNum = 1;
				m_clsSupPrecondElements[m_uiSupPreNum].m_pCond = m_pCurrCond;
				// leaf pointer to hi level preconditioner element
				m_puiSupPrecondIndex[m_ulPanelNum] = m_uiSupPreNum;
				IncrementSupPreNum();
			}
			// otherwise, just record pointer from leaf to hi level preconditioner element
			// and increment count of leaves under that hi level preconditioner element
			else {
				m_puiSupPrecondIndex[m_ulPanelNum] = m_ulTmpSupPreIndex;
				m_ulSupPreNumOfLeaves++;
			}

			// and copy the area
			if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
				m_pfSupPrecondAreae[m_ulPanelNum] = ((CAutoPanel*)element)->GetArea();
			}
			else {
				m_pfSupPrecondAreae[m_ulPanelNum] = ((CAutoSegment*)element)->GetLength();
			}
		}

		m_ulPanelNum++;
	}
	// else, element is not a leaf
	else {
		// only if using Super precond
		if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) != 0 ) {

			// if correct depth for super preconditioner
			if(m_iLevel == (int)m_uiSupPreLevel) {
				// initialize count of number of subtended leaves
				m_ulSupPreNumOfLeaves = 0;
				// and remember what is the hi level preconditioner element
				m_ulTmpSupPreIndex = m_uiSupPreNum;
				// hi level preconditioner elements array
				m_clsSupPrecondElements[m_uiSupPreNum].m_pclsHiLevLeaf = element;
				m_clsSupPrecondElements[m_uiSupPreNum].m_pCond = m_pCurrCond;
			}
		}

		// if not a leaf element, go into left and right sub-trees
		RecurseComputePrecond(element->m_pLeft);
		RecurseComputePrecond(element->m_pRight);

		// only if using Block precond, store current leaf element in the block precond elements array
		if( (m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_BLOCK) != 0 && calculateBlock == true) {
			m_bIsComputingBlock = false;
			ComputeBlockPrecond();
		}

		// only if using Super precond or Hier precond and Super precond can be used to accelerate iterations
		if(m_clsGlobalVars.m_ucPrecondType & AUTOREFINE_PRECOND_SUPER) {

			// if correct depth for super preconditioner
			if(m_iLevel == (int)m_uiSupPreLevel) {
				// store number of leaves
//				m_clsSupPrecondElements[m_uiSupPreNum].m_lLoLevLeavesNum = m_ulSupPreNumOfLeaves;
				IncrementSupPreNum();
			}
		}
	}

	m_iLevel--;
}

// safe function to increment m_uiSupPreNum
void CSolveCap::IncrementSupPreNum()
{
	if(m_uiSupPreNum < SOLVE_MAX_SUPER_PRECOND_NUM) {
		m_uiSupPreNum++;
	}
	else {
		ErrMsg("Internal error: number of two-levels preconditioner elements greater than %i\nReverting to Jacobi preconditioner\n", SOLVE_MAX_SUPER_PRECOND_NUM);
		m_clsGlobalVars.m_ucPrecondType = AUTOREFINE_PRECOND_JACOBI;
	}
}

void CSolveCap::ComputeSuperPrecond()
{
	unsigned int i, j;
	int isPotValid;
	double potestim1, potestim2;

	// compute super potential matrix
	for(i = 0; i < m_uiSupPreNum; i++) {
		for(j = i; j < m_uiSupPreNum; j++) {
			// if self potential
			if( i == j ) {
				m_clsMulthier.SetCurrentConductor(m_clsSupPrecondElements[i].m_pCond);
				if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
					m_clsMulthier.SelfPotential((CAutoPanel*)(m_clsSupPrecondElements[i].m_pclsHiLevLeaf), &potestim1, &potestim2);
				}
				else {
					m_clsMulthier.SelfPotential((CAutoSegment*)(m_clsSupPrecondElements[i].m_pclsHiLevLeaf), &potestim1, &potestim2);
				}
				// store only the real part; the super preconditioner in case of complex capacitance matrix
				// will be used to precondition only the real part of the potential matrix, since this is
				// the one compressed with the hierarchical method
				m_pclsSupPotMtx[i][i] = potestim1;
			}
			else {
				if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
					isPotValid = m_clsMulthier.PotEstimateOpt((CAutoPanel*)(m_clsSupPrecondElements[i].m_pclsHiLevLeaf), (CAutoPanel*)(m_clsSupPrecondElements[j].m_pclsHiLevLeaf), potestim1, potestim2, AUTOREFINE_PRECOND_SUPER);
				}
				else {
					isPotValid = m_clsMulthier.PotEstimateOpt((CAutoSegment*)(m_clsSupPrecondElements[i].m_pclsHiLevLeaf), (CAutoSegment*)(m_clsSupPrecondElements[j].m_pclsHiLevLeaf), potestim1, potestim2, AUTOREFINE_PRECOND_SUPER);
				}

				// if real error while calculating the precond
				if(isPotValid != AUTOREFINE_NO_ERROR) {
					ErrMsg("Error: invalid potential calculation during two-levels preconditioner calculation\n");
					ErrMsg("       Removing the pre-conditioner and continuing\n");
					m_clsGlobalVars.m_ucPrecondType &= (~AUTOREFINE_PRECOND_SUPER);
					DeallocatePrecond();
					return;
				}

				// if problems calculating the precond, and in automatic mode, prevent usage of the preconditioner
				if(	m_clsMulthier.m_clsGlobalVars.m_bWarnGivenPre == true && m_clsGlobalVars.m_bAutoPrecond == true) {
					ErrMsg("Warning: problems found when forming the two-levels preconditioner in automatic preconditioner mode\n");
					ErrMsg("         Moving to Jacobi pre-conditioner, to avoid possible worst convergence of GMRES, and continuing\n");
					m_clsGlobalVars.m_ucPrecondType = AUTOREFINE_PRECOND_JACOBI;
					DeallocatePrecond();
					return;
				}

				m_pclsSupPotMtx[i][j] = potestim1;
				m_pclsSupPotMtx[j][i] = potestim2;
			}
		}
	}

	// then invert it
	//

	InvertMatrix(m_pclsSupPotMtx, m_pclsSupPrecondMtx, m_uiSupPreNum);
}

void CSolveCap::ComputeBlockPrecond()
{
	unsigned int i, j;
	int isPotValid;
	double potestim1, potestim2;

	ASSERT(m_uiBlockPreNum <= m_clsGlobalVars.m_uiBlockPreSize);

	for(i=0; i < m_uiBlockPreNum; i++) {
		for(j = i; j < m_uiBlockPreNum; j++) {
			// if self potential
			if( i == j ) {
				m_clsMulthier.SelfPotential(m_clsBlockPrecondElements[i], &potestim1, &potestim2);
				// store only the real part; the super preconditioner in case of complex capacitance matrix
				// will be used to precondition only the real part of the potential matrix, since this is
				// the one compressed with the hierarchical method
				m_pclsSupPotMtx[i][i] = potestim1;
			}
			else {
				isPotValid = m_clsMulthier.PotEstimateOpt(m_clsBlockPrecondElements[i], m_clsBlockPrecondElements[j], potestim1, potestim2, AUTOREFINE_PRECOND_BLOCK);

				if(isPotValid != AUTOREFINE_NO_ERROR) {
					ErrMsg("Error: invalid potential calculation during two-levels preconditioner calculation\n");
					ErrMsg("       Efficiency of the preconditioner may be impacted\n");
				}

				m_pclsSupPotMtx[i][j] = potestim1;
				m_pclsSupPotMtx[j][i] = potestim2;
			}
		}
	}

	InvertMatrix(m_pclsSupPotMtx, m_pclsTmpCapMtx, m_uiBlockPreNum);

	// store block dimension for preconditioner usage
	m_pucBlockPrecondDim[m_ulBlockPreBaseNum] = (unsigned char) m_uiBlockPreNum;
	// copy inverted matrix into preconditioner block
	for(i=0; i < m_uiBlockPreNum; i++) {
		for(j=0; j < m_uiBlockPreNum; j++) {
			m_pdBlockPrecond[m_ulBlockPreBaseNum + i][j] = m_pclsTmpCapMtx[i][j];
		}
	}

	// and position index to next block
	m_ulBlockPreBaseNum += m_uiBlockPreNum;
}

// REMARK: can invert matrices with up to 64k entries (max size is 'unsigned int')
void CSolveCap::InvertMatrix(double (*matrix)[SOLVE_MAX_SUPER_PRECOND_NUM], double (*invMatrix)[SOLVE_MAX_SUPER_PRECOND_NUM], unsigned long size)
{
	unsigned int j, k, jp, ii, jj, ip;
	int i, *indx;
	double t, sum, recp, *b;
	bool flag;

	indx = new int[size];
	b = new double[size];

	// algorithm taken from LinAlg
	// LU factorization with pivoting; if mtx is not singular,
	// LU factorization is always possible

	// right-looking LU factorization algorithm (unblocked)
	//
	//   Factors matrix A into lower and upper  triangular matrices
	//   (L and U respectively)
	//
	//
	// A        (input/output) Matrix(1:n, 1:n)  In input, matrix to be
	//                  factored.  On output, overwritten with lower and
	//                  upper triangular factors.
	//
	// indx     (output) Vector(1:n)    Pivot vector. Describes how
	//                  the rows of A were reordered to increase
	//                  numerical stability.
	//
	for (j=0; j < size; j++) {
		// find pivot in column j and test for singularity.

		jp = j;
		t = fabs(matrix[j][j]);
		for (i=j+1; i<(long)size; i++) {
			if ( fabs(matrix[i][j]) > t) {
				jp = i;
				t = fabs(matrix[i][j]);
			}
		}

		indx[j] = jp;

		// jp now has the index of maximum element
		// of column j, below the diagonal

		// factorization failed because of zero pivot
		// TBC warning: maybe should add a real warning?
		ASSERT( matrix[jp][j] != 0 );

		// if pivot not already on the diagonal
		if (jp != j) {
			// swap rows j and jp
			for (k=0; k<size; k++) {
				t = matrix[j][k];
				matrix[j][k] = matrix[jp][k];
				matrix[jp][k] = t;
			}
		}

		// divide elements j+1:panelNum of jth column by pivot element
		// (for the first panelNum-1 cols, last col's below-the-diagonal element
		// is already the pivot)
		if (j < size - 1) {
			// note that potentialMtx(j,j) was previously potentialMtx(jp,p), which was
			// guaranteed not to be zero
			recp =  1.0 / matrix[j][j];

			for (k=j+1; k<size; k++)
				matrix[k][j] *= recp;
		}


		if (j < size - 1) {
			// rank-1 update to trailing submatrix:   E = E - x*y;
			//
			// E is the region potentialMtx(j+1:M, j+1:N)
			// x is the column vector potentialMtx(j+1:M,j)
			// y is row vector potentialMtx(j,j+1:N)

			for (ii=j+1; ii<size; ii++)
				for (jj=j+1; jj<size; jj++)
					matrix[ii][jj] -= matrix[ii][j]*matrix[j][jj];
		}
	}

	// now we have the LU factors of 'potentialMtx', stored in-place in 'potentialMtx'.
	// So we compute inv(potentialMtx) from the LU factors.
	//
	// Remark: to get results in the correct order, must remember pivot vector
	// 'indx', which stores rows permutation

	for(k=0; k<size; k++) {

		// compute b vector
		for (i=0; i<(long)size; i++) {
			b[i] = 0.0;
		}
		b[k] = 1.0;

		// forward substitution
		flag = false;
		sum = 0.0;

		for (i=0; i<(long)size; i++) {
			ip = indx[i];
			sum = b[ip];
			b[ip] = b[i];

			if (flag) {
				for (j=ii; j<=(unsigned long)(i-1); j++) {
					sum -= matrix[i][j]*b[j];
				}
			}
			else if(sum != 0.0) {
				ii = i;
				flag = true;
			}

			b[i] = sum;
		}

		// backward substitution
		for (i=size-1; i>=0; i--) {
			sum=b[i];
			for (j=i+1; j<size; j++) {
				sum -= matrix[i][j]*b[j];
			}
			b[i] = sum/matrix[i][i];
		}

		// copy the column 'x' resulting from solving A*x = b ('x' is stored in-place in 'b' vector)
		// into the precond matrix
		for(i=0; i < (long)size; i++) {
			_ASSERT(!isnan(b[i]));
			_ASSERT(isfinite(b[i]));

			invMatrix[i][k] = b[i];
		}
	}

	delete [] indx;
	delete [] b;
}

//   Generalized Minimum Residual Method with preconditioner
//
//   x = gmresPrecond(A,b,P) attempt to solve the preconditioned
//   system A*P*y = b & x = P*y
//
//   Enrico Di Lorenzo, 2004/07/06
//
//   'fast' refers to the direct vector multiplication loops used
//   instead of CLin_Vector class operator*, which requires
//   the creation of a temp object for every assignment operation
//
//   the 's' on the other hand refers to the scilab gmres correction,
//   to make convergence faster, that basically uses as residual
//   a weighted residual: normr = || b - Ax || / || b ||, that is
//   normr = fabs(g[i+1]) / normb;
//   and not
//   normr = fabs(g[i+1]);
//   see also the notes inline about this choice
//
//   the 'all' refers to the increase of speed not allocating memory
//   locally every time
//
//   This version has been modified to work with a complex potential matrix;
//   this is not a general-purpose gmres dealing with complex values, but an optimized
//   version that considers that the only imaginary parts are on the diagonal, and only
//   for the rows corresponding to the dielectric panels. In this case, b and x have
//   a length 2x the number of panels n, where the first n elements are the real parts,
//   while the remaining n elements are the imaginary parts
int CSolveCap::gmresPrecondSFastAll(CLin_Vector *b, CLin_Vector *x, double gmresTol)
{
	double normr, normb, tmp1, tmp2, length;
	long i, j, iteration, size, k, ii;
	bool ret;
	int retInt;
	CLin_Range bre, xre;

	// get system size
	size = b->size();
	// check consistency of size
	ASSERT(size == (long)(*x).size());

	// if there is a preconditioner
	if(m_clsGlobalVars.m_ucPrecondType != AUTOREFINE_PRECOND_NONE) {

		// initial vector: having the approx inverted matrix, it is cheap to calculate the initial vector.
		//  Simply put, to solve the preconditioned system
		//  A * P * y = b and x = P * y, the best choice for y0 initial vector is b,
		//  since a good P should be as similar as possible to inv(A), so A * P ~= I
		//  So x0 = P * y0 = P * b
		ComputePrecondVectFast(&m_clsGmres_x0, b, m_clsGlobalVars.m_ucPrecondType);

		// r = b - A * x0;
		// matrix - vector multiplication
		retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, &m_clsGmres_x0);
		if(retInt != FC_NORMAL_END) {
			return retInt;
		}
		// and get r
		//m_clsGmres_r = *b - m_clsGmres_v;
		for(k=0; k<size; k++) {
			m_clsGmres_r[k] = (*b)[k] - m_clsGmres_v[k];
		}
		normr = mod(m_clsGmres_r);
		normb = mod(*b);
	}
	else {
		// no preconditioner; assumes that initial vector x0 is an all zeros vector
		//
		m_clsGmres_r = *b;
		normr = mod(m_clsGmres_r);
		normb = normr;
	}

	// allocate and init first column of Q matrix
	//m_pclsGmres_q[0] = m_clsGmres_r / normr;
	for(k=0; k<size; k++) {
		m_pclsGmres_q[0][k] = m_clsGmres_r[k] / normr;
	}

	// init first element of first row of the
	// Q matrix of the QR factors of H; this should be 1,
	// using normr will produce the first row
	// multiplied by normr, which is used to compute the solution
	// Remark: this is not the Q matrix whose columns are q(i,j)
	m_clsGmres_g[0] = normr;

	//
	// if the norm of the residual is small enough,
	// initial vector is a good enough solution
	//
	if(normr / normb < gmresTol) {
		if(m_clsGlobalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE) {
			*x = CLin_Vector(size, 0.0);
		}
		else {
			*x = m_clsGmres_x0;
		}
		return FC_NORMAL_END;
	}

	LogMsg("GMRES Iteration: ");

	// start iteration
	for(i = 0; i < SOLVE_GMRES_ITER_MAX && i < (long)size; i++) {

		if(g_bFCContinue == false) {
			return FC_USER_BREAK;
		}

		//
		// execute the i-th Arnoldi step, thus computing
		// the i-th column of the H (upper Hessemberg) matrix
		//

		if(m_clsGlobalVars.m_bDumpResidual == true) {
			LogMsg("%.3g ", normr);
		}
		LogMsg("%d ", i);

		// allocate only if not already pre-allocated
		if((long)m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_0_LEVEL] < i) {

			ASSERT((long)m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_0_LEVEL] == i-1);

			// allocate new column of h
			// (remark: since vectors are 0-based, must
			// allocate i+2 entries and not only i+1;
			// i.e. variable 'i' starts at 0!)
			ret = m_pclsGmres_h[i].newsize(i+2);
			if(ret == false) {
				return FC_OUT_OF_MEMORY;
			}
			g_clsMemUsage.m_ulGmresMem += (i+2) * sizeof(double);
			// allocate new column of q
			ret = m_pclsGmres_q[i+1].newsize(size);
			if(ret == false) {
				return FC_OUT_OF_MEMORY;
			}
			g_clsMemUsage.m_ulGmresMem += size * sizeof(double);

			m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_0_LEVEL]++;
		}

		// compute new vector
		//
		if(m_clsGlobalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE) {
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, &m_pclsGmres_q[i]);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}
		else {
			// compute preconditioned vector
			ComputePrecondVectFast(&m_clsGmres_Pq, &m_pclsGmres_q[i], m_clsGlobalVars.m_ucPrecondType);
			// and use it in matrix - vector multiplication
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, &m_clsGmres_Pq);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}

		//m_clsGmres_w = m_clsGmres_v;
		for(k=0; k<size; k++) {
			m_clsGmres_w[k] = m_clsGmres_v[k];
		}

		// orthogonalize it
		for(j=0; j<=i; j++) {
			m_pclsGmres_h[i][j] = dot_prod(m_pclsGmres_q[j], m_clsGmres_w);
			for(ii=0; ii<size; ii++) {
				m_clsGmres_v[ii] = m_clsGmres_v[ii] - m_pclsGmres_h[i][j] * m_pclsGmres_q[j][ii];
			}
		}

		m_pclsGmres_h[i][i+1] = mod(m_clsGmres_v);

		for(ii=0; ii<size; ii++) {
			m_pclsGmres_q[i+1][ii] = m_clsGmres_v[ii] / (m_pclsGmres_h[i][i+1]);
		}

		//
		// rotate new vector (Givens rotations) to compute
		// the (R;0) matrix of the QR factors of H(j+1);
		// this is done updating (Q';tmp')*H(j) = (R(j);0)
		// to produce (Q';tmp')*H(j+1)=(R(j+1);0)
		// Remark: the first time the loop is executed,
		// this rotation does not happen
		//

		// apply all old rotations to new i-th column
		// of the H matrix to get Q*H from H, but for the
		// last two elements
		for(j=0; j<i; j++) {
			tmp1 = m_pclsGmres_h[i][j];
			tmp2 = m_pclsGmres_h[i][j+1];
			m_pclsGmres_h[i][j] = m_clsGmres_c[j] * tmp1 - m_clsGmres_s[j]* tmp2;
			m_pclsGmres_h[i][j+1] = m_clsGmres_c[j] * tmp2 + m_clsGmres_s[j] * tmp1;
		}

		// compute the new Givens rotation to annihilate h[i][i+1]
		tmp1 = m_pclsGmres_h[i][i];
		tmp2 = m_pclsGmres_h[i][i+1];
		length = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		m_clsGmres_c[i] = tmp1/length;
		m_clsGmres_s[i] = -tmp2/length;

		// apply new Givens rotation to incomplete
		// (Q';tmp')*H(j+1) matrix i-th column
		// to get the complete (R(j+1);0)
		m_pclsGmres_h[i][i] = m_clsGmres_c[i] * tmp1 - m_clsGmres_s[i] * tmp2;
		// this should be zero anyway
		//h(i+1,i) = c(i) * tmp2 + s(i) * tmp1;
		m_pclsGmres_h[i][i+1] = 0;

		//
		// Compute recursively the first column of the (Q';tmp')
		// matrix of the QR factors of H(j+1); this is used to
		// get implicitly the norm of the residual and to compute
		// the solution once the norm is small enough
		//

		tmp1 = m_clsGmres_g[i];
		m_clsGmres_g[i] = m_clsGmres_c[i] * tmp1;
		m_clsGmres_g[i+1] = m_clsGmres_s[i] * tmp1;

		// this is the norm of the residual
		// remark: || b - A*Qj*z || = normb * fabs(g[i+1])
		// (see "Matrix analysis and applied linear algebra", Meyers)
		// so fabs(g[i+1]) is already the norm of the residual
		// weighted by the norm of b (i.e the % difference)
		// However this is not true here, also verified
		// explicitly calculating the residual norm, once the solution
		// is found (see below), so the weighted formula is correct
		normr = fabs(m_clsGmres_g[i+1]) / normb;

		//
		// if the norm of the residual is small enough,
		// exit loop (and return solution)
		//
		if(normr < gmresTol) {
			break;
		}
	}

	if(i >= SOLVE_GMRES_ITER_MAX || i >= size) {
		ErrMsg("\nError: not converging after %d iterations, norm of the residual is %.3f, while targeting %.3f\n", i, normr, gmresTol);
		// either case, since we ended the 'for' loop, 'i' has been incremented of one more;
		// to avoid breaking the end of the arrays here below, we must decrease 'i'
		i--;
	}

	if(m_clsGlobalVars.m_bDumpResidual == true) {
		LogMsg("%.3f ", normr);
	}
	LogMsg("\n");

	//
	// compute the solution, solving H(i)*z = normr * e1 for z using
	// the first row of the Q matrix and the R matrix of the QR
	// factorization of H; then the solution is x = Q(i)*z
	//

	iteration = i;

	// H(i)*z = normb*e1 gives Q*R*z = normb*e1 -> R*z = Q'*normb*e1
	// Since g(i) is already the first column of Q' multiplied by normb,
	// it is the term Q'*normb*e1, so init the solution
	//m_clsGmres_z = m_clsGmres_g;
	for(k=0; k<size; k++) {
		m_clsGmres_z[k] = m_clsGmres_g[k];
	}

	// then solve R*z = Q'*normr*e1 by back substitution,
	// knowing that R is upper triangular
	for(i = iteration; i>=0; i--) {
		m_clsGmres_z[i]= m_clsGmres_z[i] / m_pclsGmres_h[i][i];
		for(j = i-1; j>=0; j--) {
			m_clsGmres_z[j] = m_clsGmres_z[j] - m_pclsGmres_h[i][j] * m_clsGmres_z[i];
		}
	}


	// multiply z by Q(i) to get y, the solution
	// of the (possibly preconditioned) system
	for(i=0; i<size; i++) {
		for(j=0, tmp1 = 0; j<=iteration; j++) {
			tmp1 += m_pclsGmres_q[j][i] * m_clsGmres_z[j];
		}
		m_clsGmres_y[i] = tmp1;
	}

	if(m_clsGlobalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE) {
		*x = m_clsGmres_y;
	}
	else {
		// undo the preconditioner to get x, the final solution,
		// also considering initial vector in the solution
		//x = x0 + P * y;
		ComputePrecondVectFast(x, &m_clsGmres_y, m_clsGlobalVars.m_ucPrecondType);
		//*x = m_clsGmres_x0 + *x;
		for(k=0; k<size; k++) {
			(*x)[k] = m_clsGmres_x0[k] + (*x)[k];
		}
	}

	// debug : explicit calculation of residual vector
	// (is equal to normr*normb in this routine; this proves
	// that normr here is weighted on normb)
	//double residual2;
	//CLin_Vector b_approx(size);
	//m_clsMulthier.MultiplyMatByVec_fast(x, &b_approx);
	//residual2 = mod(*b - b_approx);

	return FC_NORMAL_END;
}

//   Generalized Minimum Residual Method with preconditioner
//
//   x = gmresPrecond(A,b,P) attempt to solve the preconditioned
//   system A*P*y = b & x = P*y
//
//   Enrico Di Lorenzo, 2004/07/06
//
//   'fast' refers to the direct vector multiplication loops used
//   instead of CLin_Vector class operator*, which requires
//   the creation of a temp object for every assignment operation
//
//   the 's' on the other hand refers to the scilab gmres correction,
//   to make convergence faster, that basically uses as residual
//   a weighted residual: normr = || b - Ax || / || b ||, that is
//   normr = fabs(g[i+1]) / normb;
//   and not
//   normr = fabs(g[i+1]);
//   see also the notes inline about this choice
//
//   the 'all' refers to the increase of speed not allocating memory
//   locally every time
//
//   the 'upper' refers to the fact that this function is used in
//   upper iteration of hierarchical precond, so it is using
//   other global vars than 'gmresPrecondSFastAll', since it will
//   be called for each iteration step
int CSolveCap::gmresPrecondSFastAllUpper(CLin_Vector *b, CLin_Vector *x, double gmresTol, unsigned char precondType)
{
	double normr, normb, tmp1, tmp2, length;
	long i, ii, size, j, iteration;
	bool ret;
	int retInt;

	// get system size
	size = b->size();
	// check consistency of size
//	ASSERT(size == *x.size());

	// if there is a preconditioner
	if(precondType != AUTOREFINE_PRECOND_NONE) {

		// initial vector: having the approx inverted matrix, it is cheap to calculate the initial vector.
		//  Simply put, to solve the preconditioned system
		//  A * P * y = b and x = P * y, the best choice for y0 initial vector is b,
		//  since a good P should be as similar as possible to inv(A), so A * P ~= I
		//  So x0 = P * y0 = P * b
		ComputePrecondVectFast(&m_clsGmres1_x0, b, precondType);

		// r = b - A * x0;
		// matrix - vector multiplication
		retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres1_v, &m_clsGmres1_x0);
		if(retInt != FC_NORMAL_END) {
			return retInt;
		}

		// and get r
		m_clsGmres1_r = *b - m_clsGmres1_v;
		normr = mod(m_clsGmres1_r);
		normb = mod(*b);
	}
	else {
		// no preconditioner; assumes that initial vector x0 is an all zeros vector
		//
		m_clsGmres1_r = *b;
		normr = mod(m_clsGmres1_r);
		normb = normr;
	}

	// allocate and init first column of Q matrix
	m_pclsGmres1_q[0] = CLin_Vector(m_clsGmres1_r / normr);

	// init first element of first row of the
	// Q matrix of the QR factors of H; this should be 1,
	// using normr will produce the first row
	// multiplied by normr, which is used to compute the solution
	// Remark: this is not the Q matrix whose columns are q(i,j)
	m_clsGmres1_g[0] = normr;

	LogMsg("\nUpper GMRES Iteration: ");

	//
	// if the norm of the residual is small enough,
	// initial vector is a good enough solution
	//
	if(normr / normb < gmresTol) {

		// dump residual and first iteration
		if(m_clsGlobalVars.m_bDumpResidual == true) {
			LogMsg("%.3g ", normr / normb);
		}
		LogMsg("%d \n", 0);

		if(precondType == AUTOREFINE_PRECOND_NONE) {
			*x = CLin_Vector(size, 0.0);
		}
		else {
			*x = m_clsGmres1_x0;
		}
		return 0;
	}

	// start iteration
	for(i = 0; i < SOLVE_GMRES_ITER_MAX && i < size; i++) {

		if(g_bFCContinue == false) {
			return FC_USER_BREAK;
		}

		//
		// execute the i-th Arnoldi step, thus computing
		// the i-th column of the H (upper Hessemberg) matrix
		//

		if(m_clsGlobalVars.m_bDumpResidual == true) {
			LogMsg("%.3g ", normr);
		}
		LogMsg("%d ", i);

		// allocate only if not already pre-allocated
		if((long)m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_1_LEVEL] < i) {

			ASSERT((long)m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_1_LEVEL] == i-1);

			// allocate new column of h
			// (remark: since vectors are 0-based, must
			// allocate i+2 entries and not only i+1;
			// i.e. variable 'i' starts at 0!)
			ret = m_pclsGmres1_h[i].newsize(i+2);
			if(ret == false) {
				return FC_OUT_OF_MEMORY;
			}
			g_clsMemUsage.m_ulGmresMem += (i + 2) * sizeof(double);
			// allocate new column of q
			ret = m_pclsGmres1_q[i+1].newsize(size);
			if(ret == false) {
				return FC_OUT_OF_MEMORY;
			}
			g_clsMemUsage.m_ulGmresMem += size * sizeof(double);

			m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_1_LEVEL]++;
		}

		// compute new vector
		//
		if(precondType == AUTOREFINE_PRECOND_NONE) {
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres1_v, &m_pclsGmres1_q[i]);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}
		else {
			// compute preconditioned vector
			ComputePrecondVectFast(&m_clsGmres1_Pq, &m_pclsGmres1_q[i], precondType);
			// and use it in matrix - vector multiplication
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres1_v, &m_clsGmres1_Pq);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}

		m_clsGmres1_w = m_clsGmres1_v;

		// orthogonalize it
		for(j=0; j<=i; j++) {
			m_pclsGmres1_h[i][j] = dot_prod(m_pclsGmres1_q[j], m_clsGmres1_w);
			for(ii=0; ii<(long)m_clsGmres1_v.size(); ii++) {
				m_clsGmres1_v[ii] = m_clsGmres1_v[ii] - m_pclsGmres1_h[i][j] * m_pclsGmres1_q[j][ii];
			}
		}

		m_pclsGmres1_h[i][i+1] = mod(m_clsGmres1_v);

		for(ii=0; ii<(long)m_clsGmres1_v.size(); ii++) {
			m_pclsGmres1_q[i+1][ii] = m_clsGmres1_v[ii] / (m_pclsGmres1_h[i][i+1]);
		}

		//
		// rotate new vector (Givens rotations) to compute
		// the (R;0) matrix of the QR factors of H(j+1);
		// this is done updating (Q';tmp')*H(j) = (R(j);0)
		// to produce (Q';tmp')*H(j+1)=(R(j+1);0)
		// Remark: the first time the loop is executed,
		// this rotation does not happen
		//

		// apply all old rotations to new i-th column
		// of the H matrix to get Q*H from H, but for the
		// last two elements
		for(j=0; j<i; j++) {
			tmp1 = m_pclsGmres1_h[i][j];
			tmp2 = m_pclsGmres1_h[i][j+1];
			m_pclsGmres1_h[i][j] = m_clsGmres1_c[j] * tmp1 - m_clsGmres1_s[j]* tmp2;
			m_pclsGmres1_h[i][j+1] = m_clsGmres1_c[j] * tmp2 + m_clsGmres1_s[j] * tmp1;
		}

		// compute the new Givens rotation to annihilate h[i][i+1]
		tmp1 = m_pclsGmres1_h[i][i];
		tmp2 = m_pclsGmres1_h[i][i+1];
		length = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		m_clsGmres1_c[i] = tmp1/length;
		m_clsGmres1_s[i] = -tmp2/length;

		// apply new Givens rotation to incomplete
		// (Q';tmp')*H(j+1) matrix i-th column
		// to get the complete (R(j+1);0)
		m_pclsGmres1_h[i][i] = m_clsGmres1_c[i] * tmp1 - m_clsGmres1_s[i] * tmp2;
		// this should be zero anyway
		//h(i+1,i) = c(i) * tmp2 + s(i) * tmp1;
		m_pclsGmres1_h[i][i+1] = 0;

		//
		// Compute recursively the first column of the (Q';tmp')
		// matrix of the QR factors of H(j+1); this is used to
		// get implicitly the norm of the residual and to compute
		// the solution once the norm is small enough
		//

		tmp1 = m_clsGmres1_g[i];
		m_clsGmres1_g[i] = m_clsGmres1_c[i] * tmp1;
		m_clsGmres1_g[i+1] = m_clsGmres1_s[i] * tmp1;

		// this is the norm of the residual
		// remark: || b - A*Qj*z || = normb * fabs(g[i+1])
		// (see "Matrix analysis and applied linear algebra", Meyers)
		// so fabs(g[i+1]) is already the norm of the residual
		// weighted by the norm of b (i.e the % difference)
		// However this is not true here, also verified
		// explicitly calculating the residual norm, once the solution
		// is found (see below), so the weighted formula is correct
		normr = fabs(m_clsGmres1_g[i+1]) / normb;

		//
		// if the norm of the residual is small enough,
		// exit loop (and return solution)
		//
		if(normr < gmresTol) {
			break;
		}
	}

	if(i >= SOLVE_GMRES_ITER_MAX || i >= size) {
		ErrMsg("\nError: not converging after %d iterations, norm of the residual is %.3f, while targeting %.3f\n", i, normr, gmresTol);
		// either case, since we ended the 'for' loop, 'i' has been incremented of one more;
		// to avoid breaking the end of the arrays here below, we must decrease 'i'
		i--;
	}

	if(m_clsGlobalVars.m_bDumpResidual == true) {
		LogMsg("%.3f ", normr);
	}
	LogMsg("\n");

	//
	// compute the solution, solving H(i)*z = normr * e1 for z using
	// the first row of the Q matrix and the R matrix of the QR
	// factorization of H; then the solution is x = Q(i)*z
	//

	iteration = i;

	// H(i)*z = normb*e1 gives Q*R*z = normb*e1 -> R*z = Q'*normb*e1
	// Since g(i) is already the first column of Q' multiplied by normb,
	// it is the term Q'*normb*e1, so init the solution
	m_clsGmres1_z = m_clsGmres1_g;

	// then solve R*z = Q'*normr*e1 by back substitution,
	// knowing that R is upper triangular
	for(i = iteration; i>=0; i--) {
		m_clsGmres1_z[i]= m_clsGmres1_z[i] / m_pclsGmres1_h[i][i];
		for(j = i-1; j>=0; j--) {
			m_clsGmres1_z[j] = m_clsGmres1_z[j] - m_pclsGmres1_h[i][j] * m_clsGmres1_z[i];
		}
	}


	// multiply z by Q(i) to get y, the solution
	// of the (possibly preconditioned) system
	for(i=0; i<size; i++) {
		for(j=0, tmp1 = 0; j<=iteration; j++) {
			tmp1 += m_pclsGmres1_q[j][i] * m_clsGmres1_z[j];
		}
		m_clsGmres1_y[i] = tmp1;
	}

	if(precondType == AUTOREFINE_PRECOND_NONE) {
		*x = m_clsGmres1_y;
	}
	else {
		// undo the preconditioner to get x, the final solution,
		// also considering initial vector in the solution
		//x = x0 + P * y;
		ComputePrecondVectFast(x, &m_clsGmres1_y, precondType);
		*x = m_clsGmres1_x0 + *x;
	}

	// debug : explicit calculation of residual vector
	// (is equal to normr*normb in this routine; this proves
	// that normr here is weighted on normb)
	//double residual2;
	//CLin_Vector b_approx(size);
	//m_clsMulthier.MultiplyMatByVec_fast(x, &b_approx);
	//residual2 = mod(*b - b_approx);

	LogMsg("\n");
	return 0;
}

//   Flexible Generalized Minimum Residual Method with preconditioner
//
//   x = gmresPrecond(A,b,P) attempt to solve the preconditioned
//   system A*P*y = b & x = P*y
//   Flexible means that the preconditioner can change at every iteration
//   (see "Iterative methods for sparse linear systems", Yousef Saad)
//
//   Enrico Di Lorenzo, 2004/07/06
//
//   'fast' refers to the direct vector multiplication loops used
//   instead of CLin_Vector class operator*, which requires
//   the creation of a temp object for every assignment operation
//
//   the 's' on the other hand refers to the scilab gmres correction,
//   to make convergence faster, that basically uses as residual
//   a weighted residual: normr = || b - Ax || / || b ||, that is
//   normr = fabs(g[i+1]) / normb;
//   and not
//   normr = fabs(g[i+1]);
//   see also the notes inline about this choice
//
//   the 'all' refers to the increase of speed not allocating memory
//   locally every time
//
int CSolveCap::gmresFlexPrecondSFastAll(CLin_Vector *b, CLin_Vector *x, double gmresTol)
{
	double normr, normb, tmp1, tmp2, length;
	long size, i, j, iteration, ii;
	bool ret;
	int retInt;

	// get system size
	size = b->size();
	// check consistency of size
//	ASSERT(size == *x.size());

	// if there is a preconditioner
	if(m_clsGlobalVars.m_ucPrecondType != AUTOREFINE_PRECOND_NONE) {

		// initial vector: having the approx inverted matrix, it is cheap to calculate the initial vector.
		//  Simply put, to solve the preconditioned system
		//  A * P * y = b and x = P * y, the best choice for y0 initial vector is b,
		//  since a good P should be as similar as possible to inv(A), so A * P ~= I
		//  So x0 = P * y0 = P * b
		ComputePrecondVectFast(&m_clsGmres_x0, b, m_clsGlobalVars.m_ucPrecondType);

		// r = b - A * x0;
		// matrix - vector multiplication
		retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, &m_clsGmres_x0);
		if(retInt != FC_NORMAL_END) {
			return retInt;
		}
		// and get r
		m_clsGmres_r = *b - m_clsGmres_v;
		normr = mod(m_clsGmres_r);
		normb = mod(*b);
	}
	else {
		// no preconditioner; assumes that initial vector x0 is an all zeros vector
		//
		m_clsGmres_r = *b;
		normr = mod(m_clsGmres_r);
		normb = normr;
	}

	// allocate and init first column of Q matrix
	m_pclsGmres_q[0] = CLin_Vector(m_clsGmres_r / normr);

	// init first element of first row of the
	// Q matrix of the QR factors of H; this should be 1,
	// using normr will produce the first row
	// multiplied by normr, which is used to compute the solution
	// Remark: this is not the Q matrix whose columns are q(i,j)
	m_clsGmres_g[0] = normr;

	//
	// if the norm of the residual is small enough,
	// initial vector is a good enough solution
	//
	if(normr / normb < gmresTol) {
		if(m_clsGlobalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE) {
			*x = CLin_Vector(size, 0.0);
		}
		else {
			*x = m_clsGmres_x0;
		}
		return FC_NORMAL_END;
	}

	LogMsg("GMRES Iteration: ");

	// start iteration
	for(i = 0; i < SOLVE_GMRES_ITER_MAX && i < size; i++) {

		if(g_bFCContinue == false) {
			return FC_USER_BREAK;
		}

		//
		// execute the i-th Arnoldi step, thus computing
		// the i-th column of the H (upper Hessemberg) matrix
		//

		if(m_clsGlobalVars.m_bDumpResidual == true) {
			LogMsg("%.3g ", normr);
		}
		LogMsg("%d ", i);

		// allocate only if not already pre-allocated
		if((long)m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_0_LEVEL] < i) {

			ASSERT((long)m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_0_LEVEL] == i-1);

			// allocate new column of h
			// (remark: since vectors are 0-based, must
			// allocate i+2 entries and not only i+1;
			// i.e. variable 'i' starts at 0!)
			ret = m_pclsGmres_h[i].newsize(i+2);
			if(ret == false) {
				return FC_OUT_OF_MEMORY;
			}
			g_clsMemUsage.m_ulGmresMem += (i + 2) * sizeof(double);
			// allocate new column of q
			ret = m_pclsGmres_q[i+1].newsize(size);
			if(ret == false) {
				return FC_OUT_OF_MEMORY;
			}
			g_clsMemUsage.m_ulGmresMem += size * sizeof(double);
			// allocate new column of zf
			ret = m_pclsGmres_zf[i].newsize(size);
			if(ret == false) {
				return FC_OUT_OF_MEMORY;
			}
			g_clsMemUsage.m_ulGmresMem += size * sizeof(double);

			m_uiGmresPrealloc[AUTOREFINE_HIER_PRE_0_LEVEL]++;
		}

		// compute new vector
		//
		if(m_clsGlobalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE) {
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, &m_pclsGmres_q[i]);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}
		else {
			// compute preconditioned vector
			ComputePrecondVectFast(&m_pclsGmres_zf[i], &m_pclsGmres_q[i], m_clsGlobalVars.m_ucPrecondType);
			// and use it in matrix - vector multiplication
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, &m_pclsGmres_zf[i]);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}

		m_clsGmres_w = m_clsGmres_v;

		// orthogonalize it
		for(j=0; j<=i; j++) {
			m_pclsGmres_h[i][j] = dot_prod(m_pclsGmres_q[j], m_clsGmres_w);
			for(ii=0; ii<(long)m_clsGmres_v.size(); ii++) {
				m_clsGmres_v[ii] = m_clsGmres_v[ii] - m_pclsGmres_h[i][j] * m_pclsGmres_q[j][ii];
			}
		}

		m_pclsGmres_h[i][i+1] = mod(m_clsGmres_v);

		for(ii=0; ii<(long)m_clsGmres_v.size(); ii++) {
			m_pclsGmres_q[i+1][ii] = m_clsGmres_v[ii] / (m_pclsGmres_h[i][i+1]);
		}

		//
		// rotate new vector (Givens rotations) to compute
		// the (R;0) matrix of the QR factors of H(j+1);
		// this is done updating (Q';tmp')*H(j) = (R(j);0)
		// to produce (Q';tmp')*H(j+1)=(R(j+1);0)
		// Remark: the first time the loop is executed,
		// this rotation does not happen
		//

		// apply all old rotations to new i-th column
		// of the H matrix to get Q*H from H, but for the
		// last two elements
		for(j=0; j<i; j++) {
			tmp1 = m_pclsGmres_h[i][j];
			tmp2 = m_pclsGmres_h[i][j+1];
			m_pclsGmres_h[i][j] = m_clsGmres_c[j] * tmp1 - m_clsGmres_s[j]* tmp2;
			m_pclsGmres_h[i][j+1] = m_clsGmres_c[j] * tmp2 + m_clsGmres_s[j] * tmp1;
		}

		// compute the new Givens rotation to annihilate h[i][i+1]
		tmp1 = m_pclsGmres_h[i][i];
		tmp2 = m_pclsGmres_h[i][i+1];
		length = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		m_clsGmres_c[i] = tmp1/length;
		m_clsGmres_s[i] = -tmp2/length;

		// apply new Givens rotation to incomplete
		// (Q';tmp')*H(j+1) matrix i-th column
		// to get the complete (R(j+1);0)
		m_pclsGmres_h[i][i] = m_clsGmres_c[i] * tmp1 - m_clsGmres_s[i] * tmp2;
		// this should be zero anyway
		//h(i+1,i) = c(i) * tmp2 + s(i) * tmp1;
		m_pclsGmres_h[i][i+1] = 0;

		//
		// Compute recursively the first column of the (Q';tmp')
		// matrix of the QR factors of H(j+1); this is used to
		// get implicitly the norm of the residual and to compute
		// the solution once the norm is small enough
		//

		tmp1 = m_clsGmres_g[i];
		m_clsGmres_g[i] = m_clsGmres_c[i] * tmp1;
		m_clsGmres_g[i+1] = m_clsGmres_s[i] * tmp1;

		// this is the norm of the residual
		// remark: || b - A*Qj*z || = normb * fabs(g[i+1])
		// (see "Matrix analysis and applied linear algebra", Meyers)
		// so fabs(g[i+1]) is already the norm of the residual
		// weighted by the norm of b (i.e the % difference)
		// However this is not true here, also verified
		// explicitly calculating the residual norm, once the solution
		// is found (see below), so the weighted formula is correct
		normr = fabs(m_clsGmres_g[i+1]) / normb;

		//
		// if the norm of the residual is small enough,
		// exit loop (and return solution)
		//
		if(normr < gmresTol) {
			break;
		}
	}

	if(i >= SOLVE_GMRES_ITER_MAX || i >= size) {
		ErrMsg("\nError: not converging after %d iterations, norm of the residual is %.3f, while targeting %.3f\n", i, normr, gmresTol);
		// either case, since we ended the 'for' loop, 'i' has been incremented of one more;
		// to avoid breaking the end of the arrays here below, we must decrease 'i'
		i--;
	}

	if(m_clsGlobalVars.m_bDumpResidual == true) {
		LogMsg("%.3f ", normr);
	}
	LogMsg("\n");

	//
	// compute the solution, solving H(i)*z = normr * e1 for z using
	// the first row of the Q matrix and the R matrix of the QR
	// factorization of H; then the solution is x = Q(i)*z
	//

	iteration = i;

	// H(i)*z = normb*e1 gives Q*R*z = normb*e1 -> R*z = Q'*normb*e1
	// Since g(i) is already the first column of Q' multiplied by normb,
	// it is the term Q'*normb*e1, so init the solution
	m_clsGmres_z = m_clsGmres_g;

	// then solve R*z = Q'*normr*e1 by back substitution,
	// knowing that R is upper triangular
	for(i = iteration; i>=0; i--) {
		m_clsGmres_z[i]= m_clsGmres_z[i] / m_pclsGmres_h[i][i];
		for(j = i-1; j>=0; j--) {
			m_clsGmres_z[j] = m_clsGmres_z[j] - m_pclsGmres_h[i][j] * m_clsGmres_z[i];
		}
	}


	if(m_clsGlobalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE) {
		// multiply z by Q(i) to get y, the solution
		// of the (possibly preconditioned) system
		for(i=0; i<size; i++) {
			for(j=0, tmp1 = 0; j<=iteration; j++) {
				tmp1 += m_pclsGmres_q[j][i] * m_clsGmres_z[j];
			}
			m_clsGmres_y[i] = tmp1;
		}

		*x = m_clsGmres_y;
	}
	else {
		// multiply by z by zf matrix to undo all preconditioners to get the final solution,
		// and sum x0 to consider also initial vector
		//x = x0 + zf * z;

		for(i = 0; i<size; i++) {
			for(j = 0, tmp1 = 0.0; j<=iteration; j++) {
				tmp1 += m_pclsGmres_zf[j][i] * m_clsGmres_z[j];
			}
			m_clsGmres_y[i] = tmp1;
		}

		*x = m_clsGmres_x0 + m_clsGmres_y;
	}

	// debug : explicit calculation of residual vector
	// (is equal to normr*normb in this routine; this proves
	// that normr here is weighted on normb)
	//double residual2;
	//CLin_Vector b_approx(size);
	//m_clsMulthier.MultiplyMatByVec_fast(x, &b_approx);
	//residual2 = mod(*b - b_approx);

	return FC_NORMAL_END;
}


//   Generalized Minimum Residual Method with preconditioner and initial vector
//
//   x = gmresPrecond(A,b,P,x0) attempt to solve the preconditioned
//   system A*P*y = b & x = P*y
//   using x0 as initial vector of the solution iteration
//
//   Enrico Di Lorenzo, 2009/08/17
//
//   'fast' refers to the direct vector multiplication loops used
//   instead of CLin_Vector class operator*, which requires
//   the creation of a temp object for every assignment operation
//
//   the 's' on the other hand refers to the scilab gmres correction,
//   to make convergence faster, that basically uses as residual
//   a weighted residual: normr = || b - Ax || / || b ||, that is
//   normr = fabs(g[i+1]) / normb;
//   and not
//   normr = fabs(g[i+1]);
//   see also the notes inline about this choice
//
//   the 'all' refers to the increase of speed not allocating memory
//   locally every time
//
int CSolveCap::gmresPrecondSFastAllX0(CLin_Vector *b, CLin_Vector *x, CLin_Vector *x0)
{
	// old, used with no x0 as initial vector
	//CLin_Vector r(*b);
	double normr, normb, tmp1, tmp2, length;
	long size, i, j, iteration, ii;
	bool ret;
	int retInt;

	// get system size
	size = b->size();
	// check consistency of size
//	ASSERT(size == *x.size());

	// r = b - A * x0;
	// matrix - vector multiplication
	retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, x0);
	if(retInt != FC_NORMAL_END) {
		return retInt;
	}
	// and get r
	m_clsGmres_r = *b - m_clsGmres_v;
	normr = mod(m_clsGmres_r);
	normb = mod(*b);


	// allocate and init first column of Q matrix
	m_pclsGmres_q[0] = CLin_Vector(m_clsGmres_r / normr);

	// init first element of first row of the
	// Q matrix of the QR factors of H; this should be 1,
	// using normr will produce the first row
	// multiplied by normr, which is used to compute the solution
	// Remark: this is not the Q matrix whose columns are q(i,j)
	m_clsGmres_g[0] = normr;

	//
	// if the norm of the residual is small enough,
	// initial vector is a good enough solution
	//
	if(normr / normb < m_clsGlobalVars.m_dGmresTol) {
		*x = *x0;
		return 0;
	}

	LogMsg("GMRES Iteration: ");

	// start iteration
	for(i = 0; i < SOLVE_GMRES_ITER_MAX && i < size; i++) {

		//
		// execute the i-th Arnoldi step, thus computing
		// the i-th column of the H (upper Hessemberg) matrix
		//

		if(m_clsGlobalVars.m_bDumpResidual == true) {
			LogMsg("%.3g ", normr);
		}
		LogMsg("%d ", i);

		// allocate new column of h
		// (remark: since vectors are 0-based, must
		// allocate i+2 entries and not only i+1;
		// i.e. variable 'i' starts at 0!)
		ret = m_pclsGmres_h[i].newsize(i+2);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += (i+2) * sizeof(double);
		// allocate new column of q
		ret = m_pclsGmres_q[i+1].newsize(size);
		if(ret == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulGmresMem += size * sizeof(double);

		// compute new vector
		//
		if(m_clsGlobalVars.m_ucPrecondType == AUTOREFINE_PRECOND_NONE) {
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_v, &m_pclsGmres_q[i]);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}
		else {
			// compute preconditioned vector
			ComputePrecondVectFast(&m_clsGmres_Pq, &m_pclsGmres_q[i], m_clsGlobalVars.m_ucPrecondType);
			// and use it in matrix - vector multiplication
			retInt = m_clsMulthier.MultiplyMatByVec_fast(&m_clsGmres_Pq, &m_clsGmres_v);
			if(retInt != FC_NORMAL_END) {
				return retInt;
			}
		}

		m_clsGmres_w = m_clsGmres_v;

		// orthogonalize it
		for(j=0; j<=i; j++) {
			m_pclsGmres_h[i][j] = dot_prod(m_pclsGmres_q[j], m_clsGmres_w);
			for(ii=0; ii<(long)m_clsGmres_v.size(); ii++) {
				m_clsGmres_v[ii] = m_clsGmres_v[ii] - m_pclsGmres_h[i][j] * m_pclsGmres_q[j][ii];
			}
		}

		m_pclsGmres_h[i][i+1] = mod(m_clsGmres_v);

		for(ii=0; ii<(long)m_clsGmres_v.size(); ii++) {
			m_pclsGmres_q[i+1][ii] = m_clsGmres_v[ii] / (m_pclsGmres_h[i][i+1]);
		}

		//
		// rotate new vector (Givens rotations) to compute
		// the (R;0) matrix of the QR factors of H(j+1);
		// this is done updating (Q';tmp')*H(j) = (R(j);0)
		// to produce (Q';tmp')*H(j+1)=(R(j+1);0)
		// Remark: the first time the loop is executed,
		// this rotation does not happen
		//

		// apply all old rotations to new i-th column
		// of the H matrix to get Q*H from H, but for the
		// last two elements
		for(j=0; j<i; j++) {
			tmp1 = m_pclsGmres_h[i][j];
			tmp2 = m_pclsGmres_h[i][j+1];
			m_pclsGmres_h[i][j] = m_clsGmres_c[j] * tmp1 - m_clsGmres_s[j]* tmp2;
			m_pclsGmres_h[i][j+1] = m_clsGmres_c[j] * tmp2 + m_clsGmres_s[j] * tmp1;
		}

		// compute the new Givens rotation to annihilate h[i][i+1]
		tmp1 = m_pclsGmres_h[i][i];
		tmp2 = m_pclsGmres_h[i][i+1];
		length = sqrt(tmp1 * tmp1 + tmp2 * tmp2);
		m_clsGmres_c[i] = tmp1/length;
		m_clsGmres_s[i] = -tmp2/length;

		// apply new Givens rotation to incomplete
		// (Q';tmp')*H(j+1) matrix i-th column
		// to get the complete (R(j+1);0)
		m_pclsGmres_h[i][i] = m_clsGmres_c[i] * tmp1 - m_clsGmres_s[i] * tmp2;
		// this should be zero anyway
		//h(i+1,i) = c(i) * tmp2 + s(i) * tmp1;
		m_pclsGmres_h[i][i+1] = 0;

		//
		// Compute recursively the first column of the (Q';tmp')
		// matrix of the QR factors of H(j+1); this is used to
		// get implicitly the norm of the residual and to compute
		// the solution once the norm is small enough
		//

		tmp1 = m_clsGmres_g[i];
		m_clsGmres_g[i] = m_clsGmres_c[i] * tmp1;
		m_clsGmres_g[i+1] = m_clsGmres_s[i] * tmp1;

		// this is the norm of the residual
		// remark: || b - A*Qj*z || = normb * fabs(g[i+1])
		// (see "Matrix analysis and applied linear algebra", Meyers)
		// so fabs(g[i+1]) is already the norm of the residual
		// weighted by the norm of b (i.e the % difference)
		// However this is not true here, also verified
		// explicitly calculating the residual norm, once the solution
		// is found (see below), so the weighted formula is correct
		normr = fabs(m_clsGmres_g[i+1]) / normb;

		//
		// if the norm of the residual is small enough,
		// exit loop (and return solution)
		//
		if(normr < m_clsGlobalVars.m_dGmresTol) {
			break;
		}
	}

	if(i >= SOLVE_GMRES_ITER_MAX || i >= size) {
		ErrMsg("\nError: not converging after %d iterations, norm of the residual is %.3f, while targeting %.3f\n", i, normr, m_clsGlobalVars.m_dGmresTol);
		// either case, since we ended the 'for' loop, 'i' has been incremented of one more;
		// to avoid breaking the end of the arrays here below, we must decrease 'i'
		i--;
	}

	if(m_clsGlobalVars.m_bDumpResidual == true) {
		LogMsg("%.3f ", normr);
	}
	LogMsg("\n");

	//
	// compute the solution, solving H(i)*z = normr * e1 for z using
	// the first row of the Q matrix and the R matrix of the QR
	// factorization of H; then the solution is x = Q(i)*z
	//

	iteration = i;

	// H(i)*z = normb*e1 gives Q*R*z = normb*e1 -> R*z = Q'*normb*e1
	// Since g(i) is already the first column of Q' multiplied by normb,
	// it is the term Q'*normb*e1, so init the solution
	m_clsGmres_z = m_clsGmres_g;

	// then solve R*z = Q'*normr*e1 by back substitution,
	// knowing that R is upper triangular
	for(i = iteration; i>=0; i--) {
		m_clsGmres_z[i]= m_clsGmres_z[i] / m_pclsGmres_h[i][i];
		for(j = i-1; j>=0; j--) {
			m_clsGmres_z[j] = m_clsGmres_z[j] - m_pclsGmres_h[i][j] * m_clsGmres_z[i];
		}
	}


	// multiply z by Q(i) to get y, the solution
	// of the (possibly preconditioned) system
	for(i=0; i<size; i++) {
		for(j=0, tmp1 = 0; j<=iteration; j++) {
			tmp1 += m_pclsGmres_q[j][i] * m_clsGmres_z[j];
		}
		m_clsGmres_y[i] = tmp1;
	}

	if(m_clsGlobalVars.m_ucPrecondType != AUTOREFINE_PRECOND_NONE) {
		// undo the preconditioner to get x, the final solution,
		// also considering initial vector in the solution
		//x = x0 + P * y;
		ComputePrecondVectFast(x, &m_clsGmres_y, m_clsGlobalVars.m_ucPrecondType);
		*x = *x0 + *x;
	}
	else {
		// consider initial vector in the solution
		//x = x0 + y;
		*x = *x0 + m_clsGmres_y;
	}

	// debug : explicit calculation of residual vector
	// (is equal to normr*normb in this routine; this proves
	// that normr here is weighted on normb)
	//double residual2;
	//CLin_Vector b_approx(size);
	//m_clsMulthier.MultiplyMatByVec_fast(x, &b_approx);
	//residual2 = mod(*b - b_approx);

	return 0;
}

void CSolveCap::ComputePrecondVectFast(CLin_Vector *Pq, CLin_Vector *q, unsigned char precondType)
{
	long i, j, k;
	unsigned long supPreIndex;
	unsigned int leavesNum;

	// if alternate preconditioner
	if( (precondType & AUTOREFINE_PRECOND_BLOCK) != 0 &&
	        (precondType & ~(AUTOREFINE_PRECOND_BLOCK)) != 0) {
		// select block preconditioner or the other preconditioner
		// (whatever it is) in turns
		if(m_ucAlternatePrecond % 2 != 0) {
			precondType = AUTOREFINE_PRECOND_BLOCK;
		}
		else {
			precondType &= ~(AUTOREFINE_PRECOND_BLOCK);
		}

		m_ucAlternatePrecond++;
	}

	if(precondType == AUTOREFINE_PRECOND_NONE) {

		//
		// no preconditioner
		//

		for(i=0; i<(long)(*Pq).size(); i++) {
			(*Pq)[i] = (*q)[i];
		}
	}
	else if(precondType == AUTOREFINE_PRECOND_JACOBI) {

		//
		// diagonal preconditioner
		//

		// Jacobi precond is the diagonal preconditioner, so
		// multiplication of precond by vector is a simple
		// element-by-element multiplication
		if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM) {
			for(i=0; i<(long)(*Pq).size(); i++) {
				// dividing for the 'm_clsSelfPotCoeff' is multiplying for the Jacobi precond
				(*Pq)[i] = (*q)[i] / m_clsMulthier.m_clsSelfPotCoeff[i];
				ASSERT(fabs(m_clsMulthier.m_clsSelfPotCoeff[i]) < 1E20);
				ASSERT(fabs(m_clsMulthier.m_clsSelfPotCoeff[i]) != 0.0);
				ASSERT(fabs((*Pq)[i]) < 1E20);
			}
		}
		else {
			// if complex permittivity, the '*q' and '*Pq' vectors are twice as long,
			// since they contain also the imaginary parts.
			// in this case, since the pot matrix is [R -C; C R] where R is the 'old'
			// potential matrix for the real case, must multiply by [J 0; 0 J],
			// where J is the Jacobi preconditioner
			for(i=0; i<(long)(*Pq).size()/2; i++) {
				// dividing for the 'm_clsSelfPotCoeff' is multiplying for the Jacobi precond
				(*Pq)[i] = (*q)[i] / m_clsMulthier.m_clsSelfPotCoeff[i];
				ASSERT(fabs(m_clsMulthier.m_clsSelfPotCoeff[i]) < 1E20);
				ASSERT(fabs(m_clsMulthier.m_clsSelfPotCoeff[i]) != 0.0);
				ASSERT(fabs((*Pq)[i]) < 1E20);
			}
			for(i=0, j=(long)(*Pq).size()/2; i<(long)(*Pq).size()/2; i++, j++) {
				// dividing for the 'm_clsSelfPotCoeff' is multiplying for the Jacobi precond
				(*Pq)[j] = (*q)[j] / m_clsMulthier.m_clsSelfPotCoeff[i];
				ASSERT(fabs((*Pq)[j]) < 1E20);
			}
		}
	}
	else if(precondType == AUTOREFINE_PRECOND_BLOCK) {

		//
		// block preconditioner
		//

		leavesNum = m_pucBlockPrecondDim[0];
		for(k=0; k<(long)m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); k += leavesNum) {
			leavesNum = m_pucBlockPrecondDim[k];
			for(i=0; i < (long)leavesNum; i++) {
				(*Pq)[i+k] = 0.0;
				for(j = 0; j < (long)leavesNum; j++) {
					(*Pq)[i+k] += m_pdBlockPrecond[i+k][j] * (*q)[j+k];
				}
			}
		}

	}
	else if(precondType == AUTOREFINE_PRECOND_SUPER) {

		//
		// super preconditioner
		//

		if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM) {
			// clear charges and potentials of super nodes
			for(i=0; i<(long)m_uiSupPreNum; i++) {
				m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dCharge = 0.0;
				m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dPotential = 0.0;
			}
			// collecting step, calculates the elements of the vector to be multiplied
			// by the super precond matrix starting from the leaves
			for(i=0; i<(long)m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++) {
				supPreIndex = m_puiSupPrecondIndex[i];
				// TBC warning: this is in the assumption that all leaves have the same dimensions, so we can average
				// without weighting. In general, should weight against areae.
				//m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dPotential += (*q)[i] / m_clsSupPrecondElements[supPreIndex].m_lLoLevLeavesNum;
				// This is with correct weighting. Remark: 'm_lLoLevLeavesNum' may not be useful anymore
				m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dPotential += (*q)[i] * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
			}
			// multiplication step, calculates the result of super preconditioner matrix multiplication
			for(i=0; i<(long)m_uiSupPreNum; i++) {
				for(k=0; k<(long)m_uiSupPreNum; k++) {
					// if auto potential, skip: must zero off-diagonal leaves-level elements, which would be positive
					// instead of negative, when extending the precond hi-lev matrix to lo-lev leaves (are based on
					// diagonal positive elements of Maxwell capacitance matrix of the super precond)
					if(i != k) {
						m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dCharge += m_pclsSupPrecondMtx[i][k] * m_clsSupPrecondElements[k].m_pclsHiLevLeaf->m_dPotential;
					}
				}
			}
			// distribution step, distributes the results obtained at super preconditioner level to the leaves
			for(i=0; i<(long)m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++) {
				supPreIndex = m_puiSupPrecondIndex[i];
				// TBC warning: this is in the assumption that all leaves have the same dimensions, so we can average
				// without weighting. In general, should weight against areae.
				//(*Pq)[i] = m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dCharge / m_clsSupPrecondElements[supPreIndex].m_lLoLevLeavesNum;
				// This is with correct weighting. Remark: 'm_lLoLevLeavesNum' may not be useful anymore
				(*Pq)[i] = m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dCharge * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
				// and add diagonal element, multiplied by the overall number of elements of the super block,
				// to account for zeroed off-diagonal elements (remark: super precond vector element value (*q)[i] should be
				// divided twice, if using compressed matrix multiplication algorithm. First time is in the collecting step,
				// second time is in the distribution step. Dividing only once accounts for the multiplication)
				// Why the zeroed off-diagonal elements of the super block are equal to the self-elements (in the hypotesis
				// of leaves all of the same size)? Because, see the comment above, in the multiplication step, where the
				// auto potential is skipped: in the super block there are all identical elements.
				// Another way to see the below operation: we divide the auto-capacitance of the super precond matrix
				// by the number of leaves of the block, then we multiply by q. Seen in this way, there is a straightforward
				// approach to account for different leave sizes, i.e. weighting on the areae
				//(*Pq)[i] += m_pclsSupPrecondMtx[supPreIndex][supPreIndex] * (*q)[i] / m_clsSupPrecondElements[supPreIndex].m_lLoLevLeavesNum;
				(*Pq)[i] += m_pclsSupPrecondMtx[supPreIndex][supPreIndex] * (*q)[i] * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
			}
		}
		else {
			// if complex permittivity, the '*q' and '*Pq' vectors are twice as long,
			// since they contain also the imaginary parts.
			// in this case, since the pot matrix is [R -C; C R] where R is the 'old'
			// potential matrix for the real case, must multiply by [S 0; 0 S],
			// where S is the Super preconditioner

			// first step for first half of '*q' and '*Pq'
			//

			// clear charges and potentials of super nodes
			for(i=0; i<(long)m_uiSupPreNum; i++) {
				m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dCharge = 0.0;
				m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dPotential = 0.0;
			}
			// collecting step, calculates the elements of the vector to be multiplied
			// by the super precond matrix starting from the leaves
			for(i=0; i<(long)m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++) {
				supPreIndex = m_puiSupPrecondIndex[i];
				m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dPotential += (*q)[i] * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
			}
			// multiplication step, calculates the result of super preconditioner matrix multiplication
			for(i=0; i<(long)m_uiSupPreNum; i++) {
				for(k=0; k<(long)m_uiSupPreNum; k++) {
					// if auto potential, skip
					if(i != k) {
						m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dCharge += m_pclsSupPrecondMtx[i][k] * m_clsSupPrecondElements[k].m_pclsHiLevLeaf->m_dPotential;
					}
				}
			}
			// distribution step, distributes the results obtained at super preconditioner level to the leaves
			for(i=0; i<(long)m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++) {
				supPreIndex = m_puiSupPrecondIndex[i];
				(*Pq)[i] = m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dCharge * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
				// and add diagonal element
				(*Pq)[i] += m_pclsSupPrecondMtx[supPreIndex][supPreIndex] * (*q)[i] * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
			}

			// second step for second half of '*q' and '*Pq'
			//

			// clear charges and potentials of super nodes
			for(i=0; i<(long)m_uiSupPreNum; i++) {
				m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dCharge = 0.0;
				m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dPotential = 0.0;
			}
			// collecting step, calculates the elements of the vector to be multiplied
			// by the super precond matrix starting from the leaves
			for(i=0, j=(long)(*Pq).size()/2; i<(long)m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++, j++) {
				supPreIndex = m_puiSupPrecondIndex[i];
				m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dPotential += (*q)[j] * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
			}
			// multiplication step, calculates the result of super preconditioner matrix multiplication
			for(i=0; i<(long)m_uiSupPreNum; i++) {
				for(k=0; k<(long)m_uiSupPreNum; k++) {
					// if auto potential, skip
					if(i != k) {
						m_clsSupPrecondElements[i].m_pclsHiLevLeaf->m_dCharge += m_pclsSupPrecondMtx[i][k] * m_clsSupPrecondElements[k].m_pclsHiLevLeaf->m_dPotential;
					}
				}
			}
			// distribution step, distributes the results obtained at super preconditioner level to the leaves
			for(i=0, j=(long)(*Pq).size()/2; i<(long)m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++, j++) {
				supPreIndex = m_puiSupPrecondIndex[i];
				(*Pq)[j] = m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->m_dCharge * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
				// and add diagonal element
				(*Pq)[j] += m_pclsSupPrecondMtx[supPreIndex][supPreIndex] * (*q)[j] * m_pfSupPrecondAreae[i] / m_clsSupPrecondElements[supPreIndex].m_pclsHiLevLeaf->GetDimension();
			}
		}

	}
	else {

		//
		// hierarchical preconditioner
		//

		// multiplication step, calculates the result of super preconditioner matrix multiplication
		// using gmres (hierarchical super preconditioner)
		m_clsMulthier.SetInteractionLevel(AUTOREFINE_HIER_PRE_1_LEVEL);
		// leave only the lower level precond (if any)
		precondType &= ~(AUTOREFINE_PRECOND_HIER);
		gmresPrecondSFastAllUpper(q, Pq, m_clsGlobalVars.m_dHierPreGmresTol, precondType);
		// reset lower level
		m_clsMulthier.SetInteractionLevel(AUTOREFINE_HIER_PRE_0_LEVEL);
	}
}


#ifdef DEBUG_DUMP_POT

// warning: dumps the transpose of the pot matrix
void CSolveCap::DebugDumpPotMtxAndIndex()
{
	long condindex, dielindex;
	StlAutoCondDeque::iterator itc1;
	FILE *fp;

	// dump pot matrix in Scilab format
	DebugDumpPotMtx();

	// dump a vector indexing the panels belonging to each conductor and to each dielectric

	fp = fopen("pot_mat_index.txt", "w");

	ASSERT(fp != NULL);

	// scan every conductor to mark the panels belonging to each conductor
	condindex = 1;
	dielindex = -1;
	for(itc1 = m_clsMulthier.m_stlConductors.begin();
	        itc1 != m_clsMulthier.m_stlConductors.end(); itc1++) {

		// if conductor is real conductor and not dielectric
		if( (*itc1)->m_bIsDiel == false) {

			DumpCondPanelPerm(*itc1, (CAutoElement*)((*itc1)->m_uTopElement.m_pTopPanel), fp, condindex);

			condindex++;
		}
		else {

			DumpDielPanelPerm(*itc1, (CAutoElement*)((*itc1)->m_uTopElement.m_pTopPanel), fp, dielindex);

			dielindex--;
		}
	}

	fclose(fp);
}

// warning: dumps the transpose of the pot matrix
// remark: to use the matrix for dielectric i/f, must invert, multiply by V, and *multiply the entries
// by outperm* before summing the charge densities!
void CSolveCap::DebugDumpPotMtx()
{
	CLin_Vector p, c;
	long i, j, panelnum, mtxsize;
	unsigned char tmpType;
	FILE *fp;

    panelnum = m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL);

	if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
        mtxsize = panelnum;
	}
	else {
	    mtxsize = panelnum * 2;
	}

    p.newsize(mtxsize);
    c.newsize(mtxsize);
    c = 0.0;

	// dump transpose of pot matrix in SciLab format

	// avoid 'playing' with the fast multiplication for 2D solvers
	tmpType = g_ucSolverType;
//    g_ucSolverType = SOLVERGLOBAL_3DSOLVER;

	fp = fopen("pot_mat.txt", "w");

	ASSERT(fp != NULL);

	// to dump pot matrix, we find each column of A
	// by multiplying by c vectors with all zeros
	// but the i-th element, set to 1
	// p = A*c is then the i-th column of A
	for(i=0; i<mtxsize; i++) {

		if(i>0) {
			c[i-1] = 0;
		}
		c[i] = 1.0;

		m_clsMulthier.MultiplyMatByVec_fast(&p, &c);

		// actually dump to file in scilab matrix format
		// warning: dumped in this way, this is the transpose of A
		for(j=0; j<mtxsize; j++) {

			// once every 20 values, split row
			//if(j%20 = 0) {
			//	fprintf(fp, "\n");
			//}
			fprintf(fp, "%f ", p[j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	// restore solver type
	g_ucSolverType = tmpType;

}

void CSolveCap::DumpCondPanelPerm(CAutoConductor *cond, CAutoElement *panel, FILE *fp, long condindex)
{
	// visit the tree

	// termination condition
	if (panel->IsLeaf() == true) {
		// output index and outperm for each panel belonging to the 'itc1' conductor
		fprintf(fp, "%f %f %f %f %f %f\n", (float)condindex, (float)cond->m_dSurfOutperm[panel->m_ucDielIndex][0], (float)cond->m_dSurfOutperm[panel->m_ucDielIndex][1], 0.0f, 0.0f, panel->GetDimension());
	}
	else {
		// otherwise, scan the tree
		DumpCondPanelPerm(cond, panel->m_pLeft, fp, condindex);
		DumpCondPanelPerm(cond, panel->m_pRight, fp, condindex);
	}
}

void CSolveCap::DumpDielPanelPerm(CAutoConductor *cond, CAutoElement *panel, FILE *fp, long condindex)
{
	// visit the tree

	// termination condition
	if (panel->IsLeaf() == true) {
		// output index and outperm for each panel belonging to the 'itc1' conductor
		fprintf(fp, "%f %f %f %f %f %f\n", (float)condindex, (float)cond->m_dOutperm[0], (float)cond->m_dOutperm[1], (float)cond->m_dInperm[0], (float)cond->m_dInperm[1], panel->GetDimension());
	}
	else {
		// otherwise, scan the tree
		DumpDielPanelPerm(cond, panel->m_pLeft, fp, condindex);
		DumpDielPanelPerm(cond, panel->m_pRight, fp, condindex);
	}
}

#endif //DEBUG_DUMP_POT

#ifdef DEBUG_DUMP_OTHER

// warning: dumps the transpose of the precond matrix
void CSolveCap::DebugDumpPrecondMtx()
{
	CLin_Vector p(m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL));
	CLin_Vector c(m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL), 0.0);
	long i, j;
	FILE *fp;

	fp = fopen("precond_mat.txt", "w");

	ASSERT(fp != NULL);

	// to dump pot matrix, we find each column of A
	// by multiplying by c vectors with all zeros
	// but the i-th element, set to 1
	// p = A*c is then the i-th column of A

	for(i=0; i<m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++) {

		if(i>0) {
			c[i-1] = 0;
		}
		c[i] = 1.0;

		ComputePrecondVectFast(&p, &c, m_clsGlobalVars.m_ucPrecondType);

		// actually dump to file in scilab matrix format
		// warning: dumped in this way, this is the transpose of A
		for(j=0; j<m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); j++) {

			// once every 20 values, split row
			//if(j%20 = 0) {
			//	fprintf(fp, "\n");
			//}
			fprintf(fp, "%e ", p[j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}
#endif //DEBUG_DUMP_OTHER

#ifdef DEBUG_DUMP_UNCOMP_POT
// warning: dumps the transpose of the native, uncompressed pot matrix
// remark: valid ONLY for uniform dielectric
void CSolveCap::DebugDumpUncomprPotMtx()
{
	StlAutoCondDeque::iterator itc1;
	int isPotValid;
	double potestim1, potestim2, potestRe, potestIm;
	long i, j;
	FILE *fp;

	if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {

		// create leaf panels array and fill it up
		m_pLeafPanels = new CAutoPanel*[m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL)];

		m_pLeafPanIndex = 0;
		for(itc1 = m_clsMulthier.m_stlConductors.begin(); itc1 != m_clsMulthier.m_stlConductors.end(); itc1++) {
			DebugRecurseCopyPanels((*itc1)->m_uTopElement.m_pTopPanel);
		}


		fp = fopen("pot_mat_unc.txt", "w");

		ASSERT(fp != NULL);

		// to dump pot matrix, we find each column of A
		// by multiplying by c vectors with all zeros
		// but the i-th element, set to 1
		// p = A*c is then the i-th column of A

		for(i=0; i<m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++) {

			// actually dump to file in scilab matrix format
			// warning: dumped in this way, this is the transpose of A
			for(j=0; j<m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); j++) {

				if(i!=j) {
					isPotValid = m_clsMulthier.PotEstimateOpt(m_pLeafPanels[j], m_pLeafPanels[i], potestim1, potestim2);
				}
				else {
					m_clsMulthier.SelfPotential(m_pLeafPanels[i], &potestRe, &potestIm);
					potestim1 = potestRe;
					isPotValid = 0;
				}

				if(isPotValid != 0) {
					ASSERT(0);
				}

				// once every 20 values, split row
				//if(j%20 = 0) {
				//	fprintf(fp, "\n");
				//}
				fprintf(fp, "%f ", potestim1);
			}
			fprintf(fp, "\n");
		}

		fclose(fp);

		delete m_pLeafPanels;
	}
	// 2D case
	else {
		// create leaf panels array and fill it up
		m_pLeafSegments = new CAutoSegment*[m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL)];

		m_pLeafPanIndex = 0;
		for(itc1 = m_clsMulthier.m_stlConductors.begin(); itc1 != m_clsMulthier.m_stlConductors.end(); itc1++) {
			DebugRecurseCopyPanels((*itc1)->m_uTopElement.m_pTopSegment);
		}


		fp = fopen("pot_mat_unc.txt", "w");

		ASSERT(fp != NULL);

		// to dump pot matrix, we find each column of A
		// by multiplying by c vectors with all zeros
		// but the i-th element, set to 1
		// p = A*c is then the i-th column of A

		for(i=0; i<m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); i++) {

			// actually dump to file in scilab matrix format
			// warning: dumped in this way, this is the transpose of A
			for(j=0; j<m_clsMulthier.GetPanelNum(AUTOREFINE_HIER_PRE_0_LEVEL); j++) {

				if(i != j) {
					isPotValid = m_clsMulthier.PotEstimateOpt(m_pLeafSegments[j], m_pLeafSegments[i], potestim1, potestim2);
				}
				else {
					m_clsMulthier.SelfPotential(m_pLeafSegments[i], &potestRe, &potestIm);
					potestim1 = potestRe;
					isPotValid = 0;
				}

				if(isPotValid != 0) {
					ASSERT(0);
				}

				// once every 20 values, split row
				//if(j%20 = 0) {
				//	fprintf(fp, "\n");
				//}
				fprintf(fp, "%f ", potestim1);
			}
			fprintf(fp, "\n");
		}

		fclose(fp);

		delete m_pLeafSegments;
	}
}

// recursively copy leaf panels in an array, for debug
// remark: for how the function works, panels must have already
// been indexed (counted)
void CSolveCap::DebugRecurseCopyPanels(CAutoPanel* panel)
{
	if(panel->IsLeaf() == true) {
		m_pLeafPanels[m_pLeafPanIndex] = panel;
		m_pLeafPanIndex++;
	}
	else {
		// if not a leaf panel, go into left and right sub-trees
		DebugRecurseCopyPanels((CAutoPanel*)panel->m_pLeft);
		DebugRecurseCopyPanels((CAutoPanel*)panel->m_pRight);
	}
}

// recursively copy leaf panels in an array, for debug
// remark: for how the function works, panels must have already
// been indexed (counted)
void CSolveCap::DebugRecurseCopyPanels(CAutoSegment* panel)
{
	if(panel->IsLeaf() == true) {
		m_pLeafSegments[m_pLeafPanIndex] = panel;
		m_pLeafPanIndex++;
	}
	else {
		// if not a leaf panel, go into left and right sub-trees
		DebugRecurseCopyPanels((CAutoSegment*)panel->m_pLeft);
		DebugRecurseCopyPanels((CAutoSegment*)panel->m_pRight);
	}
}

#endif //DEBUG_DUMP_UNCOMP_POT

