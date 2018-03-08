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


// AutoRefine.cpp : autorefine class
//
// Automatic mesh refinement routines for capacitance extraction
// with N-body methods (Fast Multipole, Hierarchical)
// Enrico Di Lorenzo, 2002/08/21

#include "../stdafx.h"

#include <string>
#include <iostream>
#include <time.h>
// for _isnan() and _finite()
#include <float.h>
// for _chdir()
//#include <direct.h>
// for openmp
#include "omp.h"

// for wxGetFreeMemory()
#include <wx/utils.h> 

#include "SolverGlobal.h"

#include "Autorefine.h"

#include "Geometry/Triangulate.h"
#include "Geometry/Operation3D.h"

// for debug only
#include "SphereGen.h"

// to test against calcp()
//#include "mulStruct.h"

// max path length
#define AR_MAX_PATH	    4096

// maximum number of time we try to generate a new temp file name
// (actually, the number is AR_MAX_TEMP_IDS * AR_MAX_TEMP_IDS)
#define AR_MAX_TEMP_IDS     10

#define AUTOREFINE_MAX_LINE_LEN		512

// for quality triangulation of quadrialteral panels
#define AUTOREFINE_MIN_ANGLE	20

// tolerance in the difference of opposite sides lengths of a quad panel
#define AUTOREFINE_QUAD_SIDE_TOL 0.2
// cosine of the tolerance in the rectangular angle of a quad panel
#define AUTOREFINE_QUAD_ANGLE_TOL 0.1

// tolerance in the difference between two relative dielctric constants across
// a dielectric-dielectric interface
#define AUTOREFINE_REL_DIEL_TOL		0.001

// type of operations when creating panels
#define AUTOREFINE_SIMPLE_CREATE_PANEL			0
#define AUTOREFINE_TRIANGULATE_CREATE_PANEL		1


// extension for temporary files; l = links, p = potentials
#define FC_LINK_TMP_FILE_PREFIX _T("frcl")
#define FC_POT_TMP_FILE_PREFIX _T("frcp")

// cosinus of the minumum angle to declare a triangle 'thin'
// (0.99619469809174553229501040247389 means 5 degrees)
// (0.99862953475457387378449205843944 means 3 degrees)
#define AUTOREFINE_COS_MIN	0.99619469809174553229501040247389

// minimum segment length
#define AUTOREFINE_MIN_LEN  1E-10

// z coordinate to be used when dumping 2D models as 3D FastCap-compatible geometry files
#define AUTOREFINE_2D_ZCOORD    1.0

// maximum distance of any two points in a 2D model, for rescaling;
// must be less than 1.0 to avoid -log(r) to become negative,
// and possibly not so near to 1.0 to avoid getting a zero, that in spite
// of all the possible eps would never lead to a refinement
#define AUTOREFINE_MAX_2D_DIM   0.5

// init static vars
unsigned long CAutoRefine::m_ulTempFileID = 1;

CAutoRefine::CAutoRefine()
{
	unsigned char i;

	// init number of panels
	for(i=0; i<AUTOPANEL_MAX_NUM_OF_HIERARCHIES; i++) {
		m_ulPanelNum[i] = 0;
		m_ulNodeNum[i] = 0;
		m_dPotCoeffLinks[i] = NULL;
		m_pdPanelPtrLinks[i] = NULL;
		m_ulLinkChunkNum[i] = 0;
	}

	m_ulUniqueLinkChunkIDs = NULL;
	m_ulUniquePotChunkIDs = NULL;
	m_fGlobalCharges = NULL;
	m_pucDielIndex = NULL;
	m_pNodes = NULL;

	// init seed used in guessing unique file IDs
	srand(time(NULL));
}

CAutoRefine::~CAutoRefine()
{
	Clean(AUTOREFINE_DEALLMEM_ALL, m_clsGlobalVars);
}
/*
// reset the structures to the status after BuildSuperHierarchy()
void CAutoRefine::CleanOnlyUpToSupHie()
{
	unsigned int i;

	// delete only refined panels
	DeleteRefinedPanels();

	// delete link arrays and reset other vars
	for(i=0; i<AUTOPANEL_MAX_NUM_OF_HIERARCHIES; i++) {
		DeleteLinkArray(i);
		m_ulPanelNum[i] = 0;
		m_ulLinkChunkNum[i] = 0;
	}
}
*/
void CAutoRefine::Clean(int command, CAutoRefGlobalVars globalVars)
{
	unsigned int i;
	unsigned long totChunksNum, k;
	wxFileName tmpFileName;
    StlAutoCondDeque::iterator itc;

	if( (command == AUTOREFINE_DEALLMEM_ALL) ||
	        (command == AUTOREFINE_DEALLMEM_AT_END && globalVars.m_bKeepCharge == false) ||
	        (command == AUTOREFINE_DEALLMEM_AT_START && globalVars.m_bRefineCharge == false && globalVars.m_bKeepMesh == false) ) {

		// delete panel tree and conductor list
		DeletePanelsAndConductors();
		// delete nodes vector
		if(m_pNodes != NULL) {
			delete m_pNodes;
			m_pNodes = NULL;
		}
		// delete permittivity vector
		if(m_pucDielIndex != NULL) {
			delete m_pucDielIndex;
			m_pucDielIndex = NULL;
		}
		g_clsMemUsage.m_ulPanelsMem = 0;
		g_clsMemUsage.m_ulCondMem = 0;
	}

	if( (command == AUTOREFINE_DEALLMEM_ALL) ||
	        (command == AUTOREFINE_DEALLMEM_AT_END && globalVars.m_bKeepCharge == false ) ||
	        (command == AUTOREFINE_DEALLMEM_AT_START  && globalVars.m_bRefineCharge == false) ||
	        (command == AUTOREFINE_DEALLMEM_CLEAN_CHARGES)) {

		if(m_fGlobalCharges != NULL) {
			delete m_fGlobalCharges;
			m_fGlobalCharges = NULL;
		}
		g_clsMemUsage.m_ulChargesMem = 0;
	}

	// delete self-potential array
	m_clsSelfPotCoeff.destroy();
	m_clsImgSelfPotCoeff.destroy();


	// delete temporary chunks files and the chunk IDs array

	if(m_ulBlocksNum > 1) {
		// total number of chunks
		totChunksNum = m_ulBlocksNum * m_ulLinkChunkNum[m_ucInteractionLevel];
		// retrieve file name from chunk IDs
		for(k=0; k< totChunksNum; k++) {
			// get the temp file path and name from the ID
			PortableGetTempFileName(wxT(""), FC_LINK_TMP_FILE_PREFIX, m_ulUniqueLinkChunkIDs[k], tmpFileName);
			// delete the file
			if( remove((const char*)(tmpFileName.GetLongPath())) != 0 ) {
				ErrMsg("Error: cannot delete the temporary file '%s'\n", (const char*)(tmpFileName.GetLongPath()));
				ErrMsg("       Temporary files must be manually removed.\n");
				ErrMsg("       Continuing with memory deallocation\n");
			}
			// get the temp file path and name from the ID
			PortableGetTempFileName(wxT(""), FC_POT_TMP_FILE_PREFIX, m_ulUniquePotChunkIDs[k], tmpFileName);
			// delete the file
			if( remove((const char*)(tmpFileName.GetLongPath())) != 0 ) {
				ErrMsg("Error: cannot delete the temporary file '%s'\n", (const char*)(tmpFileName.GetLongPath()));
				ErrMsg("       Temporary files must be manually removed.\n");
				ErrMsg("       Continuing with memory deallocation\n");
			}
		}
		// deallocate the memory
		delete m_ulUniqueLinkChunkIDs;
		m_ulUniqueLinkChunkIDs = NULL;
		delete m_ulUniquePotChunkIDs;
		m_ulUniquePotChunkIDs = NULL;
	}

	// delete link arrays and reset other vars
	for(i=0; i<AUTOPANEL_MAX_NUM_OF_HIERARCHIES; i++) {
		DeleteLinkArray(i);
		m_ulLinkChunkNum[i] = 0;
	}

	// clean up memory usage information
	g_clsMemUsage.m_ulLinksMem = 0;
}

void CAutoRefine::DeleteLinkArray(unsigned int level)
{
	unsigned long j;

	// delete interaction link arrays
	if(m_dPotCoeffLinks[level] != NULL) {
		for(j=0; j<m_ulLinkChunkNum[level]; j++) {
			if(m_dPotCoeffLinks[level][j] != NULL) {
				delete m_dPotCoeffLinks[level][j];
			}
		}
		delete m_dPotCoeffLinks[level];
		m_dPotCoeffLinks[level] = NULL;
	}
	if(m_pdPanelPtrLinks[level] != NULL) {
		for(j=0; j<m_ulLinkChunkNum[level]; j++) {
			if(m_pdPanelPtrLinks[level][j] != NULL) {
				delete m_pdPanelPtrLinks[level][j];
			}
		}
		delete m_pdPanelPtrLinks[level];
		m_pdPanelPtrLinks[level] = NULL;
	}
}

// Automatic refinement (mesh generator), given an input file
// describing a geometry in FastCap format
//
// Composed by two functions: AutoRefinePanels(), AutoRefineLinks()
//
// Version 14.0
//
// This version:
// 1. uses Stroud numerical quadrature (mutual2.m)
// 2. stores panels interactions and potential values
//    into 'panels' array structure
// 3. uses original Appel refinement algorithm
// 4. uses potestim * panel area as refinement criteria, instead
//    of potestim * panel diameter, to avoid infinite recursion
// 5. normalize panel areas to avoid refinement dependance from
//    overall geometrical dimensions
// 6. refines also panels belonging to the same conductor, to accurately
//    compute auto capacitances
//
//
// fileinname	is the name of the input file
//
// fileoutname	is the name of the output file
//
// pepsauto     auto potential tolerance
//
// pepsmutual   mutual potential tolerance
//
// the function returns the total number of patches composing the mesh
//
int CAutoRefine::AutoRefinePanels(CAutoRefGlobalVars globalVars, unsigned char interactLevel)
{
	unsigned long i;
	int ret;
	StlAutoCondDeque::iterator itc1, itc2;
	clock_t start, finish;

// debug
//	int panelSize = sizeof(CAutoPanel);
//	int pointerSize = sizeof(CAutoConductor*);

	m_clsGlobalVars = globalVars;

	m_ucInteractionLevel = interactLevel;

	// de-normalize epsilon
//	m_clsGlobalVars.m_dEps = m_clsGlobalVars.m_dEps / FOUR_PI_TIMES_E0;

	// clean up link array at this level
	DeleteLinkArray(m_ucInteractionLevel);

	// init binary tree structure-related vars
	m_iLevel = 0;
	m_iMaxLevel = 0;
	m_ulLinksNum[m_ucInteractionLevel] = 0;

	// init complexity computation vars (for algorithm analysis)
	m_ulNumofpotest = 0;
//	m_ulNumofFastPotest = 0;


#ifdef DEBUG_DUMP_BASIC
#ifdef DEBUG_DUMP_OTHER

	for(i=0; i<32; i++) {
		for(j=0; j<32; j++) {
			m_iaLinksBtwLevels[i][j] = 0;
			m_iaPotestBtwLevels[i][j] = 0;
		}
		m_laLinkNumStatistics[i] = 0;
		m_laPanelsPerLevel[i] = 0;
		m_laLinksPerLevel[i] = 0;
	}
#endif
#endif

	// start timer
	start = clock();

	// if refinement is going to be based on the charge computed in
	// a previous run, copy the charge densitied into the panels,
	// while resetting the link counters for each panel
	if(m_clsGlobalVars.m_bRefineCharge == true) {
		// signal to init charge vars
		m_bInitCharges = true;

		m_ulCountPanelNum = 0;
		// scan all conductor groups
		for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {

			if(g_bFCContinue == false) {
				return FC_USER_BREAK;
			}

			// copy charges into leaf panels
			CopyCharges((*itc1)->m_uTopElement.m_pTopPanel);
		}

		// refine only the panels whose charge density is in the last 'x' percentile
		m_dMidSigma = (m_dMaxSigma - m_dMinSigma) * 0.5 + m_dMinSigma;
	}

	if(m_clsGlobalVars.m_bKeepMesh == false) {

		// init number of panels and nodes at this level
		m_ulPanelNum[m_ucInteractionLevel] = 0;
		m_ulNodeNum[m_ucInteractionLevel] = 0;

		// if the user wants to refine the geometry based on the charge information
		if(m_clsGlobalVars.m_bRefineCharge == true) {

			// scan all conductor groups
			for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {

				if(g_bFCContinue == false) {
					return FC_USER_BREAK;
				}

				// Discretize geometry
				ret = Discretize((*itc1)->m_uTopElement.m_pTopPanel);
				if(ret !=  FC_NORMAL_END) {
					return ret;
				}
			}
		}
		else {
			// otherwise create an adaptive mesh, in function of:
			// 1) curvature of conductors (self) 2) proximity (mutual)

			// initialize max mesh eps
			m_dMaxMeshEps = 0.0;

			for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
				// refine against every other conductor
				for(itc2 = itc1 + 1; itc2 != m_stlConductors.end(); itc2++) {
					if(g_bFCContinue == false) {
						return FC_USER_BREAK;
					}
					// refine
					if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
						DiscretizeMutual((*itc1)->m_uTopElement.m_pTopPanel, (*itc2)->m_uTopElement.m_pTopPanel, false);
					}
					else if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
						DiscretizeMutual((*itc1)->m_uTopElement.m_pTopSegment, (*itc2)->m_uTopElement.m_pTopSegment, false);
					}
					else {
						ASSERT(false);
					}
				}
			}
			for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
				if(g_bFCContinue == false) {
					return FC_USER_BREAK;
				}
				// remark: discretize self should be called only AFTER all mutual discretization
				// has been performed. This is to allow the mutual routines to refine the panels,
				// in case, so we don't have the case of conductors made of a single panel (e.g.
				// square gnd plane), not discretized, unless this discretization was not needed
				// for mutual refinement
				m_iLevel = 0;
				if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
					DiscretizeSelf((*itc1)->m_uTopElement.m_pTopPanel);
				}
				else if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
					DiscretizeSelf((*itc1)->m_uTopElement.m_pTopSegment);
				}
				else {
					ASSERT(false);
				}
			}
		}
	}
	else {
		// scan all conductor groups
		for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {

			if(g_bFCContinue == false) {
				return FC_USER_BREAK;
			}

			// perform panel initializations (to reset required panel parameters)
			(*itc1)->m_uTopElement.m_pTopPanel->InitElementsTree();
		}
	}

	// if global charges array was allocated, clear it (in this point of the code, the array is no more useful:
	// if required, it should have already been used in Discretize() function)
	Clean(AUTOREFINE_DEALLMEM_CLEAN_CHARGES, m_clsGlobalVars);

	// allocate nodes array (will be populated in RecurseIndex() called in AutoRefineLinks() )
	SAFENEW_ARRAY_RET(CAutoElement*, m_pNodes, m_ulNodeNum[m_ucInteractionLevel], g_clsMemUsage.m_ulPanelsMem)

	// allocate array of global charges, used for refinement (if requested, and if on bottom level)
	if(m_clsGlobalVars.m_bKeepCharge == true && m_ucInteractionLevel == AUTOREFINE_HIER_PRE_0_LEVEL) {

		// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
		SAFENEW_ARRAY_RET(float, m_fGlobalCharges, m_ulPanelNum[m_ucInteractionLevel], g_clsMemUsage.m_ulChargesMem)
		// clear the array
		for(i=0; i<m_ulPanelNum[m_ucInteractionLevel]; i++) {
			m_fGlobalCharges[i] = 0.0;
		}
	}

	// dry-scan all conductor groups to calculate # of links per panel and per supernode,
	// plus the self-potentials.
	// Must do it here and not in AutoRefineLinks() because we need to know how many links
	// will be generated before actually computing them, loosing time in case the refinement
	// is not enough (in automatic mode)
	m_bComputeLinks = false;
	ret = ComputeLinks();


	finish = clock();
	m_fDurationDiscretize = (float)(finish - start) / CLOCKS_PER_SEC;

	return ret;
}

int CAutoRefine::AutoRefineLinks(CAutoRefGlobalVars globalVars)
{
	int ret;
	StlAutoCondDeque::iterator itc1, itc2;
	double start, finish;
	unsigned long j, k, totChunksNum;
	// 'i' must have signed integral type due to MS OpenMP limitation (does not accept unsigned)
	long i;
	wxLongLong mem_Potest, mem_pCharge, mem_LinksTotal, mem_AvailVirtual, mem_MaxAllocVirtual, mem_AllocVirtual;
	wxLongLong freeDiskBytes;
	wxFileName tmpFileName, tmpFNObj;
	bool goOutOfCore, retBool;
	unsigned long uniqueFileID;
	unsigned long linkIndex, chunk, block, nodeIndex, nodeBlockEnd, linksPerBlock;


	// start timer
	start = omp_get_wtime();

	// index links and index panels
	m_ulBasePanelNum = 0;
	m_ulBaseLinksNum = 0;
	m_ulCountNodeNum = 0;
	m_bPopulateNodeArray = true;
	for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
		m_ulCountPanelNum = 0;
		RecurseIndex((*itc1)->m_uTopElement.m_pTopPanel);
		m_ulBasePanelNum += m_ulCountPanelNum;
		// and record number of leaf panels in the parent conductor
		(*itc1)->m_ulLeafPanelNum = m_ulCountPanelNum;
	}
	ASSERT( m_ulBasePanelNum == m_ulPanelNum[m_ucInteractionLevel]);
	ASSERT( m_ulBaseLinksNum == m_ulLinksNum[m_ucInteractionLevel]);
	ASSERT( m_ulCountNodeNum == m_ulNodeNum[m_ucInteractionLevel]);


	// calculate the memory needed for the potential estimates
	mem_Potest = ((wxLongLong)m_ulLinksNum[m_ucInteractionLevel]) * ((wxLongLong)sizeof(double));
	// calculate the memory needed for the pointer to interacting panel charge
	mem_pCharge = ((wxLongLong)m_ulLinksNum[m_ucInteractionLevel]) * ((wxLongLong)sizeof(double*));
	// total required memory
	mem_LinksTotal = mem_Potest + mem_pCharge;

	// get available memory
	mem_AvailVirtual = (wxLongLong) wxGetFreeMemory();
// debug
//mem_AvailVirtual = 50000000;

	goOutOfCore = false;
	if(mem_AvailVirtual.ToLong() == -1L) {
		ErrMsg("Error: cannot retrieve the information about the free memory quantity\n");
		ErrMsg("       Cannot go out-of-core, continuing in-core\n");
	}
	else if( mem_LinksTotal.ToDouble() * globalVars.m_dOutOfCoreRatio < mem_AvailVirtual.ToDouble() ) {
		// no need to go out of core
	}
	else {
		LogMsg("Estimated memory required for storing panel interaction links is: %lu Mbytes\n", (mem_LinksTotal/G_MEGABYTE).ToLong());
		LogMsg("Available free memory left is: %lu Mbytes\n", (mem_AvailVirtual/G_MEGABYTE).ToLong());
		LogMsg("Estimated links memory is more than %.0f%% of the available free memory. Going out-of-core.\n", 100.0f / globalVars.m_dOutOfCoreRatio);

		// try to create a temporary file, to understand if we can access the temp directory
		// and if we have the write access rights to it
		uniqueFileID = PortableGetTempFileName(wxT(""), FC_LINK_TMP_FILE_PREFIX, 0, tmpFileName);

		if(uniqueFileID == 0) {
			// get the temp path
			tmpFNObj.Assign(wxFileName::GetTempDir());

			// remark: use long path name version (if possible), i.e. without '~'
			ErrMsg("Error: cannot access the temporary folder '%s' to store out-of-core files\n", (const char*)(tmpFNObj.GetLongPath()));
			ErrMsg("       Cannot go out-of-core, continuing in-core\n");
		}
		else {
			// remove the created temp file
			remove(tmpFileName.GetFullPath().c_str());

			// get only the path part, w/o the file name
			tmpFNObj.Assign(tmpFileName.GetPath());

			LogMsg("Out-of-core temporary files folder: '%s'\n", (const char*)tmpFNObj.GetLongPath());

			// and check the free space
			// This function returns the total number of bytes and number of free bytes on the disk
			//  containing the directory path (it should exist).
			retBool = wxGetDiskSpace(tmpFileName.GetPath(), &freeDiskBytes);

			if(retBool == false) {
				ErrMsg("Error: cannot retrieve the free space on disk for directory %s\n", (const char*)tmpFNObj.GetLongPath());
				ErrMsg("       Cannot go out-of-core, continuing in-core\n");
			}
			else {
				if( freeDiskBytes < mem_LinksTotal ) {
					ErrMsg("Error: free space in the temporary directory on disk is less than the required value of %lu Mbytes\n", (mem_LinksTotal/G_MEGABYTE).ToLong());
					ErrMsg("       Cannot go out-of-core, trying to continue in-core\n");
				}

				// ok we can go out-of-core
				goOutOfCore = true;
			}
		}
	}

	if(goOutOfCore == true) {
		// 'mem_AvailVirtual / globalVars.m_dOutOfCoreRatio' is the max ram memory block we decided to allocate for storing the chunks
		mem_MaxAllocVirtual.Assign(mem_AvailVirtual.ToDouble() / globalVars.m_dOutOfCoreRatio);
		// let's calculate how many chunks fit in this block size
		// We use the wxLongLong division, which is available either with native 64 bits integers or through wxWidgets
		// implementation (remark: the class documentation does not report some operators, including '/', '<', etc.).
		// Integers division results in quotient and remainder. Ignoring the remainder results always in truncation
		// (approximation to the lower integer).
		// Here we want truncation: the number of chunks kept in memory at the same time must be such
		// to require in any case an amount of memory below the total max ram memory block allocated
		// for the purpose ('mem_MaxAllocVirtual')
		// Remark: all values are casted to wxLongLong before operations happen, to avoid implicit conversions that may lead
		// to errors, e.g. if 'AUTOREFINE_LINK_CHUNK_SIZE' is considered long, sizeof is considered long, but their product
		// needs more than 32 bits to be stored; not likely, but to be avoided anyway
		m_ulLinkChunkNum[m_ucInteractionLevel] = ( mem_MaxAllocVirtual / ( ((wxLongLong)AUTOREFINE_LINK_CHUNK_SIZE) * ((wxLongLong)(sizeof(double) + sizeof(double*))) ) ).ToLong();
		if(m_ulLinkChunkNum[m_ucInteractionLevel] == 0) {
			ErrMsg("Error: available free memory is not enough to allocate any chunk\n");
			ErrMsg("       Cannot go out-of-core, terminating process\n");
			return FC_CANNOT_GO_OOC;
		}
		// let's calculate the real size of the memory block (it has been rounded down to fit entire chunks)
		mem_AllocVirtual = ((wxLongLong)AUTOREFINE_LINK_CHUNK_SIZE) * ((wxLongLong)(sizeof(double) + sizeof(double*))) * ((wxLongLong)m_ulLinkChunkNum[m_ucInteractionLevel]);
		// now let's calculate how many blocks we need
		m_ulBlocksNum = (unsigned long) (mem_LinksTotal / mem_AllocVirtual).ToLong();
		// result of the division is rounded towards zero (truncated). So if there is a remainder, need one block more
		if(mem_LinksTotal % mem_AllocVirtual != 0) {
			m_ulBlocksNum++;
		}
		// remark: if number of blocks is one, no need to go out-of-core
	}
	else {
		// calculate how many chunks are needed
		m_ulLinkChunkNum[m_ucInteractionLevel] = m_ulLinksNum[m_ucInteractionLevel] / AUTOREFINE_LINK_CHUNK_SIZE;
		// result of the division is rounded towards zero. So if there is a remainder, need one chunk more
		if(m_ulLinksNum[m_ucInteractionLevel] % AUTOREFINE_LINK_CHUNK_SIZE != 0) {
			m_ulLinkChunkNum[m_ucInteractionLevel]++;
		}
		// and only one block (all in-core)
		m_ulBlocksNum = 1;
	}

	// allocate temporary file names for all chunks to be stored out-of-core
	//

	// allocate array of chunk file name IDs (only if we need to go out-of-core)
	//
	if(m_ulBlocksNum > 1) {
		// total number of chunks
		totChunksNum = m_ulBlocksNum * m_ulLinkChunkNum[m_ucInteractionLevel];
		// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
		SAFENEW_ARRAY_RET(unsigned long, m_ulUniqueLinkChunkIDs, totChunksNum, g_clsMemUsage.m_ulLinksMem)
		SAFENEW_ARRAY_RET(unsigned long, m_ulUniquePotChunkIDs, totChunksNum, g_clsMemUsage.m_ulLinksMem)
		// create IDs
		for(i=0, k=0; i< (long)m_ulBlocksNum; i++) {
			for(j=0; j< m_ulLinkChunkNum[m_ucInteractionLevel]; j++) {

				uniqueFileID = PortableGetTempFileName(wxT(""), FC_LINK_TMP_FILE_PREFIX, 0, tmpFileName);
				if(uniqueFileID == 0) {
					ErrMsg("Error: cannot create the temporary file '%s', to store out-of-core files\n", (const char*)tmpFileName.GetFullPath());
					ErrMsg("       Cannot go out-of-core, stopping process\n");

					return FC_CANNOT_GO_OOC;
				}
				m_ulUniqueLinkChunkIDs[k] = uniqueFileID;

				uniqueFileID = PortableGetTempFileName(wxT(""), FC_POT_TMP_FILE_PREFIX, 0, tmpFileName);
				if(uniqueFileID == 0) {
					ErrMsg("Error: cannot create the temporary file '%s', to store out-of-core files\n", (const char*)tmpFileName.GetFullPath());
					ErrMsg("       Cannot go out-of-core, stopping process\n");

					return FC_CANNOT_GO_OOC;
				}
				m_ulUniquePotChunkIDs[k] = uniqueFileID;

				k++;
			}
		}
	}


	// allocate interaction link arrays
	//


	// allocate pointers to chunks
	// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
	SAFENEW_ARRAY_RET(double*, m_dPotCoeffLinks[m_ucInteractionLevel], m_ulLinkChunkNum[m_ucInteractionLevel], g_clsMemUsage.m_ulLinksMem)
	SAFENEW_ARRAY_RET(CAutoElement**, m_pdPanelPtrLinks[m_ucInteractionLevel], m_ulLinkChunkNum[m_ucInteractionLevel], g_clsMemUsage.m_ulLinksMem)

//	LogMsg("Memory information before starting the allocation of link chunks\n");
//	DumpMemoryInfo();

	// and allocate chunks
	for(k = 0; k < m_ulLinkChunkNum[m_ucInteractionLevel]; k++) {
		// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
		SAFENEW_ARRAY_RET(double, m_dPotCoeffLinks[m_ucInteractionLevel][k], AUTOREFINE_LINK_CHUNK_SIZE, g_clsMemUsage.m_ulLinksMem)
		SAFENEW_ARRAY_RET(CAutoElement*, m_pdPanelPtrLinks[m_ucInteractionLevel][k], AUTOREFINE_LINK_CHUNK_SIZE, g_clsMemUsage.m_ulLinksMem)

//		LogMsg("Memory information after allocation of link chunk #%d\n", k);
//		DumpMemoryInfo();
	}


	// allocate self-potential array (kept separate to be used also as Jacobi precond,
	// in particular not to have it out-of-core in case we decided to go OOC)
	retBool = m_clsSelfPotCoeff.newsize(m_ulPanelNum[AUTOREFINE_HIER_PRE_0_LEVEL]);
	if(retBool == false) {
		return FC_OUT_OF_MEMORY;
	}
	g_clsMemUsage.m_ulLinksMem += m_ulPanelNum[AUTOREFINE_HIER_PRE_0_LEVEL] * sizeof(double);
	// and initialize the array to all zeroes (to be able to understand if the self-potential
	// has already been calculated or not, so not to calculate it again over multiple iterations
	// over ComputeLinks() )
	m_clsSelfPotCoeff = 0.0;
	if( globalVars.m_ucHasCmplxPerm != AUTOREFINE_REAL_PERM) {
		// allocate the imaginary part of the self-potential array
		retBool = m_clsImgSelfPotCoeff.newsize(m_ulPanelNum[AUTOREFINE_HIER_PRE_0_LEVEL]);
		if(retBool == false) {
			return FC_OUT_OF_MEMORY;
		}
		g_clsMemUsage.m_ulLinksMem += m_ulPanelNum[AUTOREFINE_HIER_PRE_0_LEVEL] * sizeof(double);
	}
	m_clsImgSelfPotCoeff = 0.0;


	// allocate array for the conductor panels permittivity values
	// SAFENEW_ARRAY_RET(TYPE, VAR, LEN, MEM)
	SAFENEW_ARRAY_RET(unsigned char, m_pucDielIndex, m_ulPanelNum[AUTOREFINE_HIER_PRE_0_LEVEL], g_clsMemUsage.m_ulPanelsMem)

	//
	// first pass in computing links. Just find and store the links, not their values
	//

	for(m_ulCurrBlock=0; m_ulCurrBlock < m_ulBlocksNum; m_ulCurrBlock++) {

		// actually scan all conductor groups and calculate interactions, stored in the links arrays
		m_bComputeLinks = true;
		ret = ComputeLinks();

		if(ret !=  FC_NORMAL_END) {
			return ret;
		}

		// and save the links out-of-core, if requested
		if(m_ulBlocksNum > 1) {

			ret = SaveLinks(false);

			if(ret !=  FC_NORMAL_END) {
				return ret;
			}
		}

		// reset start and end link pointers for each panel,
		// this is needed to run again ComputeLinks()
		// The pointers are not reset for the last iteration,
		// since in this case ComputeLinks() won't be called any more

		if(m_ulCurrBlock < m_ulBlocksNum - 1) {
			// index links and index panels
			m_ulBasePanelNum = 0;
			m_ulBaseLinksNum = 0;
			m_bPopulateNodeArray = false;
			for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
				m_ulCountPanelNum = 0;
				RecurseIndex((*itc1)->m_uTopElement.m_pTopPanel);
				m_ulBasePanelNum += m_ulCountPanelNum;
			}
			ASSERT( m_ulBasePanelNum == m_ulPanelNum[m_ucInteractionLevel]);
			ASSERT( m_ulBaseLinksNum == m_ulLinksNum[m_ucInteractionLevel]);
		}

	}
	// set current block to the latest actually in memory
	m_ulCurrBlock--;

	//
	// second pass in computing links. Now use the links list to know which are the interacting panels,
	// and calculate interactions (this is the long part, but in this way it can be done in parallel)
	// This is similar in structure to what done in ComputePanelPotentials_2fast()
	//

	// 'm_ulLinkChunkNum' is the number of chunks per block
	linksPerBlock = AUTOREFINE_LINK_CHUNK_SIZE * m_ulLinkChunkNum[m_ucInteractionLevel];

	// scan all links, in blocks
	for(linkIndex=0, nodeIndex=0; linkIndex<GetLinksNum(); linkIndex += linksPerBlock)  {

		// load data from Mass Memory device (out of core)
		//
		// determine in which block the current chunk is ('m_ulLinkChunkNum' is the number of chunks per block)
		chunk = linkIndex / AUTOREFINE_LINK_CHUNK_SIZE;
		block = chunk / m_ulLinkChunkNum[m_ucInteractionLevel];
		// if not in current block
		if(block != m_ulCurrBlock) {
			ret = LoadLinks(block, false);
			if(ret != FC_NORMAL_END) {
				return ret;
			}
		}

		// scan the nodes whose links are in the current block (of course they must be ordered)
		// but first we need to find where to stop (to avoid scanning all nodes every time)
		// note that 'nodeBlockEnd' will indicate one node more (both in cases we exit the for statement
		// normally or due to the break)
		for(nodeBlockEnd=nodeIndex; nodeBlockEnd<m_ulNodeNum[AUTOREFINE_HIER_PRE_0_LEVEL]; nodeBlockEnd++) {
			// if going into next chunk, we stop here
			if(m_pNodes[nodeBlockEnd]->m_ulLinkIndexEnd[m_ucInteractionLevel] >= linkIndex+linksPerBlock) {
				nodeBlockEnd++;
				break;
			}
		}

		// now calculate coefficients of potential

		#pragma omp parallel for
		for(i=nodeIndex; i<(long)nodeBlockEnd; i++) {
			unsigned long localLinkIndex, localChunk, localPosInChunk;
			CAutoElement *element2;
			double potestim1;

			// perform summation
			for(localLinkIndex = m_pNodes[i]->m_ulLinkIndexStart[m_ucInteractionLevel]; localLinkIndex < m_pNodes[i]->m_ulLinkIndexEnd[m_ucInteractionLevel]; localLinkIndex++) {

				// some of the links could be outside the boundary of the chunk, either on the left
				// (but not for the first node) or on the right. In this case, skip
				if(localLinkIndex >= linkIndex) {
					localChunk = localLinkIndex / AUTOREFINE_LINK_CHUNK_SIZE;
					localPosInChunk = localLinkIndex % AUTOREFINE_LINK_CHUNK_SIZE;
					// adjust chunk to position within the current block
					localChunk -= m_ulLinkChunkNum[m_ucInteractionLevel] * m_ulCurrBlock;
					// if still within the boundary
					if(localChunk < m_ulLinkChunkNum[m_ucInteractionLevel]) {
						// interacting panels have already been identified
						element2 = m_pdPanelPtrLinks[m_ucInteractionLevel][localChunk][localPosInChunk];
						// calculate coefficient of potential
						PotEstimateOpt(m_pNodes[i], element2, potestim1);
						// and store it
						m_dPotCoeffLinks[m_ucInteractionLevel][localChunk][localPosInChunk] = potestim1;
					}
				}
			}
		}
		// and move to next nodeIndex position
		nodeIndex = nodeBlockEnd - 1;

		// save the links out-of-core, if requested
		if(m_ulBlocksNum > 1) {

			ret = SaveLinks();

			if(ret !=  FC_NORMAL_END) {
				return ret;
			}
		}

	}

	// check time
	finish = omp_get_wtime();
	m_fDurationRefine = (float)(finish - start);

	return ret;
}

// can save both links and coefficients of potential, or only links
int CAutoRefine::SaveLinks(bool saveAlsoPot)
{
	FILE *stream;
	unsigned long j, k, start;
	wxFileName tmpFileName;
	int numwritten, ret;

	ret = FC_NORMAL_END;

	for(j=0; j< m_ulLinkChunkNum[m_ucInteractionLevel]; j++) {
		// current start of chunks block
		start  = m_ulCurrBlock * m_ulLinkChunkNum[m_ucInteractionLevel];

		PortableGetTempFileName(wxT(""), FC_LINK_TMP_FILE_PREFIX, m_ulUniqueLinkChunkIDs[start + j], tmpFileName);
		// open the file
		stream = fopen(tmpFileName.GetLongPath(), "wb");
		if( stream == NULL ) {
			ErrMsg("Error: cannot open the temporary file '%s' for writing\n", (const char*)tmpFileName.GetLongPath());
			return FC_FILE_ERROR;
		}

		// and write current chunk to disk
		for(k = 0; k < AUTOREFINE_LINK_CHUNK_SIZE; k++) {
			numwritten = fwrite( &m_pdPanelPtrLinks[m_ucInteractionLevel][j][k], sizeof( CAutoPanel* ), 1, stream );
			if(numwritten != 1) {
				ErrMsg("Error: cannot write to temporary out-of-core file, stopping the process\n");
				j = m_ulLinkChunkNum[m_ucInteractionLevel];
				ret = FC_FILE_ERROR;
				break;
			}
		}

		fclose(stream);

		if(saveAlsoPot == true && j < m_ulLinkChunkNum[m_ucInteractionLevel]) {
			PortableGetTempFileName(wxT(""), FC_POT_TMP_FILE_PREFIX, m_ulUniquePotChunkIDs[start + j], tmpFileName);
			// open the file
			stream = fopen(tmpFileName.GetLongPath(), "wb");
			if( stream == NULL ) {
				ErrMsg("Error: cannot open the temporary file '%s' for writing\n", (const char*)tmpFileName.GetLongPath());
				return FC_FILE_ERROR;
			}

			// and write current chunk to disk
			for(k = 0; k < AUTOREFINE_LINK_CHUNK_SIZE; k++) {
				numwritten = fwrite( &m_dPotCoeffLinks[m_ucInteractionLevel][j][k], sizeof( double ), 1, stream );
				if(numwritten != 1) {
					ErrMsg("Error: cannot write to temporary out-of-core file, stopping the process\n");
					j = m_ulLinkChunkNum[m_ucInteractionLevel];
					ret = FC_FILE_ERROR;
					break;
				}
			}

			fclose(stream);
		}
	}

	return ret;
}

int CAutoRefine::LoadLinks(unsigned long block, bool loadAlsoPot)
{
	FILE *stream;
	unsigned long j, k, start;
	wxFileName tmpFileName;
	int numwritten, ret;

	ret = FC_NORMAL_END;

	// set the current memory block to the block actually being loaded from disk
	m_ulCurrBlock = block;

	for(j=0; j< m_ulLinkChunkNum[m_ucInteractionLevel]; j++) {
		// current start of chunks block
		start  = m_ulCurrBlock * m_ulLinkChunkNum[m_ucInteractionLevel];

		PortableGetTempFileName(wxT(""), FC_LINK_TMP_FILE_PREFIX, m_ulUniqueLinkChunkIDs[start + j], tmpFileName);
		// open the file
		stream = fopen(tmpFileName.GetFullPath().c_str(), "rb");
		if( stream == NULL ) {
			ErrMsg("Error: cannot open the temporary file '%s' for reading\n", (const char*)tmpFileName.GetFullPath());
			return FC_FILE_ERROR;
		}

		// and read the current chunk from disk
		for(k = 0; k < AUTOREFINE_LINK_CHUNK_SIZE; k++) {
			numwritten = fread( &m_pdPanelPtrLinks[m_ucInteractionLevel][j][k], sizeof( CAutoPanel* ), 1, stream );
			if(numwritten != 1) {
				if(feof(stream)) {
					ErrMsg("Error: unexpected end-of-file of temporary out-of-core file, stopping the process\n");
				}
				else {
					ErrMsg("Error: cannot read from temporary out-of-core file, stopping the process\n");
				}
				ErrMsg("Error: cannot read from temporary out-of-core file, stopping the process\n");
				j = m_ulLinkChunkNum[m_ucInteractionLevel];
				ret = FC_FILE_ERROR;
				break;
			}
		}

		fclose(stream);

		if(loadAlsoPot == true && j < m_ulLinkChunkNum[m_ucInteractionLevel]) {
			PortableGetTempFileName(wxT(""), FC_POT_TMP_FILE_PREFIX, m_ulUniquePotChunkIDs[start + j], tmpFileName);
			// open the file
			stream = fopen(tmpFileName.GetFullPath().c_str(), "rb");
			if( stream == NULL ) {
				ErrMsg("Error: cannot open the temporary file '%s' for reading\n", (const char*)tmpFileName.GetFullPath());
				return FC_FILE_ERROR;
			}

			// and read the current chunk from disk
			for(k = 0; k < AUTOREFINE_LINK_CHUNK_SIZE; k++) {
				numwritten = fread( &m_dPotCoeffLinks[m_ucInteractionLevel][j][k], sizeof( double ), 1, stream );
				if(numwritten != 1) {
					if(feof(stream)) {
						ErrMsg("Error: unexpected end-of-file of temporary out-of-core file, stopping the process\n");
					}
					else {
						ErrMsg("Error: cannot read from temporary out-of-core file, stopping the process\n");
					}
					ErrMsg("Error: cannot read from temporary out-of-core file, stopping the process\n");
					j = m_ulLinkChunkNum[m_ucInteractionLevel];
					ret = FC_FILE_ERROR;
					break;
				}
			}

			fclose(stream);
		}
	}

	return ret;
}

// This function mimics GetTempFileName() by MS but in a portable way using wxWidgets functions.
// The function creates a name for a temporary file. The file name is the concatenation of specified path
// and prefix strings, a hexadecimal string formed from a specified integer ('id'), and the .tmp extension.
// The specified integer can be nonzero, in which case, the function creates the file name but does not create the file.
// If you specify zero for the integer, the function creates a unique file name and creates the file in the specified directory.
// If it was not possible to create an unique file name after a certain number of trials, the function returns zero.
// Different from GetTempFileName(), if the path 'dirName' is empty, the function uses the standard temp directory, that
// is prepended as path to the 'tempfullfilename'.
// Remark: if 'dirName' is specified, it must end with the path separator (to this purpose you can use the wxFileName::GetPathWithSep()
// function), otherwise the last directory is taken as a file name
unsigned long CAutoRefine::PortableGetTempFileName(const wxString &dirName, const wxString &prefix, unsigned long id, wxFileName &filenameobj)
{
	wxString dir, name;
	int i, j;
	FILE *fout;

	// clean up the file name object
	filenameobj.Clear();

	// extract the directory specified by 'dirName'
	wxFileName::SplitPath(dirName, &dir, &name, NULL);

	// if dirName specified an empty path
	if (dir.empty() == true) {
		dir = wxFileName::GetTempDir();
	}

	if(id != 0) {
		// must not create the file, only the filename, and do not test for uniqueness
		//
		// compose the name
		name = prefix + wxString::Format(wxT("%lX"), id) + wxT(".tmp");
		// merge the path and the filename
		filenameobj.Assign(dir, name);
	}
	else {
		id = m_ulTempFileID;
		for(id=m_ulTempFileID, i=0; i<AR_MAX_TEMP_IDS; i++) {
			for(j=0; j<AR_MAX_TEMP_IDS; j++, id++ ) {
				// compose the name
				name = prefix + wxString::Format(wxT("%lX"), id) + wxT(".tmp");
				// merge the path and the filename
				filenameobj.Assign(dir, name);
				if(filenameobj.FileExists() == false) {
					break;
				}
			}
			if(filenameobj.FileExists() == false) {
				break;
			}
			else {
				// get a new base ID chosen randomly
				id = (unsigned long) rand();
			}
		}
		// store the next id for the next time the function will be called
		m_ulTempFileID = id + 1;
		// if in the end we did NOT found the name, assign zero to 'id' to signal it to the caller
		if(filenameobj.FileExists() == true) {
			id = 0;
		}
		else {
			// now create the file to:
			// 1) test we have the right priviledges to write in this directory
			// 2) 'book' the file name to prevent others from creating the same file name,
			//    before the caller could have the chance to actually use it
			fout = fopen(filenameobj.GetFullPath().c_str(), "w");
			if(fout == NULL) {
				id = 0;
			}
			else {
				fclose(fout);
			}
		}
	}

	return id;
}

void CAutoRefine::DumpMemoryInfo()
{
	wxLongLong freeMem;

	// get available memory
	freeMem = (wxLongLong) wxGetFreeMemory();

	LogMsg("Available virtual memory: %lu kilobytes\n", (freeMem/G_KILOBYTE).ToLong());

    // Non-portable method
	/*
	    MEMORYSTATUS memstate;

	    GlobalMemoryStatus(&memstate);

	    LogMsg("Total virtual memory: %lu kilobytes\n", memstate.dwTotalVirtual/G_KILOBYTE);
	    LogMsg("Available virtual memory: %lu kilobytes\n", memstate.dwAvailVirtual/G_KILOBYTE);
	*/

}

void CAutoRefine::CopyCharges(CAutoElement *panel)
{
	double sigma;

	if(panel->IsLeaf() == true) {

		ASSERT((panel->m_ucType & AUTOPANEL_IS_SUPER_NODE) != AUTOPANEL_IS_SUPER_NODE);
		ASSERT(m_ulCountPanelNum < m_ulPanelNum[m_ucInteractionLevel]);

		// save charges in panels as charge densities

		// densities
		sigma = m_fGlobalCharges[m_ulCountPanelNum] / panel->GetDimension();
// debug
//		sigma = m_fGlobalCharges[m_ulCountPanelNum];
//		sigma = m_fGlobalCharges[m_ulCountPanelNum] * panel->GetDimension();

		panel->m_dCharge = sigma;
		panel->m_dPotential = sigma;

		if(m_bInitCharges == true) {
			m_dMaxSigma = sigma;
			m_dMinSigma = sigma;
			m_bInitCharges = false;
		}
		else {
			// find max and min charges
			if(sigma > m_dMaxSigma) {
				m_dMaxSigma = sigma;
			}
			if(sigma < m_dMinSigma) {
				m_dMinSigma = sigma;
			}
		}

		// increase panel number
		m_ulCountPanelNum++;
	}
	else {
		// then call recursively the routine for each child
		CopyCharges(panel->m_pLeft);
		CopyCharges(panel->m_pRight);
	}
}

int CAutoRefine::ComputeLinks()
{
	StlAutoCondDeque::iterator itc1, itc2;
	int ret;

	ret = FC_NORMAL_END;

	for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
		// for mutual capacitance,
		// refine against every other conductor
		for(itc2 = itc1 + 1; itc2 != m_stlConductors.end(); itc2++) {
			if(g_bFCContinue == false) {
				return FC_USER_BREAK;
			}
			// refine
			if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
				ret = RefineMutual((*itc1)->m_uTopElement.m_pTopPanel, (*itc2)->m_uTopElement.m_pTopPanel);
			}
			else if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
				ret = RefineMutual((*itc1)->m_uTopElement.m_pTopSegment, (*itc2)->m_uTopElement.m_pTopSegment);
			}
			else {
				ASSERT(false);
			}
		}
		// compute auto coefficients of potential
		SetCurrentConductor( *itc1 );
		if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
			ret = RefineSelf((*itc1)->m_uTopElement.m_pTopPanel);
		}
		else if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
			ret = RefineSelf((*itc1)->m_uTopElement.m_pTopSegment);
		}
		else {
			ASSERT(false);
		}
	}

	return ret;
}

void CAutoRefine::RecurseIndex(CAutoElement *panel)
{
	long numLinks;

	if(panel->IsLeaf() == true) {

		ASSERT((panel->m_ucType & AUTOPANEL_IS_SUPER_NODE) != AUTOPANEL_IS_SUPER_NODE);

		// store panel index number
		panel->m_lIndex[m_ucInteractionLevel] = m_ulBasePanelNum+m_ulCountPanelNum;
		// increase panel number
		m_ulCountPanelNum++;

		// count number of children panels
		panel->m_lNumOfChildren = 1;
	}
	else {
		// then call recursively RecurseIndex() for each child
		RecurseIndex(panel->m_pLeft);
		RecurseIndex(panel->m_pRight);

		// and count number of children panels
		panel->m_lNumOfChildren = panel->m_pLeft->m_lNumOfChildren + panel->m_pRight->m_lNumOfChildren;
	}

	if(m_bPopulateNodeArray == true ) {
		m_pNodes[m_ulCountNodeNum] = panel;
		m_ulCountNodeNum++;
	}

	// remember for the moment the number of links of this panel
	// (first time we call the routine, 'panel->m_ulLinkIndexStart' is zero and 'panel->m_ulLinkIndexEnd'
	// is used for counting; next times, we need to perform the difference to have the correct number)
	numLinks = panel->m_ulLinkIndexEnd[m_ucInteractionLevel] - panel->m_ulLinkIndexStart[m_ucInteractionLevel];
	// if at bottom hierarchical level (AUTOREFINE_HIER_PRE_0_LEVEL) and in a leaf,
	// reserve first position for the self potential
	// (see how RefineSelf() used this parameter)
	//if(m_ucInteractionLevel == AUTOREFINE_HIER_PRE_0_LEVEL && panel->IsLeaf() == true) {
	//	panel->m_ulLinkIndexEnd[m_ucInteractionLevel] = m_ulBaseLinksNum + 1;
	//}
	//else {
	panel->m_ulLinkIndexEnd[m_ucInteractionLevel] = m_ulBaseLinksNum;
	//}
	// store link index number
	panel->m_ulLinkIndexStart[m_ucInteractionLevel] = m_ulBaseLinksNum;
	// increase link number
	m_ulBaseLinksNum += numLinks;
}

// setting of the pointer to the conductor currently processed,
// needed for SelfPotential() routine, to get information about
// outer and inner dielectric constant (stored in the CAutoConductor structure)
void CAutoRefine::SetCurrentConductor(CAutoConductor *m_pCurrCond)
{
	m_pCurrentConductor = m_pCurrCond;
}

void CAutoRefine::OutputFastCapFile(std::string fileinname, std::string suffix, CLin_Vector *condCharges)
{
	FILE *foutlst;
	StlAutoCondDeque::iterator itc;
	std::string basefilename, fileoutname, condfilename, condfilenopath;
	std::string condfilenameext, condfilenopathext, ext;
	std::string::size_type pos;
	int i;
	wxString dielIndexStr;

	m_pLocalCondCharge = condCharges;

	// build file out name from fileinname, stripping old extension if any
	// and adding 'suffix' + ".lst"

	basefilename = fileinname;
	pos = basefilename.rfind(".");
	// if there was an extension, strip it, otherwise nothing
	if(pos != basefilename.npos ) {
		basefilename.resize(pos);
	}
	fileoutname = basefilename;
	fileoutname += "_";
	fileoutname += suffix;
	fileoutname += ".lst";

	foutlst = fopen(fileoutname.c_str(), "w");

	if(foutlst == NULL) {
        ErrMsg("  Warning: cannot open \"%s\" for writing, skipping\n", fileoutname.c_str());
		return;
	}

	LogMsg("  generating list file \"%s\"\n", fileoutname.c_str());
	
    if(m_clsGlobalVars.m_bDumpInputGeo == false) {
        fprintf(foutlst, "* Refined FastCap file derived from %s\n", fileinname.c_str());
        fprintf(foutlst, "*\n");
        fprintf(foutlst, "* Total number of panels %lu\n", m_ulPanelNum[m_ucInteractionLevel]);
    }
    else {
        fprintf(foutlst, "* Dump FasterCap file derived from %s\n", fileinname.c_str());
        fprintf(foutlst, "*\n");
        fprintf(foutlst, "* Total number of panels %lu\n", m_ulInputPanelNum);
    }
    fprintf(foutlst, "*\n");
    
	// output the patches for every conductor

	// scan all conductor groups
	for(itc = m_stlConductors.begin(); itc != m_stlConductors.end(); itc++) {

		// build conductor or dielectric panel list file name
		condfilename = basefilename;
		condfilename += "_";
		condfilename += (*itc)->m_sName;
		condfilename += "_";
		condfilename += suffix;

		pos = condfilename.rfind("\\");
		// if there was a path, strip it
		if(pos != condfilename.npos ) {
			condfilenopath =  condfilename.substr(pos+1, condfilename.size()-pos-1);
		}
		else {
			condfilenopath = condfilename;
		}

		// if conductor
		if( (*itc)->m_bIsDiel == false ) {
			// must split the conductor in 'm_ucMaxSurfOutperm' surfaces with different output permittivity
			for(i=0; i<(*itc)->m_ucMaxSurfOutperm; i++) {
				// compose the extension
				ext = "_";
				dielIndexStr = wxString::Format("%d", i);
				ext += dielIndexStr.c_str();
				ext += ".txt";
				// and build the file names
				condfilenameext = condfilename + ext;
				condfilenopathext = condfilenopath + ext;

				if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
					// if last surface, must not add '+' at the end of the 'C' statement
					if(i == (*itc)->m_ucMaxSurfOutperm-1) {
						fprintf(foutlst, "C %s  %g  0.0 0.0 0.0\n", condfilenopathext.c_str(), (*itc)->m_dSurfOutperm[i][0]);
					}
					else {
						fprintf(foutlst, "C %s  %g  0.0 0.0 0.0 +\n", condfilenopathext.c_str(), (*itc)->m_dSurfOutperm[i][0]);
					}
				}
				else {
					// if last surface, must not add '+' at the end of the 'C' statement
					if(i == (*itc)->m_ucMaxSurfOutperm-1) {
						fprintf(foutlst, "C %s  %g-j%g  0.0 0.0 0.0\n", condfilenopathext.c_str(), (*itc)->m_dSurfOutperm[i][0], -1.0 * (*itc)->m_dSurfOutperm[i][1]);
					}
					else {
						fprintf(foutlst, "C %s  %g-j%g  0.0 0.0 0.0 +\n", condfilenopathext.c_str(), (*itc)->m_dSurfOutperm[i][0], -1.0 * (*itc)->m_dSurfOutperm[i][1]);
					}
				}
				
				// output the panels in the file 'condfilenameext'
				OutputPanelFile(condfilenameext, fileinname, fileoutname, itc, i);
			}
		}
		// if dielectric
		else {
			// nothing to add in the file name; close it with the extension
			condfilename += ".txt";
			condfilenopath += ".txt";

			if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
				if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
					fprintf(foutlst, "D %s  %g %g  0.0 0.0 0.0  %g %g %g\n", condfilenopath.c_str(),
					        (*itc)->m_dOutperm[0], (*itc)->m_dInperm[0], (*itc)->m_clsDielRef3DPoint.x,
					        (*itc)->m_clsDielRef3DPoint.y, (*itc)->m_clsDielRef3DPoint.z);
				}
				else {
					fprintf(foutlst, "D %s  %g-j%g %g-j%g  0.0 0.0 0.0  %g %g %g\n", condfilenopath.c_str(),
					        (*itc)->m_dOutperm[0], -1.0 * (*itc)->m_dOutperm[1], (*itc)->m_dInperm[0], -1.0 * (*itc)->m_dInperm[1], (*itc)->m_clsDielRef3DPoint.x,
					        (*itc)->m_clsDielRef3DPoint.y, (*itc)->m_clsDielRef3DPoint.z);
				}
			}
			else {
				if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
					fprintf(foutlst, "D %s  %g %g  0.0 0.0 0.0  %g %g %g\n", condfilenopath.c_str(),
					        (*itc)->m_dOutperm[0], (*itc)->m_dInperm[0],
					        (*itc)->m_clsDielRef3DPoint.x, (*itc)->m_clsDielRef3DPoint.y, AUTOREFINE_2D_ZCOORD / 2.0);
				}
				else {
					fprintf(foutlst, "D %s  %g-j%g %g-j%g  0.0 0.0 0.0  %g %g %g\n", condfilenopath.c_str(),
					        (*itc)->m_dOutperm[0], -1.0 * (*itc)->m_dOutperm[1], (*itc)->m_dInperm[0], -1.0 * (*itc)->m_dInperm[1],
					        (*itc)->m_clsDielRef3DPoint.x, (*itc)->m_clsDielRef3DPoint.y, AUTOREFINE_2D_ZCOORD / 2.0);
				}
			}

            // output the panels in the file 'condfilenameext'
            OutputPanelFile(condfilename, fileinname, fileoutname, itc);
		}
	}

	fclose(foutlst);
}

void CAutoRefine::OutputPanelFile(std::string condfilenameext, std::string fileinname, std::string fileoutname, StlAutoCondDeque::iterator itc, int dielIndex)
{
	FILE *foutgeo;
	
    foutgeo = fopen(condfilenameext.c_str(), "w");

    if(foutgeo == NULL) {
        ErrMsg("  Warning: cannot open \"%s\" for writing, skipping\n", condfilenameext.c_str());
    }
    else {
        if(m_clsGlobalVars.m_bDumpInputGeo == false) {
            fprintf(foutgeo, "0 Refined FastCap file derived from %s, referenced in list file %s\n", fileinname.c_str(), fileoutname.c_str());
        }
        else {
            fprintf(foutgeo, "0 Dump FasterCap file derived from %s, referenced in list file %s\n", fileinname.c_str(), fileoutname.c_str());
        }
        fprintf(foutgeo, "*\n");

        if(m_clsGlobalVars.m_bDumpInputGeo == false) {
            // output all panels in tree
            if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
                OutputPanelTree((*itc)->m_sName, (*itc)->m_uTopElement.m_pTopPanel, foutgeo, dielIndex);
            }
            else {
                OutputPanelTree((*itc)->m_sName, (*itc)->m_uTopElement.m_pTopSegment, foutgeo, dielIndex);
            }
        }
        else {
            // output all panels in the list
            if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
                OutputPanelList(itc, foutgeo, dielIndex);
            }
            else {
                ErrMsg("  Warning: input file dump not implemented yet for 2D case, skipping\n");
            }
        }

        fclose(foutgeo);
    }
}

void CAutoRefine::OutputPanelTree(char *condname, CAutoPanel *panel, FILE *fout, int dielIndex)
{

	// visit the tree

// debug
	//int i = 0;
	//if(i == 1)
	//	DebugOutputInteractions(panel);

	// termination condition
	if (panel->IsLeaf() == true) {
        OutputPanel(condname, panel, fout, dielIndex);
	}
	else {
		// otherwise, scan the tree
		OutputPanelTree(condname, (CAutoPanel*)panel->m_pLeft, fout, dielIndex);
		OutputPanelTree(condname, (CAutoPanel*)panel->m_pRight, fout, dielIndex);
	}
}

void CAutoRefine::OutputPanelTree(char *condname, CAutoSegment *panel, FILE *fout, int dielIndex)
{
	C2DVector dielrefpoint;

	// visit the tree

// debug
	//int i = 0;
	//if(i == 1)
	//	DebugOutputInteractions(panel);

	// termination condition
	if (panel->IsLeaf() == true) {
		// output panel only if there is no reference to the permittivity array (i.e. dielectric)
		// or if the dielectric index of the panel matches the current diel index
		if(dielIndex == AUTOREFINE_NO_DIEL_INDEX || panel->m_ucDielIndex == dielIndex) {
			// in this case, output panel
			if(m_clsGlobalVars.m_bOutputCharge == true && m_pLocalCondCharge != NULL) {
				fprintf(fout, "Q %s  %g %g %g %g %g %g %g %g %g %g %g %g", condname,
				        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, 0.0,
				        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, 0.0,
				        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, AUTOREFINE_2D_ZCOORD,
				        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, AUTOREFINE_2D_ZCOORD);
				if((panel->m_ucType & AUTOPANEL_OUTPERM_ELEMENT_LEVEL) == AUTOPANEL_OUTPERM_ELEMENT_LEVEL) {
					// local dielectric reference point
					dielrefpoint = panel->GetCentroid() + panel->GetDielNormal();
					fprintf(fout, "  %g %g %g", dielrefpoint.x, dielrefpoint.y, AUTOREFINE_2D_ZCOORD / 2.0);
				}
				fprintf(fout, "  %e\n", (*m_pLocalCondCharge)[panel->m_lIndex[AUTOREFINE_HIER_PRE_0_LEVEL]] / panel->GetDimension());
			}
			else if(m_clsGlobalVars.m_bKeepCharge == true) {
				fprintf(fout, "Q %s  %g %g %g %g %g %g %g %g %g %g %g %g", condname,
				        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, 0.0,
				        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, 0.0,
				        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, AUTOREFINE_2D_ZCOORD,
				        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, AUTOREFINE_2D_ZCOORD);
				if((panel->m_ucType & AUTOPANEL_OUTPERM_ELEMENT_LEVEL) == AUTOPANEL_OUTPERM_ELEMENT_LEVEL) {
					// local dielectric reference point
					dielrefpoint = panel->GetCentroid() + panel->GetDielNormal();
					fprintf(fout, "  %g %g %g", dielrefpoint.x, dielrefpoint.y, AUTOREFINE_2D_ZCOORD / 2.0);
				}
				fprintf(fout, "  %e\n", m_fGlobalCharges[panel->m_lIndex[AUTOREFINE_HIER_PRE_0_LEVEL]] / panel->GetDimension());
			}
			else {
				fprintf(fout, "Q %s  %g %g %g %g %g %g %g %g %g %g %g %g", condname,
				        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, 0.0,
				        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, 0.0,
				        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, AUTOREFINE_2D_ZCOORD,
				        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, AUTOREFINE_2D_ZCOORD);
				if((panel->m_ucType & AUTOPANEL_OUTPERM_ELEMENT_LEVEL) == AUTOPANEL_OUTPERM_ELEMENT_LEVEL) {
					// local dielectric reference point
					dielrefpoint = panel->GetCentroid() + panel->GetDielNormal();
					fprintf(fout, "  %g %g %g\n", dielrefpoint.x, dielrefpoint.y, AUTOREFINE_2D_ZCOORD / 2.0);
				}
				else {
					fprintf(fout, "\n");
				}
			}
		}
	}
	else {
		// otherwise, scan the tree
		OutputPanelTree(condname, (CAutoSegment*)panel->m_pLeft, fout, dielIndex);
		OutputPanelTree(condname, (CAutoSegment*)panel->m_pRight, fout, dielIndex);
	}
}

void CAutoRefine::OutputPanelList(StlAutoCondDeque::iterator itc, FILE *fout, int dielIndex)
{
	StlAutoPanelDeque::iterator itp;

	// scan every panel inside the conductor
    for(itp=(*itc)->m_stlPanels.begin(); itp!=(*itc)->m_stlPanels.end(); itp++) {
        OutputPanel((*itc)->m_sName, *itp, fout, dielIndex);
    }	
}

void CAutoRefine::OutputPanel(char *condname, CAutoPanel *panel, FILE *fout, int dielIndex)
{
   	C3DVector dielrefpoint;
   	CAutoQPanel *qpanel;

    // output panel only if there is no reference to the permittivity array (i.e. dielectric)
    // or if the dielectric index of the panel matches the current diel index
    if(dielIndex == AUTOREFINE_NO_DIEL_INDEX || panel->m_ucDielIndex == dielIndex) {
        // in this case, output panel
        //
        // if triangular panel
        if(panel->GetClass() == AUTOELEMENT_PANEL) {
            fprintf(fout, "T %s  %g %g %g  %g %g %g  %g %g %g", condname,
                        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
                        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z,
                        panel->m_clsVertex[2].x, panel->m_clsVertex[2].y, panel->m_clsVertex[2].z);
        }
        // if quadrilateral panel
        else {
            if (panel->GetClass() == AUTOELEMENT_QPANEL) {
                // re-cast to correct class (quick way)
                qpanel = (CAutoQPanel*)panel;
				fprintf(fout, "Q %s  %g %g %g  %g %g %g  %g %g %g  %g %g %g", condname,
                        qpanel->m_clsQVertex[0].x, qpanel->m_clsQVertex[0].y, qpanel->m_clsQVertex[0].z,
                        qpanel->m_clsQVertex[1].x, qpanel->m_clsQVertex[1].y, qpanel->m_clsQVertex[1].z,
                        qpanel->m_clsQVertex[2].x, qpanel->m_clsQVertex[2].y, qpanel->m_clsQVertex[2].z,
                        qpanel->m_clsQVertex[3].x, qpanel->m_clsQVertex[3].y, qpanel->m_clsQVertex[3].z);
            }
            else {
                // unknown class
                ErrMsg("Internal error: OutputPanel() trying to output a panel of an unknown class\n");
            }
        }

        if((panel->m_ucType & AUTOPANEL_OUTPERM_ELEMENT_LEVEL) == AUTOPANEL_OUTPERM_ELEMENT_LEVEL) {
            // local dielectric reference point
            dielrefpoint = panel->GetCentroid() + panel->GetDielNormal();
            fprintf(fout, "  %g %g %g", dielrefpoint.x, dielrefpoint.y, dielrefpoint.z);
        }

        if(m_clsGlobalVars.m_bOutputCharge == true && m_pLocalCondCharge != NULL && m_clsGlobalVars.m_bDumpInputGeo == false) {
            fprintf(fout, "  %e\n", (*m_pLocalCondCharge)[panel->m_lIndex[AUTOREFINE_HIER_PRE_0_LEVEL]] / panel->GetDimension());
        }
        else if(m_clsGlobalVars.m_bKeepCharge == true && m_clsGlobalVars.m_bDumpInputGeo == false) {
            fprintf(fout, "  %e\n", m_fGlobalCharges[panel->m_lIndex[AUTOREFINE_HIER_PRE_0_LEVEL]] / panel->GetDimension());
        }
        else {
            fprintf(fout, "\n");
        }
    }
}

#ifdef DEBUG_DUMP_BASIC
#ifdef DEBUG_DUMP_OTHER

void CAutoRefine::DebugOutputTree()
{
	FILE *fout;
	StlAutoCondDeque::iterator itc;

	fout = fopen("pot_tree.txt", "w");

	if(fout == NULL)
		return;

	// scan all conductor groups
	for(itc = m_stlConductors.begin(); itc != m_stlConductors.end(); itc++) {

		// output all panels in tree
		DebugOutputPanels((*itc)->m_uTopElement.m_pTopPanel, fout);
	}

	fclose(fout);
}

void CAutoRefine::DebugOutputPanels(CAutoPanel *panel, FILE *fout)
{
	// visit the tree

	InteractionC3DList::iterator iti;
	char depth1[256], depth2[256];

	for(iti = panel->m_stlInteractions[m_ucInteractionLevel].begin(); iti != panel->m_stlInteractions[m_ucInteractionLevel].end(); iti++) {

		DebugDepth(panel, depth1);
		DebugDepth((*iti).m_clsPanel, depth2);

		fprintf(fout, "C:%s:%s C:%s:%s   %e\n", panel->m_pCond->m_sName, depth1,
		        (*iti).m_clsPanel->m_pCond->m_sName, depth2, (*iti).m_dPotCoeff);

	}

	if (panel->IsLeaf() != true) {
		// otherwise, scan the tree
		DebugOutputPanels(panel->m_pLeft, fout);
		DebugOutputPanels(panel->m_pRight, fout);
	}
}

// remark: L and R information is printed in reverse order (first is innermost)
void CAutoRefine::DebugDepth(CAutoPanel *panel, char depth[])
{
	CAutoPanel *ppanel, *pppanel;
	int i;

	ppanel = panel;
	for(i=0; i<256 && ppanel->m_pParent != NULL; i++) {

		pppanel = ppanel;
		ppanel = ppanel->m_pParent;

		if(ppanel->m_pLeft == pppanel) {
			depth[i] = 'L';
		}
		else {
			depth[i] = 'R';
		}
	}
	depth[i] = '\0';
}

// warning: this function requires leaf panels to be indexed
void CAutoRefine::DebugOutputTree2()
{
	StlAutoCondDeque::iterator itc;
	long i, j, blockNo;
	FILE *fp;

	m_pBlockMtx = new long* [m_ulPanelNum[m_ucInteractionLevel]];
	m_pElemMtx = new double* [m_ulPanelNum[m_ucInteractionLevel]];
	for(i=0; i<m_ulPanelNum[m_ucInteractionLevel]; i++) {
		m_pBlockMtx[i] = new long [m_ulPanelNum[m_ucInteractionLevel]];
		m_pElemMtx[i] = new double [m_ulPanelNum[m_ucInteractionLevel]];
		for(j=0; j<m_ulPanelNum[m_ucInteractionLevel[m_ucInteractionLevel]]; j++) {
			m_pBlockMtx[i][j] = 0;
		}
	}

	fp = fopen("pot_tree2.txt", "w");

	ASSERT(fp != NULL);

	blockNo = 1;

	// scan all conductor groups
	for(itc = m_stlConductors.begin(); itc != m_stlConductors.end(); itc++) {

		// output all panels in tree
		DebugOutputPanels2((*itc)->m_uTopElement.m_pTopPanel, blockNo);
	}

	for(i=0; i<m_ulPanelNum[m_ucInteractionLevel]; i++) {
		for(j=0; j<m_ulPanelNum[m_ucInteractionLevel]; j++) {
			fprintf(fp, "%d:%e ", m_pBlockMtx[i][j], m_pElemMtx[i][j]);
		}
		fprintf(fp, "\n");
	}


	fclose(fp);

	for(i=0; i<m_ulPanelNum[m_ucInteractionLevel]; i++) {
		delete m_pBlockMtx[i];
		delete m_pElemMtx[i];
	}
	delete m_pBlockMtx;
	delete m_pElemMtx;
}

void CAutoRefine::DebugOutputPanels2(CAutoPanel *panel, long &blockNo)
{
	// visit the tree

	InteractionC3DList::iterator iti;

	for(iti = panel->m_stlInteractions[m_ucInteractionLevel].begin(); iti != panel->m_stlInteractions[m_ucInteractionLevel].end(); iti++) {

		DebugMarkBlock(panel, (*iti).m_clsPanel, blockNo, (*iti).m_dPotCoeff);
		blockNo++;
	}

	if (panel->IsLeaf() != true) {
		// scan the tree
		DebugOutputPanels2(panel->m_pLeft, blockNo);
		DebugOutputPanels2(panel->m_pRight, blockNo);
	}
}

void CAutoRefine::DebugMarkBlock(CAutoPanel *panel1, CAutoPanel *panel2, long blockNo, double potCoeff)
{

	if(panel1->IsLeaf() == true && panel2->IsLeaf() == true) {
		// if position already marked, there some serious error
		if(m_pBlockMtx[panel1->m_lIndex[m_ucInteractionLevel]][panel2->m_lIndex[m_ucInteractionLevel]] != 0) {
			ErrMsg("Debug block matrix dump error:\nBlock %d already marked with %d\n",
			       blockNo, m_pBlockMtx[panel1->m_lIndex[m_ucInteractionLevel]][panel2->m_lIndex[m_ucInteractionLevel]]);

			m_pBlockMtx[panel1->m_lIndex[m_ucInteractionLevel]][panel2->m_lIndex[m_ucInteractionLevel]] = 9.9999E9;
		}
		else {
			m_pBlockMtx[panel1->m_lIndex[m_ucInteractionLevel]][panel2->m_lIndex[m_ucInteractionLevel]] = blockNo;
			m_pElemMtx[panel1->m_lIndex[m_ucInteractionLevel]][panel2->m_lIndex[m_ucInteractionLevel]] = potCoeff;
		}
	}
	// first go to leaves of panel1
	else if(panel1->IsLeaf() != true) {
		DebugMarkBlock(panel1->m_pLeft, panel2, blockNo, potCoeff);
		DebugMarkBlock(panel1->m_pRight, panel2, blockNo, potCoeff);
	}
	// then only when on one of panel1's leaves, go to leaves of panel2
	else if(panel2->IsLeaf() != true) {
		DebugMarkBlock(panel1, panel2->m_pLeft, blockNo, potCoeff);
		DebugMarkBlock(panel1, panel2->m_pRight, blockNo, potCoeff);
	}
}

void CAutoRefine::DebugOutputHierarchy(CAutoPanel *panel)
{
	FILE *foutlst, *foutgeo;
	StlAutoCondDeque::iterator itc;
	char fileoutname[32], condfilename[32];
	C3DVector normal;

	sprintf(fileoutname, "_%s_hierarchy.lst", panel->m_pCond->m_sName);

	foutlst = fopen(fileoutname, "w");

	if(foutlst == NULL)
		return;


	fprintf(foutlst, "* Hierarchical FastCap file of a conductor\n");
	fprintf(foutlst, "*\n");


	// output bounding sphere
	//

	sprintf(condfilename, "_%s_lev%d_bsphere.qui", panel->m_pCond->m_sName, panel->m_iLevel);

	fprintf(foutlst, "C %s 1.0  %f %f %f\n", condfilename, panel->GetCentroid().x, panel->GetCentroid().y, panel->GetCentroid().z);

	foutgeo = fopen(condfilename, "w");

	if(foutgeo == NULL)
		return;

	gensphere(panel->GetMaxSideLen() / 2.0, 1, foutgeo);

	fclose(foutgeo);



	// output bounding box
	//

	sprintf(condfilename, "_%s_lev%d_bbox.qui", panel->m_pCond->m_sName, panel->m_iLevel);

	fprintf(foutlst, "C %s 1.0  %f %f %f\n", condfilename, 0.0, 0.0, 0.0);

	foutgeo = fopen(condfilename, "w");

	if(foutgeo == NULL)
		return;

	fprintf(foutgeo, "Q bbox  %f %f %f  %f %f %f  %f %f %f  %f %f %f\n",
	        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[0].x, panel->m_clsVertex[1].y, panel->m_clsVertex[0].z);
	fprintf(foutgeo, "Q bbox  %f %f %f  %f %f %f  %f %f %f  %f %f %f\n",
	        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[1].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[0].y, panel->m_clsVertex[1].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z,
	        panel->m_clsVertex[0].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z);
	fprintf(foutgeo, "Q bbox  %f %f %f  %f %f %f  %f %f %f  %f %f %f\n",
	        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[0].x, panel->m_clsVertex[1].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[0].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z,
	        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[1].z);
	fprintf(foutgeo, "Q bbox  %f %f %f  %f %f %f  %f %f %f  %f %f %f\n",
	        panel->m_clsVertex[1].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[0].y, panel->m_clsVertex[1].z);
	fprintf(foutgeo, "Q bbox  %f %f %f  %f %f %f  %f %f %f  %f %f %f\n",
	        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[0].y, panel->m_clsVertex[1].z,
	        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[1].z);
	fprintf(foutgeo, "Q bbox  %f %f %f  %f %f %f  %f %f %f  %f %f %f\n",
	        panel->m_clsVertex[0].x, panel->m_clsVertex[1].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[0].z,
	        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z,
	        panel->m_clsVertex[0].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z);

	fclose(foutgeo);


	// output normal
	//

	sprintf(condfilename, "_%s_lev%d_normal.qui", panel->m_pCond->m_sName, panel->m_iLevel);

	fprintf(foutlst, "C %s 1.0  0.0 0.0 0.0\n", condfilename);

	foutgeo = fopen(condfilename, "w");

	if(foutgeo == NULL)
		return;

	fprintf(foutgeo, "0 Normal vector for conductor %s, level %d\n", panel->m_pCond->m_sName, panel->m_iLevel);
	fprintf(foutgeo, "*\n");

	if(Mod(panel->m_clsNormal) != 0)
		normal = panel->m_clsNormal / Mod(panel->m_clsNormal);
	else
		normal = C3DVector(0,0,0);

	// output triangle pointing as normal
	fprintf(foutgeo, "T normal  0 0 0  0.1 0.1 0.1  %f %f %f\n", normal.x * 10, normal.y * 10, normal.z * 10);

	fclose(foutgeo);


	// output the node sub-patches
	//

	// build conductor or dielectric panel list file name
	sprintf(condfilename, "_%s_lev%d.qui", panel->m_pCond->m_sName, panel->m_iLevel);

	// if conductor
	if( panel->m_pCond->m_bIsDiel == false ) {
		fprintf(foutlst, "C %s  %f  0.0 0.0 0.0\n", condfilename, panel->m_pCond->m_dOutperm[0]);
	}
	// if dielectric
	else {
		fprintf(foutlst, "D %s  %f %f  0.0 0.0 0.0  %f %f %f\n", condfilename,
		        panel->m_pCond->m_dOutperm[0], panel->m_pCond->m_dInperm[0], panel->m_pCond->m_clsDielRef3DPoint.x,
		        panel->m_pCond->m_clsDielRef3DPoint.y, panel->m_pCond->m_clsDielRef3DPoint.z);
	}

	foutgeo = fopen(condfilename, "w");

	fprintf(foutgeo, "0 Debug FastCap file, referenced in list file %s\n", fileoutname);
	fprintf(foutgeo, "*\n");

	// output all panels in tree
	OutputPanelTree(panel->m_pCond->m_sName, panel, foutgeo);

	fclose(foutgeo);

	fclose(foutlst);
}

void CAutoRefine::DebugOutputInteractions(CAutoPanel *panel)
{
	FILE *foutlst, *foutgeo;
	StlAutoCondDeque::iterator itc;
	char fileoutname[32], condfilename[32];
	C3DVector normal;
	InteractionC3DList::iterator iti;
	double scale = 1.1;

	sprintf(fileoutname, "_%s_lev%d.lst", panel->m_pCond->m_sName, panel->m_iLevel);

	foutlst = fopen(fileoutname, "w");

	if(foutlst == NULL)
		return;


	fprintf(foutlst, "* Interactions FastCap file\n");
	fprintf(foutlst, "*\n");


	// output all structures
	//

	// build conductor or dielectric panel list file name
	sprintf(condfilename, "_%s_lev%d_all.qui", panel->m_pCond->m_sName, panel->m_iLevel);

	fprintf(foutlst, "C %s  %f  0.0 0.0 0.0\n", condfilename, panel->m_pCond->m_dOutperm[0]);

	foutgeo = fopen(condfilename, "w");

	fprintf(foutgeo, "0 Debug FastCap file, referenced in list file %s\n", fileoutname);
	fprintf(foutgeo, "*\n");

	// scan all conductor groups
	for(itc = m_stlConductors.begin(); itc != m_stlConductors.end(); itc++) {
		// output all panels in tree
		DebugOutputPanelTree((*itc)->m_sName, (*itc)->m_uTopElement.m_pTopPanel, foutgeo);
	}
	fclose(foutgeo);


	// output the panel interactions
	//

	// build conductor or dielectric panel list file name
	sprintf(condfilename, "_%s_lev%d_interact.qui", panel->m_pCond->m_sName, panel->m_iLevel);

	fprintf(foutlst, "C %s  %f  0.0 0.0 0.0\n", condfilename, panel->m_pCond->m_dOutperm[0]);

	foutgeo = fopen(condfilename, "w");

	fprintf(foutgeo, "0 Debug FastCap file, referenced in list file %s\n", fileoutname);
	fprintf(foutgeo, "*\n");

	// this panel
	fprintf(foutgeo, "T %s %f %f %f %f %f %f %f %f %f\n", "panels",
	        panel->m_clsVertex[0].x * scale, panel->m_clsVertex[0].y * scale, panel->m_clsVertex[0].z * scale,
	        panel->m_clsVertex[1].x * scale, panel->m_clsVertex[1].y * scale, panel->m_clsVertex[1].z * scale,
	        panel->m_clsVertex[2].x * scale, panel->m_clsVertex[2].y * scale, panel->m_clsVertex[2].z * scale);

	for(iti = panel->m_stlInteractions[m_ucInteractionLevel].begin(); iti != panel->m_stlInteractions[m_ucInteractionLevel].end(); iti++) {

		fprintf(foutgeo, "T %s %f %f %f %f %f %f %f %f %f\n", "panels",
		        (*iti).m_clsPanel->m_clsVertex[0].x * scale, (*iti).m_clsPanel->m_clsVertex[0].y * scale, (*iti).m_clsPanel->m_clsVertex[0].z * scale,
		        (*iti).m_clsPanel->m_clsVertex[1].x * scale, (*iti).m_clsPanel->m_clsVertex[1].y * scale, (*iti).m_clsPanel->m_clsVertex[1].z * scale,
		        (*iti).m_clsPanel->m_clsVertex[2].x * scale, (*iti).m_clsPanel->m_clsVertex[2].y * scale, (*iti).m_clsPanel->m_clsVertex[2].z * scale);
	}

	fclose(foutgeo);

	fclose(foutlst);
}

// output interacting panels / superpanels in terms of their leaves
void CAutoRefine::DebugOutputInteracting(CAutoPanel *panel1, CAutoPanel *panel2)
{
	FILE *foutlst, *foutgeo;
	StlAutoCondDeque::iterator itc;
	char fileoutname[32], condfilename[32];
	C3DVector normal;
	InteractionC3DList::iterator iti;
	double scale = 1.1;

	sprintf(fileoutname, "_%s_interacting.lst", panel1->m_pCond->m_sName);

	foutlst = fopen(fileoutname, "w");

	if(foutlst == NULL)
		return;


	fprintf(foutlst, "* Interacting FastCap file\n");
	fprintf(foutlst, "*\n");


	// output first panel structure
	//

	// build conductor or dielectric panel list file name
	sprintf(condfilename, "_%s_%x.qui", panel1->m_pCond->m_sName, panel1);

	fprintf(foutlst, "C %s  %f  0.0 0.0 0.0\n", condfilename, panel1->m_pCond->m_dOutperm[0]);

	foutgeo = fopen(condfilename, "w");

	fprintf(foutgeo, "0 Debug FastCap file, referenced in list file %s\n", fileoutname);
	fprintf(foutgeo, "*\n");

	// output all panels in tree
	DebugOutputPanelTree("one", panel1, foutgeo);

	fclose(foutgeo);

	// output second panel structure
	//

	// build conductor or dielectric panel list file name
	sprintf(condfilename, "_%s_%x.qui", panel2->m_pCond->m_sName, panel2);

	fprintf(foutlst, "C %s  %f  0.0 0.0 0.0\n", condfilename, panel2->m_pCond->m_dOutperm[0]);

	foutgeo = fopen(condfilename, "w");

	fprintf(foutgeo, "0 Debug FastCap file, referenced in list file %s\n", fileoutname);
	fprintf(foutgeo, "*\n");

	// output all panels in tree
	DebugOutputPanelTree("two", panel2, foutgeo);

	fclose(foutgeo);


	fclose(foutlst);
}

void CAutoRefine::DebugOutputPanelTree(char *condname, CAutoPanel *panel, FILE *fout)
{
	// visit the tree

	// termination condition
	if (panel->IsLeaf() == true) {
		// in this case, output panel
		fprintf(fout, "T %s %f %f %f %f %f %f %f %f %f\n", condname,
		        panel->m_clsVertex[0].x, panel->m_clsVertex[0].y, panel->m_clsVertex[0].z,
		        panel->m_clsVertex[1].x, panel->m_clsVertex[1].y, panel->m_clsVertex[1].z,
		        panel->m_clsVertex[2].x, panel->m_clsVertex[2].y, panel->m_clsVertex[2].z);
	}
	else {
		// otherwise, scan the tree
		DebugOutputPanelTree(condname, panel->m_pLeft, fout);
		DebugOutputPanelTree(condname, panel->m_pRight, fout);
	}
}

#endif // #ifdef DEBUG_DUMP_OTHER
#endif  // #ifdef DEBUG_DUMP_BASIC

void CAutoRefine::DebugDumpInteractions()
{
	unsigned long linkIndex, chunk, posInChunk, block;
	FILE *fp;

	fp = fopen("links.txt", "w");

	ASSERT(fp != NULL);

	m_ulCurrBlock = 0;
	// only if we went out-of-core, pre-load first set of chunks
	if(m_ulBlocksNum > 1) {
		LoadLinks(m_ulCurrBlock);
	}

	for(linkIndex=0; linkIndex<GetLinksNum(AUTOREFINE_HIER_PRE_0_LEVEL); linkIndex++) {

		chunk = linkIndex / AUTOREFINE_LINK_CHUNK_SIZE;
		posInChunk = linkIndex % AUTOREFINE_LINK_CHUNK_SIZE;
		// determine in which block the current chunk is ('m_ulLinkChunkNum' is the number of chunks per block)
		block = chunk / m_ulLinkChunkNum[m_ucInteractionLevel];
		// if not in current block
		if(block != m_ulCurrBlock) {
			LoadLinks(block);
		}
		// adjust chunk to position within the current block
		chunk -= m_ulLinkChunkNum[m_ucInteractionLevel] * m_ulCurrBlock;

		// actually dump to file in scilab matrix format
		fprintf(fp, "%f %p\n", m_dPotCoeffLinks[m_ucInteractionLevel][chunk][posInChunk], m_pdPanelPtrLinks[m_ucInteractionLevel][chunk][posInChunk]);
	}

	fclose(fp);
}

void CAutoRefine::DeletePanelsAndConductors()
{
	StlAutoCondDeque::iterator itc;

	// scan all conductor groups
	for(itc = m_stlConductors.begin(); itc != m_stlConductors.end(); itc++) {
		// delete panels tree
		if((*itc)->m_uTopElement.m_pTopPanel != NULL) {
			((*itc)->m_uTopElement.m_pTopPanel)->DeleteElementsTree();
		}
        // empty panel list
        (*itc)->m_stlPanels.clear();
		// delete conductor
		delete *itc;
	}
    // empty conductor list
	m_stlConductors.clear();
}


/*
void CAutoRefine::DeleteRefinedPanels()
{
	StlAutoCondDeque::iterator itc;

	// scan all conductor groups
	for(itc = m_stlConductors.begin(); itc != m_stlConductors.end(); itc++) {
		// delete panels tree
		DeleteRefinedPanelTree((*itc)->m_uTopElement.m_pTopPanel, AUTOPANEL_IS_SUPER_NODE);
	}
}

void CAutoRefine::DeleteRefinedPanelTree(CAutoPanel *panel, unsigned char prevpanelType)
{
	// visit the tree

	if (panel->IsLeaf() != true) {
		// scan the subtree, only after delete panel
		DeleteRefinedPanelTree(panel->m_pLeft, panel->m_ucType);
		DeleteRefinedPanelTree(panel->m_pRight, panel->m_ucType);
	}

	// then delete panel, if it is a refined one (i.e. previous panel is not a super node)
	if( (prevpanelType & AUTOPANEL_IS_SUPER_NODE) == AUTOPANEL_IS_SUPER_NODE) {
		delete panel;
	}
}
*/


// Recursive refinement
//
// Remark: before calling, use SetCurrentConductor() to set the pointer
// to the current conductor being processed (the one to which the 'panel' belongs)
// for correct handling of SelfPotential()
int CAutoRefine::RefineSelf(CAutoPanel *panel)
{
	double potestRe, potestIm;
//	unsigned long chunk, posInChunk, block;

	if(g_bFCContinue == false) {
		return FC_USER_BREAK;
	}

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	if(panel->IsLeaf() == true) {

		ASSERT((panel->m_ucType & AUTOPANEL_IS_SUPER_NODE) != AUTOPANEL_IS_SUPER_NODE);

		// store interaction for auto potential only if at bottom level, but
		// will use it for all levels (optimization to compute self-potential only once)
		if(m_ucInteractionLevel == AUTOREFINE_HIER_PRE_0_LEVEL) {

			// if just counting
			if(m_bComputeLinks == false) {
				// do nothing; number of links does not include self-potentials, since
				// self potentials are stored in 'm_clsSelfPotCoeff' and 'm_clsImgSelfPotCoeff'

#ifdef DEBUG_DUMP_BASIC
				m_iaLinksBtwLevels[panel->m_iLevel][panel->m_iLevel]++;
#endif
			}
			else {
				// calculate potential and store interaction

				// if not already calculated, compute the self potential and store interaction
				// (happens when RefineSelf() is called multiple times due to OOC division in chunks)
				if(m_clsSelfPotCoeff[panel->m_lIndex[m_ucInteractionLevel]] == 0.0) {

					// actually calculate potential estimate
					SelfPotential(panel, &potestRe, &potestIm);
					_ASSERT(!isnan(potestRe));
					_ASSERT(isfinite(potestRe));
					if(isnan(potestRe) || !isfinite(potestRe)) {
						if(m_clsGlobalVars.m_bWarnGivenNaN == false) {
							// signal we already warned the user
							m_clsGlobalVars.m_bWarnGivenNaN = true;
							ErrMsg("Error: self-potential calculation failed.\n");
							ErrMsg("       Remark: the precision of the result is affected.\n");
						}
						if(m_clsGlobalVars.m_bVerboseOutput == true) {
							if(!isfinite(potestRe)) {
								ErrMsg("Error: self-potential calculation gave infinite value\n");
							}
							else {
								ErrMsg("Error: self-potential calculation gave 'not a number' value (NaN)\n");
							}
							panel->ErrorPrintCoords();
						}
					}

					m_clsSelfPotCoeff[panel->m_lIndex[m_ucInteractionLevel]] = potestRe;

					// if dielectric
					if( (panel->m_ucType & AUTOPANEL_IS_DIEL) != 0 ) {
						// and if complex permittivity
						if( m_clsGlobalVars.m_ucHasCmplxPerm != AUTOREFINE_REAL_PERM ) {
							_ASSERT(!isnan(potestIm));
							_ASSERT(isfinite(potestIm));
							if(isnan(potestIm) || !isfinite(potestIm)) {
								if(m_clsGlobalVars.m_bWarnGivenNaN == false) {
									// signal we already warned the user
									m_clsGlobalVars.m_bWarnGivenNaN = true;
									ErrMsg("Error: electric field discontinuity calculation failed on a dielectric panel.\n");
									ErrMsg("       Remark: the precision of the result is affected.\n");
								}
								if(m_clsGlobalVars.m_bVerboseOutput == true) {
									if(!isfinite(potestRe)) {
										ErrMsg("Error: electric field discontinuity calculation gave imaginary part infinite value\n");
									}
									else {
										ErrMsg("Error: electric field discontinuity calculation gave imaginary part 'not a number' value (NaN)\n");
									}
									panel->ErrorPrintCoords();
								}
							}

							m_clsImgSelfPotCoeff[panel->m_lIndex[m_ucInteractionLevel]] = potestIm;
						}
					}

					// also store the specific permittivity index in the array (useful for the conductors,
					// since the charges on conductor panels will be multiplied by the permittivity of the
					// surrounding dielectric at the end of the solve pass)
					m_pucDielIndex[panel->m_lIndex[m_ucInteractionLevel]] = panel->m_ucDielIndex;
				}

			}
		}

	}
	else {
		// refine children one against the other
		RefineMutual((CAutoPanel*)(panel->m_pLeft), (CAutoPanel*)(panel->m_pRight));
		// then call recursively RefineSelf() for each child
		RefineSelf((CAutoPanel*)(panel->m_pLeft));
		RefineSelf((CAutoPanel*)(panel->m_pRight));
	}

	// return to upper level
	m_iLevel--;

	return FC_NORMAL_END;
}


// Recursive refinement
//
// Remark: before calling, use SetCurrentConductor() to set the pointer
// to the current conductor being processed (the one to which the 'panel' belongs)
// for correct handling of SelfPotential()
int CAutoRefine::RefineSelf(CAutoSegment *panel)
{
	double potestRe, potestIm;
//	unsigned long chunk, posInChunk, block;

	if(g_bFCContinue == false) {
		return FC_USER_BREAK;
	}

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	if(panel->IsLeaf() == true) {

		ASSERT((panel->m_ucType & AUTOPANEL_IS_SUPER_NODE) != AUTOPANEL_IS_SUPER_NODE);

		// store interaction for auto potential only if at bottom level, but
		// will use it for all levels (optimization to compute self-potential only once)
		if(m_ucInteractionLevel == AUTOREFINE_HIER_PRE_0_LEVEL) {

			// if just counting
			if(m_bComputeLinks == false) {
				// do nothing; number of links does not include self-potentials, since
				// self potentials are stored in 'm_clsSelfPotCoeff' and 'm_clsImgSelfPotCoeff'

#ifdef DEBUG_DUMP_BASIC
				m_iaLinksBtwLevels[panel->m_iLevel][panel->m_iLevel]++;
#endif
			}
			else {
				// calculate potential and store interaction

				// if not already calculated, compute the self potential and store interaction
				// (happens when RefineSelf() is called multiple times due to OOC division in chunks)
				if(m_clsSelfPotCoeff[panel->m_lIndex[m_ucInteractionLevel]] == 0.0) {

					// actually calculate potential estimate
					SelfPotential(panel, &potestRe, &potestIm);
					_ASSERT(!isnan(potestRe));
					_ASSERT(isfinite(potestRe));
					if(isnan(potestRe) || !isfinite(potestRe)) {
						if(m_clsGlobalVars.m_bWarnGivenNaN == false) {
							// signal we already warned the user
							m_clsGlobalVars.m_bWarnGivenNaN = true;
							ErrMsg("Error: self-potential calculation failed.\n");
							ErrMsg("       Remark: the precision of the result is affected.\n");
						}
						if(m_clsGlobalVars.m_bVerboseOutput == true) {
							if(!isfinite(potestRe)) {
								ErrMsg("Error: self-potential calculation gave infinite value\n");
							}
							else {
								ErrMsg("Error: self-potential calculation gave 'not a number' value (NaN)\n");
							}
							panel->ErrorPrintCoords();
						}
					}

					m_clsSelfPotCoeff[panel->m_lIndex[m_ucInteractionLevel]] = potestRe;

					// if dielectric
					if( (panel->m_ucType & AUTOPANEL_IS_DIEL) != 0 ) {
						// and if complex permittivity
						if( m_clsGlobalVars.m_ucHasCmplxPerm != AUTOREFINE_REAL_PERM ) {
							_ASSERT(!isnan(potestIm));
							_ASSERT(isfinite(potestIm));
							if(isnan(potestIm) || !isfinite(potestIm)) {
								if(m_clsGlobalVars.m_bWarnGivenNaN == false) {
									// signal we already warned the user
									m_clsGlobalVars.m_bWarnGivenNaN = true;
									ErrMsg("Error: electric field discontinuity calculation failed on a dielectric panel.\n");
									ErrMsg("       Remark: the precision of the result is affected.\n");
								}
								if(m_clsGlobalVars.m_bVerboseOutput == true) {
									if(!isfinite(potestRe)) {
										ErrMsg("Error: electric field discontinuity calculation gave imaginary part infinite value\n");
									}
									else {
										ErrMsg("Error: electric field discontinuity calculation gave imaginary part 'not a number' value (NaN)\n");
									}
									panel->ErrorPrintCoords();
								}
							}

							m_clsImgSelfPotCoeff[panel->m_lIndex[m_ucInteractionLevel]] = potestIm;
						}
					}

					// also store the specific permittivity index in the array (useful for the conductors,
					// since the charges on conductor panels will be multiplied by the permittivity of the
					// surrounding dielectric at the end of the solve pass)
					m_pucDielIndex[panel->m_lIndex[m_ucInteractionLevel]] = panel->m_ucDielIndex;
				}

			}
		}

	}
	else {
		// refine children one against the other
		RefineMutual((CAutoSegment*)(panel->m_pLeft), (CAutoSegment*)(panel->m_pRight));
		// then call recursively RefineSelf() for each child
		RefineSelf((CAutoSegment*)(panel->m_pLeft));
		RefineSelf((CAutoSegment*)(panel->m_pRight));
	}

	// return to upper level
	m_iLevel--;

	return FC_NORMAL_END;
}

// Auto potential of a panel
//
// Remark: before calling, use SetCurrentConductor() to set the pointer
// to the current conductor being processed (the one to which the 'panel' belongs)
// for correct handling of dielectric constants
void CAutoRefine::SelfPotential(CAutoPanel *panel, double *potestRe, double *potestIm)
{
	double eaRe, eaIm, ebRe, ebIm, lowerterm;
	double ebR_m_eaR, ebI_m_eaI, ebR_p_eaR, ebI_p_eaI;

	// zero the imaginary part for consistency (not used for conductor panels
	// or for real perm options
	*potestIm = 0.0f;

	// if panel is dielectric
	if( (panel->m_ucType & AUTOPANEL_IS_DIEL) != 0 ) {

		// determine epsilon on outperm side;
		// it is our convention that eb is on outperm side
//		eb = panel->m_pCond->m_dOutperm;
//		ea = panel->m_pCond->m_dInperm;
		ebRe = m_pCurrentConductor->m_dOutperm[0];
		eaRe = m_pCurrentConductor->m_dInperm[0];

		if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
			*potestRe = (eaRe + ebRe) / (TWO_TIMES_E0 * (ebRe - eaRe) * panel->GetDimension());
		}
		else {
			ebIm = m_pCurrentConductor->m_dOutperm[1];
			eaIm = m_pCurrentConductor->m_dInperm[1];

			ebR_m_eaR = ebRe - eaRe;
			ebR_p_eaR = ebRe + eaRe;
			ebI_m_eaI = ebIm - eaIm;
			ebI_p_eaI = ebIm + eaIm;

			lowerterm = TWO_TIMES_E0 * panel->GetDimension() * (ebR_m_eaR*ebR_m_eaR + ebI_m_eaI*ebI_m_eaI);

			*potestRe = (ebR_m_eaR * ebR_p_eaR + ebI_m_eaI * ebI_p_eaI) / lowerterm;
			*potestIm = (ebR_m_eaR * ebI_p_eaI - ebI_m_eaI * ebR_p_eaR) / lowerterm;
		}
	}
	// if panel is conductor
	else {
		// if not supernode
//        if( (panel->m_ucType & AUTOPANEL_IS_SUPER_NODE) != AUTOPANEL_IS_SUPER_NODE) {
		if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
			*potestRe = m_clsPotential.Auto(panel->m_clsVertex);
			*potestIm = 0.0f;
		}
		else {
			// collocation (using m_clsPotential.Auto() mixed with collocation for mutual potential gives instabilities)
			*potestRe = m_clsPotential.PotentialOpt(panel->GetCentroid(), panel->m_clsVertex) / FOUR_PI_TIMES_E0;
			if( *potestRe == 0.0) {
				if(m_clsGlobalVars.m_bWarnGivenSelfPot == false) {
					// signal we already warned the user
					m_clsGlobalVars.m_bWarnGivenSelfPot = true;
					ErrMsg("Warning: self-potential calculation equal to zero found during potential calculation\n");
					ErrMsg("         This may be caused by the presence of very small panels.\n");
					ErrMsg("         Remark: the precision of the result is affected.\n");
				}
				if(m_clsGlobalVars.m_bVerboseOutput == true) {

					ErrMsg("Warning: zero self-potential found for panel %lx.\n", panel);
					panel->ErrorPrintCoords();
				}
			}
		}
	}

	_ASSERT(!isnan(*potestRe));
	_ASSERT(isfinite(*potestRe));
	_ASSERT(!isnan(*potestIm));
	_ASSERT(isfinite(*potestIm));
}


// Auto potential of a panel
//
// Remark: before calling, use SetCurrentConductor() to set the pointer
// to the current conductor being processed (the one to which the 'panel' belongs)
// for correct handling of dielectric constants
void CAutoRefine::SelfPotential(CAutoSegment *panel, double *potestRe, double *potestIm)
{
	double eaRe, eaIm, ebRe, ebIm, lowerterm;
	double ebR_m_eaR, ebI_m_eaI, ebR_p_eaR, ebI_p_eaI;

	// zero the imaginary part for consistency (not used for conductor panels
	// or for real perm options
	*potestIm = 0.0f;

	// if panel is dielectric
	if( (panel->m_ucType & AUTOPANEL_IS_DIEL) != 0 ) {

		// determine epsilon on outperm side;
		// it is our convention that eb is on outperm side
//		eb = panel->m_pCond->m_dOutperm;
//		ea = panel->m_pCond->m_dInperm;
		ebRe = m_pCurrentConductor->m_dOutperm[0];
		eaRe = m_pCurrentConductor->m_dInperm[0];

		if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
			// must multiply by PI: formula in MKS should be (ea+eb)/((eb-ea)*2*E0)
			// so since everything is scaled by TWO_PI_TIMES_E0, here we multiply by PI
			*potestRe = PI * (eaRe + ebRe) / ((ebRe - eaRe) * panel->GetDimension());
		}
		else {
			ebIm = m_pCurrentConductor->m_dOutperm[1];
			eaIm = m_pCurrentConductor->m_dInperm[1];

			ebR_m_eaR = ebRe - eaRe;
			ebR_p_eaR = ebRe + eaRe;
			ebI_m_eaI = ebIm - eaIm;
			ebI_p_eaI = ebIm + eaIm;

			lowerterm = panel->GetDimension() * (ebR_m_eaR*ebR_m_eaR + ebI_m_eaI*ebI_m_eaI);

			*potestRe = PI * (ebR_m_eaR * ebR_p_eaR + ebI_m_eaI * ebI_p_eaI) / lowerterm;
			*potestIm = PI * (ebR_m_eaR * ebI_p_eaI - ebI_m_eaI * ebR_p_eaR) / lowerterm;
		}
	}
	// if panel is conductor
	else {

		/*		if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
					*potestRe = m_clsPotential.Auto(panel->m_clsVertex);
					*potestIm = 0.0f;
				}
				else {
		*/
		// collocation (using m_clsPotential.Auto() mixed with collocation for mutual potential gives instabilities)
		*potestRe = m_clsPotential.PotentialOpt(panel->GetCentroid(), panel->m_clsVertex);
		if( *potestRe == 0.0) {
			if(m_clsGlobalVars.m_bWarnGivenSelfPot == false) {
				// signal we already warned the user
				m_clsGlobalVars.m_bWarnGivenSelfPot = true;
				ErrMsg("Warning: self-potential calculation equal to zero found during potential calculation\n");
				ErrMsg("         This may be caused by the presence of very small panels.\n");
				ErrMsg("         Remark: the precision of the result is affected.\n");
			}
			if(m_clsGlobalVars.m_bVerboseOutput == true) {

				ErrMsg("Warning: zero self-potential found for panel %lx.\n", panel);
				panel->ErrorPrintCoords();
			}
		}
//		}
	}

	_ASSERT(!isnan(*potestRe));
	_ASSERT(isfinite(*potestRe));
	_ASSERT(!isnan(*potestIm));
	_ASSERT(isfinite(*potestIm));
}

/*
double CAutoRefine::RecurseSelfPotential(CAutoPanel *panel)
{
	double potestim;

	// can stop at the first 'real' panel, no need to go down to the leaves
	if((panel->m_ucType & AUTOPANEL_IS_SUPER_NODE) != AUTOPANEL_IS_SUPER_NODE ) {
//		if(m_pCentroid == NULL) {
//			m_pCentroid = &(panel->GetCentroid());
//		}
//		potestim = m_clsPotential.PotentialOpt(*m_pCentroid, panel->m_clsVertex, false);
		potestim = m_clsPotential.PotentialOpt(panel->GetCentroid(), panel->m_clsVertex, false);
	}
	else {
		potestim = RecurseSelfPotential(panel->m_pLeft);
		potestim += RecurseSelfPotential(panel->m_pRight);
	}

	return potestim;
}

double CAutoRefine::RecurseMutualPotential(CAutoPanel *panel)
{
	double potestim;

	if((panel->m_ucType & AUTOPANEL_IS_SUPER_NODE) != AUTOPANEL_IS_SUPER_NODE ) {
		potestim = m_clsPotential.PotentialOpt(*m_pCentroid, panel->m_clsVertex, false);
	}
	else {
		potestim = RecurseMutualPotential(panel->m_pLeft);
		potestim += RecurseMutualPotential(panel->m_pRight);
	}

	return potestim;
}
*/

// Recursive refinement
int CAutoRefine::Discretize(CAutoPanel *panel)
{
	double rmax;
	bool discretize;
	int ret;

	ret = FC_NORMAL_END;

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	// reset link counter
	panel->m_ulLinkIndexStart[m_ucInteractionLevel] = 0;
	panel->m_ulLinkIndexEnd[m_ucInteractionLevel] = 0;

	// get max lenght of panel side
	rmax = panel->GetMaxSideLen();

	if(panel->IsLeaf() == true) {

		// criteria for discretization
		discretize = false;
		// if side too large, candidate for discretization
		if(rmax > m_clsGlobalVars.m_dMaxDiscSide) {
			if(m_clsGlobalVars.m_bRefineCharge == true) {
				if(panel->m_dCharge >= m_dMidSigma) {
					discretize = true;
				}
			}
			else {
				discretize = true;
			}
		}

		// if decided to discretize, subdivide and recurse
		if(discretize == true) {
			ret = panel->Subdivide();

			if(ret != FC_NORMAL_END) {
				return ret;
			}

			if(m_clsGlobalVars.m_bRefineCharge == true) {
				// avoid recursive levels of discretization, that could be a problem
				// when an 'old' large panel falls within the discretization criteria
				// and therefore gets refined in one shot to the finest level of mesh
				panel->m_pLeft->m_dCharge = m_dMinSigma;
				panel->m_pRight->m_dCharge = m_dMinSigma;
			}

			ret = Discretize((CAutoPanel*)(panel->m_pLeft));

			if(ret != FC_NORMAL_END) {
				return ret;
			}

			ret = Discretize((CAutoPanel*)(panel->m_pRight));

			if(ret != FC_NORMAL_END) {
				return ret;
			}
		}
		else {
			// increment count of leaf panels
			m_ulPanelNum[m_ucInteractionLevel]++;
		}
	}
	else {
		ret = Discretize((CAutoPanel*)(panel->m_pLeft));

		if(ret != FC_NORMAL_END) {
			return ret;
		}

		ret = Discretize((CAutoPanel*)(panel->m_pRight));

		if(ret != FC_NORMAL_END) {
			return ret;
		}
	}

	// return to upper level
	m_iLevel--;

	return ret;
}

// Recursive refinement
int CAutoRefine::DiscretizeSelf(CAutoPanel *panel)
{
	int ret;

	ret = FC_NORMAL_END;

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	// reset link counter
	panel->m_ulLinkIndexStart[m_ucInteractionLevel] = 0;
	panel->m_ulLinkIndexEnd[m_ucInteractionLevel] = 0;

	// if leaf panel at top level, must start discretization
	if(panel->IsLeaf() == true && m_iLevel == 1) {
		ret = panel->Subdivide();

		if(ret != FC_NORMAL_END) {
			return ret;
		}
	}

	if(panel->IsLeaf() != true) {
		DiscretizeMutual((CAutoPanel*)(panel->m_pLeft), (CAutoPanel*)(panel->m_pRight), true);
		// then call recursively DiscretizeSelf() for each child
		DiscretizeSelf((CAutoPanel*)(panel->m_pLeft));
		DiscretizeSelf((CAutoPanel*)(panel->m_pRight));
	}
	else {
		// increment count of leaf panels
		m_ulPanelNum[m_ucInteractionLevel]++;
	}

	// increment count of nodes
	m_ulNodeNum[m_ucInteractionLevel]++;

	// return to upper level
	m_iLevel--;

	return ret;
}

// Recursive refinement
int CAutoRefine::DiscretizeSelf(CAutoSegment *panel)
{
	int ret;

	ret = FC_NORMAL_END;

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	// reset link counter
	panel->m_ulLinkIndexStart[m_ucInteractionLevel] = 0;
	panel->m_ulLinkIndexEnd[m_ucInteractionLevel] = 0;

	// if leaf panel at top level, must start discretization
	if(panel->IsLeaf() == true && m_iLevel == 1) {
		ret = panel->Subdivide();

		if(ret != FC_NORMAL_END) {
			return ret;
		}
	}

	if(panel->IsLeaf() != true) {
		DiscretizeMutual((CAutoSegment*)(panel->m_pLeft), (CAutoSegment*)(panel->m_pRight), true);
		// then call recursively DiscretizeSelf() for each child
		DiscretizeSelf((CAutoSegment*)(panel->m_pLeft));
		DiscretizeSelf((CAutoSegment*)(panel->m_pRight));
	}
	else {
		// increment count of leaf panels
		m_ulPanelNum[m_ucInteractionLevel]++;
	}

	// increment count of nodes
	m_ulNodeNum[m_ucInteractionLevel]++;

	// return to upper level
	m_iLevel--;

	return ret;
}

// Recursive refinement
int CAutoRefine::DiscretizeMutual(CAutoPanel *panel1, CAutoPanel *panel2, bool selfCond)
{
	bool forcerefinement;
	C3DVector dist;
	double rdist, r1, r2, rmax, ccoeff;
	double panel1crit, panel2crit;
	int ret;

	ret = FC_NORMAL_END;

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	dist = panel1->GetCentroid() - panel2->GetCentroid();

	rdist = Mod(dist);

	// get max lenght of panel side
	r1 = panel1->GetMaxSideLen();
	r2 = panel2->GetMaxSideLen();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	// check interaction condition against the threshold
	//

	// use different criteria for conductor and dielectric

	// if mutual interaction between panels belonging to the same conductor,
	// calculate the curvature coefficient, to give a different weight according
	// to the relative angle between the two panels. Parallel panels are treated
	// as if they were farther from each other w.r.t. opposite panels.
	// Remark: we are assuming the correct orientation of the normal of the panels,
	// all pointing outwards from a solid's interior
	if(selfCond == true) {
		ccoeff = (m_clsGlobalVars.m_dMeshCurvCoeff - 1.0) * (DotProd(panel1->GetGeoNormal(), panel2->GetGeoNormal()) + 1.0) + 1.0;
	}
	else {
		ccoeff = 1.0f;
	}

	forcerefinement = RefineCriteria(&panel1crit, &panel2crit, panel1, panel2, &dist, rdist, rmax, m_clsGlobalVars.m_dMeshEps, ccoeff);

	if(forcerefinement == true) {

		if(panel1->GetDimension() > panel2->GetDimension()) {
			ret = panel1->Subdivide();

			if(ret != FC_NORMAL_END) {
				return ret;
			}

			DiscretizeMutual((CAutoPanel*)(panel1->m_pLeft), panel2, selfCond);
			DiscretizeMutual((CAutoPanel*)(panel1->m_pRight), panel2, selfCond);
		}
		else {
			ret = panel2->Subdivide();;

			if(ret != FC_NORMAL_END) {
				return ret;
			}
			DiscretizeMutual(panel1, (CAutoPanel*)(panel2->m_pLeft), selfCond);
			DiscretizeMutual(panel1, (CAutoPanel*)(panel2->m_pRight), selfCond);
		}
	}
	else {
		// if panels are standard panels, record max mesh eps that did not lead to refinement
		// (recording max mesh eps for supernodes is not useful,
		// since a supernode can never be a leaf, and so will never lead to an increasingly refined mesh)

		if( !(panel1->m_ucType & AUTOPANEL_IS_SUPER_NODE) ) {
			// if not forced to refine, check and record the max mesh eps
			if(panel1crit > m_dMaxMeshEps) {
				m_dMaxMeshEps = panel1crit;
			}
		}
		if( !(panel2->m_ucType & AUTOPANEL_IS_SUPER_NODE) ) {
			// if not forced to refine, check and record the max mesh eps
			if(panel2crit > m_dMaxMeshEps) {
				m_dMaxMeshEps = panel2crit;
			}
		}
	}

	// return to upper level
	m_iLevel--;

	return ret;
}

// Recursive refinement
int CAutoRefine::DiscretizeMutual(CAutoSegment *panel1, CAutoSegment *panel2, bool selfCond)
{
	bool forcerefinement;
	C2DVector dist;
	double rdist, r1, r2, rmax, ccoeff;
	double panel1crit, panel2crit;
	int ret;

	ret = FC_NORMAL_END;

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	dist = panel1->GetCentroid() - panel2->GetCentroid();

	rdist = Mod(dist);

	// get max lenght of panel side
	r1 = panel1->GetLength();
	r2 = panel2->GetLength();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	// check interaction condition against the threshold
	//

	// use different criteria for conductor and dielectric

	// if mutual interaction between panels belonging to the same conductor,
	// calculate the curvature coefficient, to give a different weight according
	// to the relative angle between the two panels. Parallel panels are treated
	// as if they were farther from each other w.r.t. opposite panels.
	// Remark: we are assuming the correct orientation of the normal of the panels,
	// all pointing outwards from a solid's interior
	if(selfCond == true) {
		ccoeff = (m_clsGlobalVars.m_dMeshCurvCoeff - 1.0) * (DotProd(panel1->GetGeoNormal(), panel2->GetGeoNormal()) + 1.0) + 1.0;
	}
	else {
		ccoeff = 1.0f;
	}

	forcerefinement = RefineCriteria(&panel1crit, &panel2crit, panel1, panel2, &dist, rdist, rmax, m_clsGlobalVars.m_dMeshEps, ccoeff);

	if(forcerefinement == true) {

		if(panel1->GetDimension() > panel2->GetDimension()) {
			ret = panel1->Subdivide();

			if(ret != FC_NORMAL_END) {
				return ret;
			}

			DiscretizeMutual((CAutoSegment*)(panel1->m_pLeft), panel2, selfCond);
			DiscretizeMutual((CAutoSegment*)(panel1->m_pRight), panel2, selfCond);
		}
		else {
			ret = panel2->Subdivide();;

			if(ret != FC_NORMAL_END) {
				return ret;
			}
			DiscretizeMutual(panel1, (CAutoSegment*)(panel2->m_pLeft), selfCond);
			DiscretizeMutual(panel1, (CAutoSegment*)(panel2->m_pRight), selfCond);
		}
	}
	else {
		// if panels are standard panels, record max mesh eps that did not lead to refinement
		// (recording max mesh eps for supernodes is not useful,
		// since a supernode can never be a leaf, and so will never lead to an increasingly refined mesh)

		if( !(panel1->m_ucType & AUTOPANEL_IS_SUPER_NODE) ) {
			// if not forced to refine, check and record the max mesh eps
			if(panel1crit > m_dMaxMeshEps) {
				m_dMaxMeshEps = panel1crit;
			}
		}
		if( !(panel2->m_ucType & AUTOPANEL_IS_SUPER_NODE) ) {
			// if not forced to refine, check and record the max mesh eps
			if(panel2crit > m_dMaxMeshEps) {
				m_dMaxMeshEps = panel2crit;
			}
		}
	}

	// return to upper level
	m_iLevel--;

	return ret;
}

// common criteria for refinement
bool CAutoRefine::RefineCriteria(double *panel1crit, double *panel2crit, CAutoPanel *panel1, CAutoPanel *panel2, C3DVector *dist,
                                 double rdist, double rmax, double eps, double ccoeff)
{
	bool forcerefinement;
	C3DVector normal1, normal2;
	double potestim1, potestim2, dotprod;

	// init the return status
	forcerefinement = false;

	// must distinguish between geometric normal and dielectric normal
	normal1 = panel1->GetDielNormal();
	normal2 = panel2->GetDielNormal();

	// must use numerical potential estimate because using only 1/r estimates gives bad results
	// in the refinement near curvatures w.r.t. flat zones (it is inhomogeneous)
	if(rdist / rmax < 2.0) {
		if( (panel1->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
			potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 2);
		}
		else {
			potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 2);
		}
		if( (panel2->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
			potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 2);
		}
		else {
			potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 2);
		}
	}
	else {
		if( (panel1->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod = DotProd(*dist, normal1);
			potestim1 = dotprod / (rdist * rdist * rdist);
		}
		else {
			potestim1 = 1.0f/rdist;
		}
		if( (panel2->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
			// this is cos((x2-x1),n2) * mod(x1-x2)
			dotprod = DotProd(-(*dist), normal2);
			potestim2 = dotprod / (rdist * rdist * rdist);
		}
		else {
			potestim2 = 1.0f/rdist;
		}
	}

	if( (panel1->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
		*panel1crit = fabs(potestim1 * panel2->GetDimension() * panel2->GetMaxSideLen() / (m_dMaxSide * m_dMaxArea * ccoeff));
		if( *panel1crit > eps) {
			forcerefinement = true;
		}
	}
	else {
		*panel1crit = fabs(potestim1 * panel2->GetDimension() / (m_dMaxArea * ccoeff) );
		if( *panel1crit > eps) {
			forcerefinement = true;
		}
	}
	if((panel2->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
		*panel2crit = fabs(potestim2 * panel1->GetDimension() * panel1->GetMaxSideLen() / (m_dMaxSide * m_dMaxArea * ccoeff));
		if( *panel2crit > eps) {
			forcerefinement = true;
		}
	}
	else {
		*panel2crit = fabs(potestim2 * panel1->GetDimension() / (m_dMaxArea * ccoeff) );
		if( *panel2crit > eps) {
			forcerefinement = true;
		}
	}

	return forcerefinement;
}


// common criteria for refinement
bool CAutoRefine::RefineCriteria(double *panel1crit, double *panel2crit, CAutoSegment *panel1, CAutoSegment *panel2, C2DVector *dist,
                                 double rdist, double rmax, double eps, double ccoeff)
{
	bool forcerefinement;
	C2DVector normal1, normal2;
	double potestim1, potestim2, dotprod;

	// init the return status
	forcerefinement = false;

	// must distinguish between geometric normal and dielectric normal
	normal1 = panel1->GetDielNormal();
	normal2 = panel2->GetDielNormal();

	// must use numerical potential estimate because using only 1/r estimates gives bad results
	// in the refinement near curvatures w.r.t. flat zones (it is inhomogeneous)
	if(rdist / rmax < 2.0) {
		if( (panel1->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
//			potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 2);
			potestim1 = m_clsPotential.EnField(panel1->GetCentroid(), panel2->m_clsVertex, normal1);
		}
		else {
//			potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 2);
//			potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 2);
			potestim1 = m_clsPotential.PotentialOpt(panel1->GetCentroid(), panel2->m_clsVertex);
		}
		if( (panel2->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
//			potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 2);
			potestim2 = m_clsPotential.EnField(panel2->GetCentroid(), panel1->m_clsVertex, normal2);
		}
		else {
//			potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 2);
//			potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 2);
			potestim2 = m_clsPotential.PotentialOpt(panel2->GetCentroid(), panel1->m_clsVertex);
		}
	}
	else {
		if( (panel1->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod = DotProd(*dist, normal1);
			potestim1 = dotprod / (rdist * rdist);
		}
		else {
			// log(1/rdist)
			potestim1 = -log(rdist);
		}
		if( (panel2->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
			// this is cos((x2-x1),n2) * mod(x1-x2)
			dotprod = DotProd(-(*dist), normal2);
			potestim2 = dotprod / (rdist * rdist);
		}
		else {
			// log(1/rdist)
			potestim2 = -log(rdist);
		}
	}

	if( (panel1->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
		*panel1crit = fabs(potestim1 * panel2->GetLength() * panel2->GetLength() / (m_dMaxLength * m_dMaxLength * ccoeff));
		if( *panel1crit > eps) {
			forcerefinement = true;
		}
	}
	else {

		*panel1crit = fabs( (potestim1+0.0) * panel2->GetLength() / (m_dMaxLength * ccoeff) );
		if( *panel1crit > eps) {
			forcerefinement = true;
		}
	}
	if((panel2->m_ucType & AUTOPANEL_IS_DIEL) == AUTOPANEL_IS_DIEL) {
		*panel2crit = fabs(potestim2 * panel1->GetLength() * panel1->GetLength() / (m_dMaxLength * m_dMaxLength * ccoeff));
		if( *panel2crit > eps) {
			forcerefinement = true;
		}
	}
	else {

		*panel2crit = fabs( (potestim2+0.0) * panel1->GetLength() / (m_dMaxLength * ccoeff) );
		if( *panel2crit > eps) {
			forcerefinement = true;
		}
	}

	return forcerefinement;
}


// Recursive refinement
int CAutoRefine::RefineMutual(CAutoPanel *panel1, CAutoPanel *panel2)
{
	bool forcerefinement;
	char refinePanel;
	C3DVector dist;
	double rdist, r1, r2, rmax;
	double panel1crit, panel2crit;
	unsigned long chunk1, posInChunk1, block1;
	unsigned long chunk2, posInChunk2, block2;


	if(g_bFCContinue == false) {
		return FC_USER_BREAK;
	}

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	dist = panel1->GetCentroid() - panel2->GetCentroid();

	rdist = Mod(dist);

	// get max lenght of panel side
	r1 = panel1->GetMaxSideLen();
	r2 = panel2->GetMaxSideLen();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	forcerefinement = false;

	// if both panels were leaf panels, no way: must interact at this level
	if(panel1->IsLeaf() == false || panel2->IsLeaf() == false) {

		// check interaction condition against the threshold;
		// use different criteria for conductor and dielectric
		forcerefinement = RefineCriteria(&panel1crit, &panel2crit, panel1, panel2, &dist, rdist, rmax, m_clsGlobalVars.m_dEps);

// debug only, if you want to force refinement (i.e. creating full link matrix n^2 elements ( - n auto potentials) )
//forcerefinement = true;
	}

	// if must go down in the hierarchy
	if(forcerefinement == true) {

		// if any of the two is a leaf, refine the other one
		if(panel1->IsLeaf() == false && panel2->IsLeaf() == true) {
			refinePanel = 1;
		}
		else if(panel1->IsLeaf() == true && panel2->IsLeaf() == false) {
			refinePanel = 2;
		}
		// otherwise, use panel area as the only criteria to decide which panel gets to be refined
		else  {
			if(panel1->GetDimension() > panel2->GetDimension()) {
				refinePanel = 1;
			}
			else {
				refinePanel = 2;
			}
		}

		// refine the chosen panel (ideally the one producing the biggest coeff. of pot.)
		ASSERT(refinePanel == 1 || refinePanel == 2);
		if(refinePanel == 1) {
			RefineMutual((CAutoPanel*)(panel1->m_pLeft), panel2);
			RefineMutual((CAutoPanel*)(panel1->m_pRight), panel2);
		}
		else {
			RefineMutual(panel1, (CAutoPanel*)(panel2->m_pLeft));
			RefineMutual(panel1, (CAutoPanel*)(panel2->m_pRight));
		}
	}
	// this terminates recursion and stores 'potestim'
	else {

		// if just counting
		if(m_bComputeLinks == false) {

			// increase number of links
			m_ulLinksNum[m_ucInteractionLevel] += 2;

#ifdef DEBUG_DUMP_BASIC
			m_iaLinksBtwLevels[panel1->m_iLevel][panel2->m_iLevel]++;
//			m_ulNumofFastPotest++;
#endif
		}
		else {

			// store interaction information for both panels

			// panel1
			//
			chunk1 = panel1->m_ulLinkIndexEnd[m_ucInteractionLevel] / AUTOREFINE_LINK_CHUNK_SIZE;
			posInChunk1 = panel1->m_ulLinkIndexEnd[m_ucInteractionLevel] % AUTOREFINE_LINK_CHUNK_SIZE;
			// determine in which block the current chunk is ('m_ulLinkChunkNum' is the number of chunks per block)
			block1 = chunk1 / m_ulLinkChunkNum[m_ucInteractionLevel];

			// panel2
			//
			chunk2 = panel2->m_ulLinkIndexEnd[m_ucInteractionLevel] / AUTOREFINE_LINK_CHUNK_SIZE;
			posInChunk2 = panel2->m_ulLinkIndexEnd[m_ucInteractionLevel] % AUTOREFINE_LINK_CHUNK_SIZE;
			// determine in which block the current chunk is ('m_ulLinkChunkNum' is the number of chunks per block)
			block2 = chunk2 / m_ulLinkChunkNum[m_ucInteractionLevel];

			// store interaction
			//

			// only if in current block
			if(block1 == m_ulCurrBlock) {
				// adjust chunk to position within the current block
				chunk1 -= m_ulLinkChunkNum[m_ucInteractionLevel] * m_ulCurrBlock;
				m_pdPanelPtrLinks[m_ucInteractionLevel][chunk1][posInChunk1] = panel2;
			}

			// only if in current block
			if(block2 == m_ulCurrBlock) {
				// adjust chunk to position within the current block
				chunk2 -= m_ulLinkChunkNum[m_ucInteractionLevel] * m_ulCurrBlock;
				m_pdPanelPtrLinks[m_ucInteractionLevel][chunk2][posInChunk2] = panel1;
			}
		}

		// in any case, increment position pointer
		panel1->m_ulLinkIndexEnd[m_ucInteractionLevel]++;
		panel2->m_ulLinkIndexEnd[m_ucInteractionLevel]++;

	}

	// return to upper level
	m_iLevel--;

	return FC_NORMAL_END;
}


// Recursive refinement
int CAutoRefine::RefineMutual(CAutoSegment *panel1, CAutoSegment *panel2)
{
	bool forcerefinement;
	char refinePanel;
	C2DVector dist;
	double potestim1, potestim2;
	double rdist, r1, r2, rmax;
	double panel1crit, panel2crit;
	int potstatus;
	unsigned long chunk1, posInChunk1, block1;
	unsigned long chunk2, posInChunk2, block2;


	if(g_bFCContinue == false) {
		return FC_USER_BREAK;
	}

	// increase depth level
	m_iLevel++;

	// store max recursion level
	if(m_iLevel > m_iMaxLevel)
		m_iMaxLevel = m_iLevel;

	dist = panel1->GetCentroid() - panel2->GetCentroid();

	rdist = Mod(dist);

	// get max lenght of panel side
	r1 = panel1->GetLength();
	r2 = panel2->GetLength();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	forcerefinement = false;

	// if both panels were leaf panels, no way: must interact at this level
	if(panel1->IsLeaf() == false || panel2->IsLeaf() == false) {

		// check interaction condition against the threshold;
		// use different criteria for conductor and dielectric
		forcerefinement = RefineCriteria(&panel1crit, &panel2crit, panel1, panel2, &dist, rdist, rmax, m_clsGlobalVars.m_dEps);

// debug only, if you want to force refinement (i.e. creating full link matrix n^2 elements ( - n auto potentials) )
//forcerefinement = true;

	}

	// if must go down in the hierarchy
	if(forcerefinement == true) {

		// if any of the two is a leaf, refine the other one
		if(panel1->IsLeaf() == false && panel2->IsLeaf() == true) {
			refinePanel = 1;
		}
		else if(panel1->IsLeaf() == true && panel2->IsLeaf() == false) {
			refinePanel = 2;
		}
		// otherwise, use panel area as the only criteria to decide which panel gets to be refined
		else  {
			if(r1 > r2) {
				refinePanel = 1;
			}
			else {
				refinePanel = 2;
			}
		}

		// refine the chosen panel (ideally the one producing the biggest coeff. of pot.)
		ASSERT(refinePanel == 1 || refinePanel == 2);
		if(refinePanel == 1) {
			RefineMutual((CAutoSegment*)(panel1->m_pLeft), panel2);
			RefineMutual((CAutoSegment*)(panel1->m_pRight), panel2);
		}
		else {
			RefineMutual(panel1, (CAutoSegment*)(panel2->m_pLeft));
			RefineMutual(panel1, (CAutoSegment*)(panel2->m_pRight));
		}
	}
	// this terminates recursion and stores 'potestim'
	else {

		// if just counting
		if(m_bComputeLinks == false) {

			// increase number of links
			m_ulLinksNum[m_ucInteractionLevel] += 2;

#ifdef DEBUG_DUMP_BASIC
			m_iaLinksBtwLevels[panel1->m_iLevel][panel2->m_iLevel]++;
//			m_ulNumofFastPotest++;
#endif

			panel1->m_ulLinkIndexEnd[m_ucInteractionLevel]++;
			panel2->m_ulLinkIndexEnd[m_ucInteractionLevel]++;
		}
		else {

// debug
//if(panel1->IsLeaf() == false || panel2->IsLeaf() == false) {
//		forcerefinement = RefineCriteria(&panel1crit, &panel2crit, panel1, panel2, &dist, rdist, rmax, m_clsGlobalVars.m_dEps);
//}

			// store interaction information for both panels

			// panel1
			//
			chunk1 = panel1->m_ulLinkIndexEnd[m_ucInteractionLevel] / AUTOREFINE_LINK_CHUNK_SIZE;
			posInChunk1 = panel1->m_ulLinkIndexEnd[m_ucInteractionLevel] % AUTOREFINE_LINK_CHUNK_SIZE;
			// determine in which block the current chunk is ('m_ulLinkChunkNum' is the number of chunks per block)
			block1 = chunk1 / m_ulLinkChunkNum[m_ucInteractionLevel];

			// panel2
			//
			chunk2 = panel2->m_ulLinkIndexEnd[m_ucInteractionLevel] / AUTOREFINE_LINK_CHUNK_SIZE;
			posInChunk2 = panel2->m_ulLinkIndexEnd[m_ucInteractionLevel] % AUTOREFINE_LINK_CHUNK_SIZE;
			// determine in which block the current chunk is ('m_ulLinkChunkNum' is the number of chunks per block)
			block2 = chunk2 / m_ulLinkChunkNum[m_ucInteractionLevel];

			// if at least one of the two estimates falls into a block, actually calculate the potential estimate
			if(block1 == m_ulCurrBlock || block2 == m_ulCurrBlock) {
				// TBC warning: should change PotEstimateOpt to return always valid potential
				if(m_ucInteractionLevel == AUTOREFINE_HIER_PRE_0_LEVEL) {
					potstatus = PotEstimateOpt(panel1, panel2, potestim1, potestim2);
				}
				else {
					// in this case, we are computing the preconditioner, so don't be too picky about approximation
					potstatus = PotEstimateOpt(panel1, panel2, potestim1, potestim2, AUTOREFINE_PRECOND_HIER);
				}

				_ASSERT(potstatus == 0);

				_ASSERT(!isnan(potestim1));
				_ASSERT(!isnan(potestim2));
				_ASSERT(isfinite(potestim1));
				_ASSERT(isfinite(potestim2));

				if(isnan(potestim1) || !isfinite(potestim1) || isnan(potestim2) || !isfinite(potestim2)) {
					if(m_clsGlobalVars.m_bWarnGivenNaN == false) {
						// signal we already warned the user
						m_clsGlobalVars.m_bWarnGivenNaN = true;
						ErrMsg("Error: mutual-potential calculation failed.\n");
						ErrMsg("       Remark: the precision of the result is affected.\n");
// debug
//						potstatus = PotEstimateOpt(panel1, panel2, potestim1, potestim2, AUTOREFINE_PRECOND_HIER);
					}
					if(m_clsGlobalVars.m_bVerboseOutput == true) {
						if(!isfinite(potestim1) || !isfinite(potestim2)) {
							ErrMsg("Error: mutual-potential calculation gave infinite value\n");
						}
						else {
							ErrMsg("Error: mutual-potential calculation gave 'not a number' value (NaN)\n");
						}
						ErrMsg("       segment 1 data\n");
						panel1->ErrorPrintCoords();
						if(panel1->IsLeaf() == true) {
							ErrMsg("       segment is leaf\n");
						}
						else {
							ErrMsg("       segment is not leaf\n");
						}
						ErrMsg("       segment 2 data\n");
						panel2->ErrorPrintCoords();
						if(panel2->IsLeaf() == true) {
							ErrMsg("       segment is leaf\n");
						}
						else {
							ErrMsg("       segment is not leaf\n");
						}
					}
				}

				_ASSERT(panel1 != NULL);
				_ASSERT(panel2 != NULL);
			}

			// store interaction
			//

			// only if in current block
			if(block1 == m_ulCurrBlock) {
				// adjust chunk to position within the current block
				chunk1 -= m_ulLinkChunkNum[m_ucInteractionLevel] * m_ulCurrBlock;
				m_dPotCoeffLinks[m_ucInteractionLevel][chunk1][posInChunk1] = potestim1;
				m_pdPanelPtrLinks[m_ucInteractionLevel][chunk1][posInChunk1] = panel2;
			}
			// anyway increment position pointer
			panel1->m_ulLinkIndexEnd[m_ucInteractionLevel]++;

			// only if in current block
			if(block2 == m_ulCurrBlock) {
				// adjust chunk to position within the current block
				chunk2 -= m_ulLinkChunkNum[m_ucInteractionLevel] * m_ulCurrBlock;
				m_dPotCoeffLinks[m_ucInteractionLevel][chunk2][posInChunk2] = potestim2;
				m_pdPanelPtrLinks[m_ucInteractionLevel][chunk2][posInChunk2] = panel1;
			}
			// anyway increment position pointer
			panel2->m_ulLinkIndexEnd[m_ucInteractionLevel]++;
		}
	}

	// return to upper level
	m_iLevel--;

	return FC_NORMAL_END;
}

// Potential estimate of element2 on element1
int CAutoRefine::PotEstimateOpt(CAutoElement *element1, CAutoElement *element2, double &potestim1)
{
	int ret;

	if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {
		ret = PotEstimateOpt((CAutoPanel*)element1, (CAutoPanel*)element2, potestim1);
	}
	else if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
		ret = PotEstimateOpt((CAutoSegment*)element1, (CAutoSegment*)element2, potestim1);
	}
	else {
		ASSERT(false);
		ret = FC_GENERIC_ERROR;
	}

	_ASSERT(!isnan(potestim1));
	_ASSERT(isfinite(potestim1));

	return ret;
}


// Potential estimate between two triangular panels
// 'computePrecond' == true will accept also near super-panels
int CAutoRefine::PotEstimateOpt(CAutoPanel *panel1, CAutoPanel *panel2, double &potestim1, unsigned char computePrecond)
{
	C3DVector panel1center, panel2center, dist, normal1, normal2;
	C3DVector normal1l, normal1r, normal2l, normal2r;
	C3DVector dist1l2l, dist1l2r, dist1r2l, dist1r2r;
	C3DVector evalPlus, evalMinus;
	double r1, r2, rmax, rdist;
	double dotprod1;
	double h, potestim_p, potestim_n;
	unsigned char isdiel1;
	bool nearpanels;

	// if auto potential, this is an error
	if(panel1 == panel2) {
		ErrMsg("Internal Error: Trying to calculate auto capacitance with mutual capacitance routine PotEstimateOpt() on panel %x\n", panel1);
		return AUTOREFINE_ERROR_AUTOCAP;
	}

	dist = panel1->GetCentroid() - panel2->GetCentroid();
	rdist = Mod(dist);

	// get max length of panel side
	r1 = panel1->GetMaxSideLen();
	r2 = panel2->GetMaxSideLen();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	isdiel1 = panel1->m_ucType & AUTOPANEL_IS_DIEL;

	// if panel distance is too small, something strange here: they should have
	// already been refined. Possible causes are:
	// - thin panels or overlapping ones
	// - super panels (in case of preconditioner calculation, where refinement is prevented below
	//   a certain level)
	if(rdist < AUTOPANEL_EPS) {
		if(computePrecond == AUTOREFINE_PRECOND_SUPER) {
			if(m_clsGlobalVars.m_bWarnGivenPre == false) {
				// signal we already warned the user
				m_clsGlobalVars.m_bWarnGivenPre = true;
				ErrMsg("Warning: panel distance too small found during potential calculation when forming the preconditioner\n");
				ErrMsg("         In this case, convergence rate may be negatively impacted. If so, do not use the preconditioner.\n");
				ErrMsg("         Remark: the precision of the result is not affected.\n");
			}
			if(m_clsGlobalVars.m_bVerboseOutput == true) {

				ErrMsg("Warning: panel distance %f too small between panel %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}


#ifdef DEBUG_DUMP_OTHER
			if( (issupern1 != AUTOPANEL_IS_SUPER_NODE) || (issupern2 != AUTOPANEL_IS_SUPER_NODE) ) {
				DebugOutputInteracting(panel1, panel2);
			}
#endif
		}
		else {
			ErrMsg("Warning: panel distance too small found during potential calculation\n");
			ErrMsg("         Possible causes are thin panels or, more probably, overlapping panels.\n");

			if(m_clsGlobalVars.m_bVerboseOutput == true) {
				ErrMsg("Warning: panel distance %f too small between panel %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}
		}
	}

	// check if panels are near to each other, or far apart enough to use approx formula 1/r or 1/r^2
	nearpanels = false;

	// if panel1 is conductor
	if( isdiel1 != AUTOPANEL_IS_DIEL ) {

		// let calculate the normal component of the elctric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r^2, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// check if panels are too near to each other
		nearpanels = false;

		if(rdist / rmax < 1.0) {
			// if panels are very near to each other, cannot use a numerical formula for the electric field;
			// since the panels are (at least in part) almost coincident, the integral is almost singular,
			// giving rise to numerical instabilities
			// To solve the issue, we revert to the analytical formula for the potential, and we use
			// divided differences
			// divided differences formula is as per "Multipole-accelerated 3-D capacitance extraction algorithms for
			// structures with conformal dielectrics", formula (6)

			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 17);
			}
			else {
				potestim1 = m_clsPotential.PotentialOpt(panel1->GetCentroid(), panel2->m_clsVertex) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 2.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 10);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 5) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.Mutual_4thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 3) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 10.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.Mutual_2thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 2) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			potestim1 = 1/(FOUR_PI_TIMES_E0 * rdist);
		}

	}
	// if panel1 is dielectric
	else {

		// let calculate the normal component of the elctric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r^2, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal1 = panel1->GetDielNormal();

		// check if panels are too near to each other
		nearpanels = false;

		if(rdist / rmax < 1.0) {
			// if panels are very near to each other, cannot use a numerical formula for the electric field;
			// since the panels are (at least in part) almost coincident, the integral is almost singular,
			// giving rise to numerical instabilities
			// To solve the issue, we revert to the analytical formula for the potential, and we use
			// divided differences
			// divided differences formula is as per "Multipole-accelerated 3-D capacitance extraction algorithms for
			// structures with conformal dielectrics", formula (6)

			// rdist is smaller than rmax, so calculate h from rdist
			h = rdist/20;

			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, h, 17);
			}
			else {

				// evaluation points above and under the centroid
				evalPlus = panel1->GetCentroid() + normal1*h;
				evalMinus = panel1->GetCentroid() - normal1*h;
				// calculate the potentials (remark: division by FOUR_PI_TIMES_E0 done only one time later on)
				potestim_p = m_clsPotential.PotentialOpt(evalPlus, panel2->m_clsVertex);
				potestim_n = m_clsPotential.PotentialOpt(evalMinus, panel2->m_clsVertex);
				// and get the potential
				potestim1 = (potestim_n - potestim_p) / (2*h*FOUR_PI_TIMES_E0);
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 2.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 10);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 5) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 5);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 3) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 10.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualD_2thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex, normal1);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 2) / FOUR_PI_TIMES_E0;
			}

			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod1 = DotProd(dist, normal1);
			potestim1 = dotprod1 / (FOUR_PI_TIMES_E0 * rdist * rdist * rdist);
		}
	}

	// keep track of number of potential estimations
	m_ulNumofpotest += 1;

	_ASSERT(!isnan(potestim1));
	_ASSERT(isfinite(potestim1));

	if(isnan(potestim1) || !isfinite(potestim1) ) {
		if(m_clsGlobalVars.m_bWarnGivenNaN == false) {
			// signal we already warned the user
			m_clsGlobalVars.m_bWarnGivenNaN = true;
			ErrMsg("Error: mutual-potential calculation failed.\n");
			ErrMsg("       Remark: the precision of the result is affected.\n");
		}
		if(m_clsGlobalVars.m_bVerboseOutput == true) {
			if(!isfinite(potestim1) ) {
				ErrMsg("Error: mutual-potential calculation gave infinite value\n");
			}
			else {
				ErrMsg("Error: mutual-potential calculation gave 'not a number' value (NaN)\n");
			}
			ErrMsg("       panel 1 data\n");
			panel1->ErrorPrintCoords();
			if(panel1->IsLeaf() == true) {
				ErrMsg("       panel is leaf\n");
			}
			else {
				ErrMsg("       panel is not leaf\n");
			}
			ErrMsg("       panel 2 data\n");
			panel2->ErrorPrintCoords();
			if(panel2->IsLeaf() == true) {
				ErrMsg("       panel is leaf\n");
			}
			else {
				ErrMsg("       panel is not leaf\n");
			}
		}
		return AUTOREFINE_ERROR_NAN_OR_INF;
	}

	return AUTOREFINE_NO_ERROR;
}

// Potential estimate between two triangular panels
// 'computePrecond' == true will accept also near super-panels
int CAutoRefine::PotEstimateOpt(CAutoPanel *panel1, CAutoPanel *panel2, double &potestim1, double &potestim2, unsigned char computePrecond)
{
	C3DVector panel1center, panel2center, dist, normal1, normal2;
	C3DVector normal1l, normal1r, normal2l, normal2r;
	C3DVector dist1l2l, dist1l2r, dist1r2l, dist1r2r;
	C3DVector evalPlus, evalMinus;
	double r1, r2, rmax, rdist;
	double dotprod1, dotprod2;
	double h, potestim_p, potestim_n;
	unsigned char isdiel1, isdiel2;
	bool nearpanels;

	// if auto potential, this is an error
	if(panel1 == panel2) {
		ErrMsg("Internal Error: Trying to calculate auto capacitance with mutual capacitance routine PotEstimateOpt() on panel %x\n", panel1);
		return AUTOREFINE_ERROR_AUTOCAP;
	}

	dist = panel1->GetCentroid() - panel2->GetCentroid();
	rdist = Mod(dist);

	// get max length of panel side
	r1 = panel1->GetMaxSideLen();
	r2 = panel2->GetMaxSideLen();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	isdiel1 = panel1->m_ucType & AUTOPANEL_IS_DIEL;
	isdiel2 = panel2->m_ucType & AUTOPANEL_IS_DIEL;

	// if panel distance is too small, something strange here: they should have
	// already been refined. Possible causes are:
	// - thin panels or overlapping ones
	// - super panels (in case of preconditioner calculation, where refinement is prevented below
	//   a certain level)
	if(rdist < AUTOPANEL_EPS) {
		if(computePrecond == AUTOREFINE_PRECOND_SUPER) {
			if(m_clsGlobalVars.m_bWarnGivenPre == false) {
				// signal we already warned the user
				m_clsGlobalVars.m_bWarnGivenPre = true;
				ErrMsg("Warning: panel distance too small found during potential calculation when forming the preconditioner\n");
				ErrMsg("         In this case, convergence rate may be negatively impacted. If so, do not use the preconditioner.\n");
				ErrMsg("         Remark: the precision of the result is not affected.\n");
			}
			if(m_clsGlobalVars.m_bVerboseOutput == true) {

				ErrMsg("Warning: panel distance %f too small between panel %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}


#ifdef DEBUG_DUMP_OTHER
			if( (issupern1 != AUTOPANEL_IS_SUPER_NODE) || (issupern2 != AUTOPANEL_IS_SUPER_NODE) ) {
				DebugOutputInteracting(panel1, panel2);
			}
#endif
		}
		else {
			ErrMsg("Warning: panel distance too small found during potential calculation\n");
			ErrMsg("         Possible causes are thin panels or, more probably, overlapping panels.\n");

			if(m_clsGlobalVars.m_bVerboseOutput == true) {
				ErrMsg("Warning: panel distance %f too small between panel %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}
		}
	}

	// check if panels are near to each other, or far apart enough to use approx formula 1/r or 1/r^2
	nearpanels = false;

	// if both panels are conductors
	if( (isdiel1 != AUTOPANEL_IS_DIEL) && (isdiel2 != AUTOPANEL_IS_DIEL) ) {

		if(rdist / rmax < 1.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 17);
				potestim2 = m_clsPotential.MutualHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, 17);
			}
			else {
				potestim1 = m_clsPotential.PotentialOpt(panel1->GetCentroid(), panel2->m_clsVertex) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.PotentialOpt(panel2->GetCentroid(), panel1->m_clsVertex) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 2.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 10);
				potestim2 = m_clsPotential.MutualHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, 10);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 5) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 5) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.Mutual_4thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex);
				potestim2 = m_clsPotential.Mutual_4thOrd_FullNum(panel2->m_clsVertex, panel1->m_clsVertex);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 3) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 3) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 10.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.Mutual_2thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex);
				potestim2 = m_clsPotential.Mutual_2thOrd_FullNum(panel2->m_clsVertex, panel1->m_clsVertex);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 2) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 2) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			potestim1 = 1/(FOUR_PI_TIMES_E0 * rdist);
			potestim2 = potestim1;
		}

#ifdef DEBUG_DUMP_BASIC
		m_iaPotestBtwLevels[panel1->m_iLevel][panel2->m_iLevel]++;
#endif
	}
	// if both panels are dielectrics
	else if( (isdiel1 == AUTOPANEL_IS_DIEL) && (isdiel2 == AUTOPANEL_IS_DIEL) ) {


		// let calculate the normal component of the elctric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r^2, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal1 = panel1->GetDielNormal();
		normal2 = panel2->GetDielNormal();

		if(rdist / rmax < 1.0) {
			// if panels are very near to each other, cannot use a numerical formula for the electric field;
			// since the panels are (at least in part) almost coincident, the integral is almost singular,
			// giving rise to numerical instabilities
			// To solve the issue, we revert to the analytical formula for the potential, and we use
			// divided differences
			// divided differences formula is as per "Multipole-accelerated 3-D capacitance extraction algorithms for
			// structures with conformal dielectrics", formula (6)

			// rdist is smaller than rmax, so calculate h from rdist
			h = rdist/20;

			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, h, 17);
				potestim2 = m_clsPotential.MutualDHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, h, 17);
			}
			else {
				// evaluation points above and under the centroid
				evalPlus = panel1->GetCentroid() + normal1*h;
				evalMinus = panel1->GetCentroid() - normal1*h;
				// calculate the potentials (remark: division by FOUR_PI_TIMES_E0 done only one time later on)
				potestim_p = m_clsPotential.PotentialOpt(evalPlus, panel2->m_clsVertex);
				potestim_n = m_clsPotential.PotentialOpt(evalMinus, panel2->m_clsVertex);
				// and get the potential
				potestim1 = (potestim_n - potestim_p) / (2*h*FOUR_PI_TIMES_E0);

				// evaluation points above and under the centroid
				evalPlus = panel2->GetCentroid() + normal2*h;
				evalMinus = panel2->GetCentroid() - normal2*h;
				// calculate the potentials (remark: division by FOUR_PI_TIMES_E0 done only one time later on)
				potestim_p = m_clsPotential.PotentialOpt(evalPlus, panel1->m_clsVertex);
				potestim_n = m_clsPotential.PotentialOpt(evalMinus, panel1->m_clsVertex);
				// and get the potential
				potestim2 = (potestim_n - potestim_p) / (2*h*FOUR_PI_TIMES_E0);
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 2.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 10);
				potestim2 = m_clsPotential.MutualDFullNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, 10);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 5) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 5) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 5);
				potestim2 = m_clsPotential.MutualDFullNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, 5);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 3) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 3) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 10.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualD_2thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex, normal1);
				potestim2 = m_clsPotential.MutualD_2thOrd_FullNum(panel2->m_clsVertex, panel1->m_clsVertex, normal2);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 2) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 2) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod1 = DotProd(dist, normal1);
			potestim1 = dotprod1 / (FOUR_PI_TIMES_E0 * rdist * rdist * rdist);

			// this is cos((x2-x1),n2) * mod(x1-x2)
			dotprod2 = DotProd(-dist, normal2);
			potestim2 = dotprod2 / (FOUR_PI_TIMES_E0 * rdist * rdist * rdist);
		}
	}
	// if panel1 is conductor and panel2 is dielectric
	else if( (isdiel1 != AUTOPANEL_IS_DIEL) && (isdiel2 == AUTOPANEL_IS_DIEL) ) {

		// let calculate the normal component of the elctric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r^2, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal2 = panel2->GetDielNormal();

		// check if panels are too near to each other
		nearpanels = false;

		if(rdist / rmax < 1.0) {
			// if panels are very near to each other, cannot use a numerical formula for the electric field;
			// since the panels are (at least in part) almost coincident, the integral is almost singular,
			// giving rise to numerical instabilities
			// To solve the issue, we revert to the analytical formula for the potential, and we use
			// divided differences
			// divided differences formula is as per "Multipole-accelerated 3-D capacitance extraction algorithms for
			// structures with conformal dielectrics", formula (6)

			// rdist is smaller than rmax, so calculate h from rdist
			h = rdist/20;

			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 17);
				potestim2 = m_clsPotential.MutualDHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, h, 17);
			}
			else {
				potestim1 = m_clsPotential.PotentialOpt(panel1->GetCentroid(), panel2->m_clsVertex) / FOUR_PI_TIMES_E0;

				// evaluation points above and under the centroid
				evalPlus = panel2->GetCentroid() + normal2*h;
				evalMinus = panel2->GetCentroid() - normal2*h;
				// calculate the potentials (remark: division by FOUR_PI_TIMES_E0 done only one time later on)
				potestim_p = m_clsPotential.PotentialOpt(evalPlus, panel1->m_clsVertex);
				potestim_n = m_clsPotential.PotentialOpt(evalMinus, panel1->m_clsVertex);
				// and get the potential
				potestim2 = (potestim_n - potestim_p) / (2*h*FOUR_PI_TIMES_E0);
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 2.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 10);
				potestim2 = m_clsPotential.MutualDFullNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, 10);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 5) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 5) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.Mutual_4thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex);
				potestim2 = m_clsPotential.MutualDFullNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, 5);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 3) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 3) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 10.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.Mutual_2thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex);
				potestim2 = m_clsPotential.MutualD_2thOrd_FullNum(panel2->m_clsVertex, panel1->m_clsVertex, normal2);
			}
			else {
				potestim1 = m_clsPotential.PotentialNumerical(panel1->GetCentroid(), panel2->m_clsVertex, 2) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.EnFieldNumerical(panel2->GetCentroid(), panel1->m_clsVertex, normal2, 2) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			potestim1 = 1/(FOUR_PI_TIMES_E0 * rdist);

			// this is cos((x2-x1),n2) * mod(x1-x2)
			dotprod2 = DotProd(-dist, normal2);
			potestim2 = dotprod2 / (FOUR_PI_TIMES_E0 * rdist * rdist * rdist);
		}

	}
	// if panel1 is dielectric and panel2 is conductor
	else {

		// let calculate the normal component of the elctric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r^2, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal1 = panel1->GetDielNormal();

		// check if panels are too near to each other
		nearpanels = false;

		if(rdist / rmax < 1.0) {
			// if panels are very near to each other, cannot use a numerical formula for the electric field;
			// since the panels are (at least in part) almost coincident, the integral is almost singular,
			// giving rise to numerical instabilities
			// To solve the issue, we revert to the analytical formula for the potential, and we use
			// divided differences
			// divided differences formula is as per "Multipole-accelerated 3-D capacitance extraction algorithms for
			// structures with conformal dielectrics", formula (6)

			// rdist is smaller than rmax, so calculate h from rdist
			h = rdist / 20.0;

			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, h, 17);
				potestim2 = m_clsPotential.MutualHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, 17);
			}
			else {

				// evaluation points above and under the centroid
				evalPlus = panel1->GetCentroid() + normal1*h;
				evalMinus = panel1->GetCentroid() - normal1*h;
				// calculate the potentials (remark: division by FOUR_PI_TIMES_E0 done only one time later on)
				potestim_p = m_clsPotential.PotentialOpt(evalPlus, panel2->m_clsVertex);
				potestim_n = m_clsPotential.PotentialOpt(evalMinus, panel2->m_clsVertex);
				// and get the potential
				potestim1 = (potestim_n - potestim_p) / (2*h*FOUR_PI_TIMES_E0);

				potestim2 = m_clsPotential.PotentialOpt(panel2->GetCentroid(), panel1->m_clsVertex) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 2.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 10);
				potestim2 = m_clsPotential.MutualHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, 10);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 5) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 5) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 5);
				potestim2 = m_clsPotential.Mutual_4thOrd_FullNum(panel2->m_clsVertex, panel1->m_clsVertex);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 3) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 3) / FOUR_PI_TIMES_E0;
			}
			nearpanels = true;
		}
		else if(rdist / rmax < 10.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
				potestim1 = m_clsPotential.MutualD_2thOrd_FullNum(panel1->m_clsVertex, panel2->m_clsVertex, normal1);
				potestim2 = m_clsPotential.Mutual_2thOrd_FullNum(panel2->m_clsVertex, panel1->m_clsVertex);
			}
			else {
				potestim1 = m_clsPotential.EnFieldNumerical(panel1->GetCentroid(), panel2->m_clsVertex, normal1, 2) / FOUR_PI_TIMES_E0;
				potestim2 = m_clsPotential.PotentialNumerical(panel2->GetCentroid(), panel1->m_clsVertex, 2) / FOUR_PI_TIMES_E0;
			}

			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod1 = DotProd(dist, normal1);
			potestim1 = dotprod1 / (FOUR_PI_TIMES_E0 * rdist * rdist * rdist);

			potestim2 = 1/(FOUR_PI_TIMES_E0 * rdist);
		}
	}

	// keep track of number of potential estimations
	m_ulNumofpotest += 1;

	return AUTOREFINE_NO_ERROR;
}


// Potential estimate between two linear segments
// 'computePrecond' == true will accept also near super-panels
int CAutoRefine::PotEstimateOpt(CAutoSegment *panel1, CAutoSegment *panel2, double &potestim1, unsigned char computePrecond)
{
	C2DVector panel1center, panel2center, dist, normal1;
	C2DVector normal1l, normal1r, normal2l, normal2r;
	C2DVector dist1l2l, dist1l2r, dist1r2l, dist1r2r;
	C2DVector evalPlus, evalMinus;
	double r1, r2, rmax, rdist;
	double dotprod1;
	unsigned char isdiel1;
	bool nearpanels;

	// if auto potential, this is an error
	if(panel1 == panel2) {
		ErrMsg("Internal Error: Trying to calculate auto capacitance with mutual capacitance routine PotEstimateOpt() on segment %x\n", panel1);
		return AUTOREFINE_ERROR_AUTOCAP;
	}

	dist = panel1->GetCentroid() - panel2->GetCentroid();
	rdist = Mod(dist);

	// get length of segments
	r1 = panel1->GetLength();
	r2 = panel2->GetLength();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	isdiel1 = panel1->m_ucType & AUTOPANEL_IS_DIEL;

	// if panel distance is too small, something strange here: they should have
	// already been refined. Possible causes are:
	// - thin panels or overlapping ones
	// - super panels (in case of preconditioner calculation, where refinement is prevented below
	//   a certain level)
	if(rdist < AUTOPANEL_EPS) {
		if(computePrecond == AUTOREFINE_PRECOND_SUPER) {
			if(m_clsGlobalVars.m_bWarnGivenPre == false) {
				// signal we already warned the user
				m_clsGlobalVars.m_bWarnGivenPre = true;
				ErrMsg("Warning: panel distance too small found during potential calculation when forming the preconditioner\n");
				ErrMsg("         In this case, convergence rate may be negatively impacted. If so, do not use the preconditioner.\n");
				ErrMsg("         Remark: the precision of the result is not affected.\n");
			}
			if(m_clsGlobalVars.m_bVerboseOutput == true) {

				ErrMsg("Warning: segment distance %f too small between segments %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}


#ifdef DEBUG_DUMP_OTHER
			if( (issupern1 != AUTOPANEL_IS_SUPER_NODE) || (issupern2 != AUTOPANEL_IS_SUPER_NODE) ) {
				DebugOutputInteracting(panel1, panel2);
			}
#endif
		}
		else {
			ErrMsg("Warning: panel distance too small found during potential calculation\n");
			ErrMsg("         Possible causes are thin panels or, more probably, overlapping panels.\n");

			if(m_clsGlobalVars.m_bVerboseOutput == true) {
				ErrMsg("Warning: segment distance %f too small between segments %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}
		}
	}

	// check if panels are near to each other, or far apart enough to use approx formula 1/r or 1/r^2
	nearpanels = false;

	// if panel1 is a conductor
	if( isdiel1 != AUTOPANEL_IS_DIEL ) {

		if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
//				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 17);
//				potestim2 = m_clsPotential.MutualHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, 17);
			}
			else {
				potestim1 = m_clsPotential.PotentialOpt(panel1->GetCentroid(), panel2->m_clsVertex);
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			potestim1 = -log(rdist);
		}

#ifdef DEBUG_DUMP_BASIC
		m_iaPotestBtwLevels[panel1->m_iLevel][panel2->m_iLevel]++;
#endif
	}
	// if panel1 is a dielectric
	else {

		// let calculate the normal component of the elctric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal1 = panel1->GetDielNormal();

		if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
//					potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 10);
//					potestim2 = m_clsPotential.MutualDFullNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, 10);
			}
			else {
				potestim1 = m_clsPotential.EnField(panel1->GetCentroid(), panel2->m_clsVertex, normal1);
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod1 = DotProd(dist, normal1);
			potestim1 = dotprod1 / (rdist * rdist);
		}
	}

	// keep track of number of potential estimations
	m_ulNumofpotest += 1;

	_ASSERT(!isnan(potestim1));
	_ASSERT(isfinite(potestim1));

	if(isnan(potestim1) || !isfinite(potestim1) ) {
		if(m_clsGlobalVars.m_bWarnGivenNaN == false) {
			// signal we already warned the user
			m_clsGlobalVars.m_bWarnGivenNaN = true;
			ErrMsg("Error: mutual-potential calculation failed.\n");
			ErrMsg("       Remark: the precision of the result is affected.\n");
		}
		if(m_clsGlobalVars.m_bVerboseOutput == true) {
			if(!isfinite(potestim1) ) {
				ErrMsg("Error: mutual-potential calculation gave infinite value\n");
			}
			else {
				ErrMsg("Error: mutual-potential calculation gave 'not a number' value (NaN)\n");
			}
			ErrMsg("       segment 1 data\n");
			panel1->ErrorPrintCoords();
			if(panel1->IsLeaf() == true) {
				ErrMsg("       segment is leaf\n");
			}
			else {
				ErrMsg("       segment is not leaf\n");
			}
			ErrMsg("       segment 2 data\n");
			panel2->ErrorPrintCoords();
			if(panel2->IsLeaf() == true) {
				ErrMsg("       segment is leaf\n");
			}
			else {
				ErrMsg("       segment is not leaf\n");
			}
		}
		return AUTOREFINE_ERROR_NAN_OR_INF;
	}

	return AUTOREFINE_NO_ERROR;
}

// Potential estimate between two linear segments
// 'computePrecond' == true will accept also near super-segments
int CAutoRefine::PotEstimateOpt(CAutoSegment *panel1, CAutoSegment *panel2, double &potestim1, double &potestim2, unsigned char computePrecond)
{
	C2DVector panel1center, panel2center, dist, normal1, normal2;
	C2DVector normal1l, normal1r, normal2l, normal2r;
	C2DVector dist1l2l, dist1l2r, dist1r2l, dist1r2r;
	C2DVector evalPlus, evalMinus;
	double r1, r2, rmax, rdist;
	double dotprod1, dotprod2;
	unsigned char isdiel1, isdiel2;
	bool nearpanels;

	// if auto potential, this is an error
	if(panel1 == panel2) {
		ErrMsg("Internal Error: Trying to calculate auto capacitance with mutual capacitance routine PotEstimateOpt() on segment %x\n", panel1);
		return AUTOREFINE_ERROR_AUTOCAP;
	}

	dist = panel1->GetCentroid() - panel2->GetCentroid();
	rdist = Mod(dist);

	// get length of segments
	r1 = panel1->GetLength();
	r2 = panel2->GetLength();

	if(r1 > r2)
		rmax = r1;
	else
		rmax = r2;

	isdiel1 = panel1->m_ucType & AUTOPANEL_IS_DIEL;
	isdiel2 = panel2->m_ucType & AUTOPANEL_IS_DIEL;

	// if panel distance is too small, something strange here: they should have
	// already been refined. Possible causes are:
	// - thin panels or overlapping ones
	// - super panels (in case of preconditioner calculation, where refinement is prevented below
	//   a certain level)
	if(rdist < AUTOPANEL_EPS) {
		if(computePrecond == AUTOREFINE_PRECOND_SUPER) {
			if(m_clsGlobalVars.m_bWarnGivenPre == false) {
				// signal we already warned the user
				m_clsGlobalVars.m_bWarnGivenPre = true;
				ErrMsg("Warning: panel distance too small found during potential calculation when forming the preconditioner\n");
				ErrMsg("         In this case, convergence rate may be negatively impacted. If so, do not use the preconditioner.\n");
				ErrMsg("         Remark: the precision of the result is not affected.\n");
			}
			if(m_clsGlobalVars.m_bVerboseOutput == true) {

				ErrMsg("Warning: segment distance %f too small between segments %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}


#ifdef DEBUG_DUMP_OTHER
			if( (issupern1 != AUTOPANEL_IS_SUPER_NODE) || (issupern2 != AUTOPANEL_IS_SUPER_NODE) ) {
				DebugOutputInteracting(panel1, panel2);
			}
#endif
		}
		else {
			ErrMsg("Warning: panel distance too small found during potential calculation\n");
			ErrMsg("         Possible causes are thin panels or, more probably, overlapping panels.\n");

			if(m_clsGlobalVars.m_bVerboseOutput == true) {
				ErrMsg("Warning: segment distance %f too small between segments %lx and %lx.\n", rdist, panel1, panel2);
				panel1->ErrorPrintCoords();
				panel2->ErrorPrintCoords();
			}
		}
	}

	// check if panels are near to each other, or far apart enough to use approx formula 1/r or 1/r^2
	nearpanels = false;

	// if both panels are conductors
	if( (isdiel1 != AUTOPANEL_IS_DIEL) && (isdiel2 != AUTOPANEL_IS_DIEL) ) {

		if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
//				potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 17);
//				potestim2 = m_clsPotential.MutualHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, 17);
			}
			else {
				potestim1 = m_clsPotential.PotentialOpt(panel1->GetCentroid(), panel2->m_clsVertex);
				potestim2 = m_clsPotential.PotentialOpt(panel2->GetCentroid(), panel1->m_clsVertex);
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			potestim1 = -log(rdist);
			potestim2 = potestim1;
		}

#ifdef DEBUG_DUMP_BASIC
		m_iaPotestBtwLevels[panel1->m_iLevel][panel2->m_iLevel]++;
#endif
	}
	// if both panels are dielectrics
	else if( (isdiel1 == AUTOPANEL_IS_DIEL) && (isdiel2 == AUTOPANEL_IS_DIEL) ) {

		// let calculate the normal component of the elctric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal1 = panel1->GetDielNormal();
		normal2 = panel2->GetDielNormal();

		if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
//					potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 10);
//					potestim2 = m_clsPotential.MutualDFullNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, 10);
			}
			else {
				potestim1 = m_clsPotential.EnField(panel1->GetCentroid(), panel2->m_clsVertex, normal1);
				potestim2 = m_clsPotential.EnField(panel2->GetCentroid(), panel1->m_clsVertex, normal2);
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod1 = DotProd(dist, normal1);
			potestim1 = dotprod1 / (rdist * rdist);

			// this is cos((x2-x1),n2) * mod(x1-x2)
			dotprod2 = DotProd(-dist, normal2);
			potestim2 = dotprod2 / (rdist * rdist);
		}
	}
	// if panel1 is conductor and panel2 is dielectric
	else if( (isdiel1 != AUTOPANEL_IS_DIEL) && (isdiel2 == AUTOPANEL_IS_DIEL) ) {

		// let calculate the normal component of the electric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r^2, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal2 = panel2->GetDielNormal();

		// check if panels are too near to each other
		nearpanels = false;

		if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
//					potestim1 = m_clsPotential.MutualHalfNumerical(panel1->m_clsVertex, panel2->m_clsVertex, 10);
//					potestim2 = m_clsPotential.MutualDFullNumerical(panel2->m_clsVertex, panel1->m_clsVertex, normal2, 10);
			}
			else {
				potestim1 = m_clsPotential.PotentialOpt(panel1->GetCentroid(), panel2->m_clsVertex);
				potestim2 = m_clsPotential.EnField(panel2->GetCentroid(), panel1->m_clsVertex, normal2);
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			potestim1 = -log(rdist);

			// this is cos((x2-x1),n2) * mod(x1-x2)
			dotprod2 = DotProd(-dist, normal2);
			potestim2 = dotprod2 / (rdist * rdist);
		}

	}
	// if panel1 is dielectric and panel2 is conductor
	else {

		// let calculate the normal component of the electric field on the dielectric
		// panels due to the charge on the generating panels. At infinity, electric
		// field is simply 1/r^2, directed along the line joining the source point
		// and the evaluation point. However, we have to take the component normal to
		// the dielectric panel. This is the reason of the dot product.

		// get normal pointing to outperm side
		normal1 = panel1->GetDielNormal();

		// check if panels are too near to each other
		nearpanels = false;

		if(rdist / rmax < 4.0) {
			if(m_clsGlobalVars.m_cScheme == AUTOREFINE_GALERKIN) {
//				potestim1 = m_clsPotential.MutualDFullNumerical(panel1->m_clsVertex, panel2->m_clsVertex, normal1, 10);
//				potestim2 = m_clsPotential.MutualHalfNumerical(panel2->m_clsVertex, panel1->m_clsVertex, 10);
			}
			else {
				potestim1 = m_clsPotential.EnField(panel1->GetCentroid(), panel2->m_clsVertex, normal1);
				potestim2 = m_clsPotential.PotentialOpt(panel2->GetCentroid(), panel1->m_clsVertex);
			}
			nearpanels = true;
		}

		// if panels are far from each other, use simple potential estimate formula
		if(nearpanels == false) {
			// this is cos((x1-x2),n1) * mod(x1-x2)
			dotprod1 = DotProd(dist, normal1);
			potestim1 = dotprod1 / (rdist * rdist);

			potestim2 = -log(rdist);
		}
	}

	// keep track of number of potential estimations
	m_ulNumofpotest += 1;

	return AUTOREFINE_NO_ERROR;
}

// used only for quadrilateral panles (old version)
double CAutoRefine::MaxSideQ(CAutoPanel *panel)
{
	double lside1, lside2;

	// compute side vectors length
	lside1 = Mod2(panel->m_clsVertex[1]-panel->m_clsVertex[0]);

	lside2 = Mod2(panel->m_clsVertex[2]-panel->m_clsVertex[1]);

	if(lside1 >= lside2)
		return sqrt(lside1);
	else
		return sqrt(lside2);
}

// Reads a FastCap input file, returning a cell array
// of triangular panels
//
// filename	    is the name of the input file
//
// the function returns the number of conductors
//
int CAutoRefine::ReadFastCapFile(CAutoRefGlobalVars *globalVars)
{
	int i;
	long ret;
	clock_t start, finish;
	char fileinname[AR_MAX_PATH], line[AUTOREFINE_MAX_LINE_LEN];
	wxFileName inputFileName;
	bool fileret;
	FILE *fid;
	StlFilePosMap filePosMap;

	// start timer
	start = clock();

	// change the working directory to the one containing
	// the input file, to solve output directory bug when
	// called from FastModel (i.e. the output files, like
	// the refined geometry, is output to the FasterCap launch
	// directory and not to the current directory)
	inputFileName.Assign(globalVars->m_sFileIn);
	fileret = inputFileName.IsOk();
	if(fileret == false) {
		ErrMsg("Filename '%s' is not well defined, aborting\n", globalVars->m_sFileIn.c_str());
		return FC_FILE_ERROR;
	}
	fileret = inputFileName.MakeAbsolute();
	if(fileret == false) {
		ErrMsg("Cannot convert the input file path '%s' to absolute, aborting\n", globalVars->m_sFileIn.c_str());
		return FC_FILE_ERROR;
	}
	// test for existence, and if we have the priviledges to read it
	fileret = inputFileName.IsFileReadable();
	if(fileret == false) {
		ErrMsg("File '%s' is not readable, either because not existing, of you don't have the rights to access the file\n", (const char*)inputFileName.GetFullPath());
		return FC_CANNOT_OPEN_FILE;
	}
	wxSetWorkingDirectory(inputFileName.GetPath());
	// if we now set the working directory, must delete the path part of the filename stored in the globalVars,
	// otherwise next time, if we run multiple iterations, we'll always add the relative part of the path every time
	globalVars->m_sFileIn = inputFileName.GetFullName();


	strncpy(fileinname, (const char*)inputFileName.GetFullPath(), AR_MAX_PATH);
	fileinname[AR_MAX_PATH-1] = '\0';
	// strncpy() gives no way to know if the string to be copied was too long, so the only way to know is to compare the two
	if(strcmp(fileinname, inputFileName.GetFullPath().c_str()) != 0) {
		ErrMsg("Filename '%s' is too long, max number of allowed char in path+filename is %d\n", (const char*)inputFileName.GetFullPath(), AR_MAX_PATH);
		return FC_FILE_ERROR;
	}

	// init member variables
	//
	m_lDielNum = 0;
	m_lCondNum = 0;
	m_iParseLevel = -1;
	for(i=0; i < AUTOREFINE_MAX_PARSE_LEVEL; i++) {
		m_iGroupNum[i] = 1;
	}
	m_lGroupDielNum = 0;
	// init number of input panels
	m_ulInputPanelNum = 0;
	m_ulInputRawCondPanelNum = 0;
	m_ulInputRawDielPanelNum = 0;

	// set default value for dielectric permittivity real or complex value;
	// this setting is reported in the 'globalVars' variable of the caller (passed by reference)
	globalVars->m_ucHasCmplxPerm = AUTOREFINE_REAL_PERM;

	// understand the type of input file (2D or 3D)
	fid = fopen(inputFileName.GetFullPath(), "r");

	if(fid == NULL) {
		ErrMsg("Cannot open file '%s', aborting\n", (const char*)inputFileName.GetFullPath());
		ret = FC_CANNOT_OPEN_FILE;
	}
	else {
		fgets(line, AUTOREFINE_MAX_LINE_LEN, fid);

		fclose(fid);

		filePosMap.clear();

		if(strstr(line, "2D") == NULL && strstr(line, "2d") == NULL) {
			g_ucSolverType = SOLVERGLOBAL_3DSOLVER;
			// Parse3DInputFile can change the 'm_ucHasCmplxPerm', and this must be globally visible;
			// also passing dummy 'fid' and 'filePosMap' (that must be empty) since
			// there is no parent file
			ret = Parse3DInputFile(fileinname, fid, &filePosMap, globalVars);
		}
		else {
			g_ucSolverType = SOLVERGLOBAL_2DSOLVER;
			// init the global model bounding box
			m_clsGlobal2DBbox.Clear();
			// Parse2DInputFile can change the 'm_ucHasCmplxPerm', and this must be globally visible
			// also passing dummy 'fid' and 'filePosMap' (that must be empty) since
			// there is no parent file
			ret = Parse2DInputFile(fileinname, fid, &filePosMap, globalVars);
		}
	}

	m_clsGlobalVars = *globalVars;

	finish = clock();
	m_fDurationReadFile = (float)(finish - start) / CLOCKS_PER_SEC;

	return ret;
}

int CAutoRefine::CreateFileMap(char *fileinname, FILE *fid, StlFilePosMap *filePosMap)
{
	long linenum;
	char line[AUTOREFINE_MAX_LINE_LEN];
	char name[AUTOCONDUCTOR_MAX_NAME_LEN], tmpname[AUTOCONDUCTOR_MAX_NAME_LEN];
	int status, ret;
	StlPosLinenumPair posLinenumPair;
	fpos_t startPos;

	// scan all the file with the file descriptor 'fid' and create the map
	// of the sub-files contained inside

	ret = FC_NORMAL_END;

	if(fid == NULL) {
		ErrMsg("Internal error: CreateFileMap() called with NULL file id\n");
		return FC_CANNOT_OPEN_FILE;
	}

	linenum = 1;
	// get the initial position, to resore at the end of the search
	fgetpos(fid, &startPos);
	fgets(line, AUTOREFINE_MAX_LINE_LEN, fid);

	// read all fastcap file storing panel info
	// according to the their conductors
	while( !feof(fid) && ! ferror(fid) && ret == FC_NORMAL_END) {

		if(g_bFCContinue == false) {
			ret = FC_USER_BREAK;
			break;
		}

		// if beginning of a new sub-file
		if(line[0] == 'F') {

			// read file name
			status = sscanf(line, "%s %s", tmpname, name);
			if(status == EOF || status != 2) {
				ErrMsg("Error: sub-file 'F' definition found, but wrong file name at line %d\n", linenum);
				ErrMsg("       ('%s') of file '%s'\n", name, fileinname);
			}
			else {
				// store the reference
				status = fgetpos(fid, &(posLinenumPair.first));
				if(status != 0)  {
					ErrMsg("Error: cannot get sub-file '%s' position at line %d\n", name, linenum);
					ErrMsg("       of file '%s'\n", fileinname);
					ret = FC_FILE_ERROR;
				}
				posLinenumPair.second = linenum + 1;
				filePosMap->insert(StlFilePosMap::value_type(std::string(name), posLinenumPair));
			}
		}

		fgets(line, AUTOREFINE_MAX_LINE_LEN, fid);
		linenum++;
	}

	if( ferror(fid) ) {
		ErrMsg("Error: file error #%d while reading file '%s', aborting\n", ferror(fid), fileinname);
		ret = FC_FILE_ERROR;
	}

	// restore the initial position
	fsetpos(fid, &startPos);

	return ret;
}

int CAutoRefine::Parse3DInputFile(char *fileinname, FILE *parentFid, StlFilePosMap *parentFilePosMap, CAutoRefGlobalVars *globalVars, bool isdiel, C3DVector offset, double outpermRe, double outpermIm,
                                  const char *groupname, double inpermRe, double inpermIm, C3DVector dielrefpoint)
{
	FILE *fid;
	long linenum, i, j, numsegs;
	int skip, ret, status;
	char line[AUTOREFINE_MAX_LINE_LEN], *point;
	char name[AUTOCONDUCTOR_MAX_NAME_LEN], tmpname[AUTOCONDUCTOR_MAX_NAME_LEN], newname[AUTOCONDUCTOR_MAX_NAME_LEN];
	char condname[AUTOCONDUCTOR_MAX_NAME_LEN], localGroupname[AUTOCONDUCTOR_MAX_NAME_LEN];
	float x, y, z, localOutperm[2], localInperm[2], tmpperm;
	StlAutoCondDeque::iterator itc, itc_old;
	C3DVector localOffset, localDielrefpoint;
	double mod, dist, dir, diag1, diag2, dsides[4], cosangle;
	C3DVector_float vertex[3];
	C3DVector qvertex[4], sides[4], nodes[4], span01, span32;
	C2DVector vertex2d[3];
	bool is_planar, dummyDiel, uselocaldiel;
	C3DPlane plane;
	C3DOperation op;
	char c1, c2, c3;
	C2DTriangulate tri;
	unsigned char dielIndex;
	StlFilePosMap filePosMap, *filePosMapPtr;
	StlFilePosMap::iterator itp;
	fpos_t startPos;

	// increment the recursion level counter
	m_iParseLevel++;

	ret = FC_NORMAL_END;

	itc = m_stlConductors.begin();

	// check if the conductor file is a sub-file
	itp = parentFilePosMap->find(std::string(fileinname));
	if(itp != parentFilePosMap->end()) {
		// if it is a sub-file, copy parent data
		filePosMapPtr = parentFilePosMap;
		fid = parentFid;
		// store current file position to restore it at the end
		status = fgetpos(fid, &startPos);
		if(status != 0)  {
			ErrMsg("Error: cannot get current file '%s' position, aborting\n", fileinname);
			// return without closing the 'fid' since it was opened by the caller
			return FC_FILE_ERROR;
		}
		//  and get linenum and position
		linenum = (*itp).second.second;
		status = fsetpos(fid, &((*itp).second.first));
		if(status != 0)  {
			ErrMsg("Error: cannot set sub-file '%s' position, aborting\n", fileinname);
			// return without closing the 'fid' since it was opened by the caller
			return FC_FILE_ERROR;
		}
	}
	else {
		// open the sub file id
		fid = fopen(fileinname, "rb");
		linenum = 1;

		if(fid == NULL) {
			ErrMsg("Error: cannot open file '%s', aborting\n", fileinname);
			return FC_CANNOT_OPEN_FILE;
		}

		// build the file map (for single input file)
		ret = CreateFileMap(fileinname, fid, &filePosMap);
		if(ret != FC_NORMAL_END) {
			return ret;
		}
		filePosMapPtr = &filePosMap;
	}

	fgets(line, AUTOREFINE_MAX_LINE_LEN, fid);

	// read all fastcap file storing panel info
	// according to the their conductors
	while( !feof(fid) && line[0] != 'E' && line[0] != 'e' && !ferror(fid) &&
	        line[0] != 'F' && line[0] != 'f' && ret == FC_NORMAL_END) {

		if(g_bFCContinue == false) {
			ret = FC_USER_BREAK;
			break;
		}

		// if conductor file
		if(line[0] == 'C') {

			// read file name
			point = line+1;
			sscanf(point, "%s%n", name, &skip);
			point += skip;

			// read outer permittivity

			sscanf(point, "%f%n", &localOutperm[0], &skip);
			point += skip;
			// check if it is a complex value
			status = sscanf(point, " - j %f%n", &localOutperm[1], &skip);
			if(status == 1) {
				// must add the sign
				localOutperm[1] = -localOutperm[1];
				// signal that we have at least one complex permittivity value; must solve as complex
				globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
				point += skip;
			}
			else {
				status = sscanf(point, " + j %f%n", &localOutperm[1], &skip);
				if(status == 1) {
					globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
					point += skip;
				}
				else {
					// local permittivity imaginary part is zero,
					// in case are permittivities are real and not complex values
					localOutperm[1] = 0.0;
				}
			}

			// read offset coordinates
			sscanf(point, "%f%f%f%n", &x, &y, &z, &skip);
			localOffset.pos(x, y, z);
			localOffset += offset;
			point += skip;

			// compute group name (to distinguish between panels with the same
			// conductor name because in the same file called more than once)
			if(m_iParseLevel == 0) {
				strcpy(localGroupname, "g");
			}
			else {
				strcpy(localGroupname, groupname);
			}
			sprintf(tmpname, "%d_", m_iGroupNum[m_iParseLevel]);
			strcat(localGroupname, tmpname);

			// read optional '+'
			tmpname[0] = '\0';
			sscanf(point, "%s%n", tmpname, &skip);
			point += skip;
			// if '+', do not increase group name
			// in this way, panels with the same conductor name in different files
			// will belong to the same conductor
			if(tmpname[0] != '+') {
				// increase group name
				m_iGroupNum[m_iParseLevel]++;
			}

			// recurse into new conductor file
			m_iGroupNum[m_iParseLevel+1] = 1;

			ret = Parse3DInputFile(name, fid, filePosMapPtr, globalVars, false, localOffset, localOutperm[0], localOutperm[1], localGroupname);

			// if any error in Parse3DInputFile() (e.g. out of memory, user break, etc.)
			if(ret != FC_NORMAL_END) {
				break;
			}

		}
		// if dielectric file
		else if(line[0] == 'D') {

			// read file name
			point = line+1;
			sscanf(point, "%s%n", name, &skip);
			point += skip;

			// read outer permittivity
			sscanf(point, "%f%n", &localOutperm[0], &skip);
			point += skip;
			// check if it is a complex value
			status = sscanf(point, " - j %f%n", &localOutperm[1], &skip);
			if(status == 1) {
				// must add the sign
				localOutperm[1] = -localOutperm[1];
				globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
				point += skip;
			}
			else {
				status = sscanf(point, " + j %f%n", &localOutperm[1], &skip);
				if(status == 1) {
					globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
					point += skip;
				}
				else {
					// local permittivity imaginary part is zero,
					// in case are permittivities are real and not complex values
					localOutperm[1] = 0.0;
				}
			}

			// read inner permittivity
			sscanf(point, "%f%n", &localInperm[0], &skip);
			point += skip;
			// check if it is a complex value
			status = sscanf(point, " - j %f%n", &localInperm[1], &skip);
			if(status == 1) {
				// must add the sign
				localInperm[1] = -localInperm[1];
				globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
				point += skip;
			}
			else {
				status = sscanf(point, " + j %f%n", &localInperm[1], &skip);
				if(status == 1) {
					globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
					point += skip;
				}
				else {
					// local permittivity imaginary part is zero,
					// in case are permittivities are real and not complex values
					localInperm[1] = 0.0;
				}
			}

			// read offset coordinates
			sscanf(point, "%f%f%f%n", &x, &y, &z, &skip);
			localOffset.pos(x, y, z);
			localOffset += offset;
			point += skip;

			// read dielectric reference point coordinates
			sscanf(point, "%f%f%f%n", &x, &y, &z, &skip);
			localDielrefpoint.pos(x, y, z);
			localDielrefpoint += offset;
			point += skip;

			// read optional '-'
			tmpname[0] = '\0';
			sscanf(point, "%s%n", tmpname, &skip);
			point += skip;
			// if '-', reverse outperm and inperm;
			// in this way, the reference point is always on the outperm side
			if(tmpname[0] == '-') {
				tmpperm = localOutperm[0];
				localOutperm[0] = localInperm[0];
				localInperm[0] = tmpperm;
				tmpperm = localOutperm[1];
				localOutperm[1] = localInperm[1];
				localInperm[1] = tmpperm;
			}

			dummyDiel = false;
			if(localInperm[0]<= 0.0 || localOutperm[0]<=0.0) {
				dummyDiel = true;
			}
			if(globalVars->m_ucHasCmplxPerm == AUTOREFINE_CPLX_PERM) {
				if(fabs(localInperm[0] - localOutperm[0])/localInperm[0] < AUTOREFINE_REL_DIEL_TOL &&
				        fabs(localInperm[1] - localOutperm[1])/localInperm[1] < AUTOREFINE_REL_DIEL_TOL) {
					dummyDiel = true;
				}
			}
			else {
				if(fabs(localInperm[0] - localOutperm[0])/localInperm[0] < AUTOREFINE_REL_DIEL_TOL) {
					dummyDiel = true;
				}
			}

			if(dummyDiel == true) {
				ErrMsg("Warning: dummy dielectric-dielectric interface found in the input file, skipping\n");
				if(globalVars->m_bVerboseOutput == true) {
					ErrMsg("         Dummy interface defined at line %d of the input file \"%s\"\n", linenum, fileinname);
					ErrMsg("         Dielectric constants are: inperm %f-j%f, outperm %f-j%f; differing less than %f%%\n", localInperm[0], localInperm[1], localOutperm[0], localOutperm[1], AUTOREFINE_REL_DIEL_TOL * 100.0f);
				}
			}
			else {
				// compute dielectric name (to distinguish between panels
				// in the same file called more than once)
				sprintf(localGroupname, "diel%ld", m_lGroupDielNum);
				// increase group name
				m_lGroupDielNum++;

				// recurse into new dielectric file
				ret = Parse3DInputFile(name, fid, filePosMapPtr, globalVars, true, localOffset,  localOutperm[0], localOutperm[1],
				                       localGroupname, localInperm[0], localInperm[1], localDielrefpoint);

				// if any error in Parse3DInputFile() (e.g. out of memory, user break, etc.)
				if(ret != FC_NORMAL_END) {
					break;
				}
			}

		}
		// if triangular patch
		else if(line[0] == 'T') {

			// read conductor name to which the patch belongs
			point = line+1;
			sscanf(point, "%s%n", tmpname, &skip);
			point += skip;

			// read panel coordinates
			for(i=0; i<3; i++) {
				sscanf(point, "%f%f%f%n", &x, &y, &z, &skip);
				vertex[i].pos(x,y,z);
				vertex[i] += offset;

				point += skip;
			}

			// read optional reference point
			status = sscanf(point, "%f%f%f%n", &x, &y, &z, &skip);
			if(status == 3) {
				localDielrefpoint.pos(x,y,z);
				localDielrefpoint += offset;
				uselocaldiel = true;
			}
			else {
				uselocaldiel = false;
			}

			strcpy(name, groupname);
			// if this is a dielectric interface, ignore specific conductor names
			if(isdiel == false) {
				// concat name with group name
				strcat(name, tmpname);
			}

			ret = (long)GetConductor(&(itc), &dielIndex, name, isdiel, outpermRe, outpermIm, inpermRe, inpermIm, dielrefpoint);
			if(ret != FC_NORMAL_END) {
				break;
			}

			// create (triangular) panel
			ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_SIMPLE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
			if(ret != FC_NORMAL_END) {
				break;
			}

			// counting 'raw' panels (i.e. the panel # in the input file, no input refinement, e.g. Q panels
			// split into triangles)
			if(isdiel == false) {
				m_ulInputRawCondPanelNum++;
			}
			else {
				m_ulInputRawDielPanelNum++;
			}
		}
		// if quadrilateral patch
		else if (line[0] == 'Q') {

			// a quadrilateral panel is broken into triangular ones

			// read conductor name to which the patch belongs
			point = line+1;
			sscanf(point, "%s%n", tmpname, &skip);
			point += skip;

			// read panel coordinates
			// Remark: the quadrilateral must be convex!
			for(i=0; i<4; i++) {
				sscanf(point, "%f%f%f%n", &x, &y, &z, &skip);
				qvertex[i].pos(x,y,z);
				qvertex[i] += offset;

				point += skip;
			}

			// read optional reference point
			status = sscanf(point, "%f%f%f%n", &x, &y, &z, &skip);
			if(status == 3) {
				localDielrefpoint.pos(x,y,z);
				localDielrefpoint += offset;
				uselocaldiel = true;
			}
			else {
				uselocaldiel = false;
			}

			// find the quadrilateral polygon plane
			status = op.PlaneFromPolygon(plane, mod, qvertex, 4);
			if(status != C3D_OK) {
				if(status == C3D_POINTS_ARE_COLLINEAR) {
					ErrMsg("Warning: degenerate quadrilateral panel 'Q' found (co-linear vertexes), corner coordinates are:\n");
				}
				else {
					ErrMsg("Warning: generic error found during quadrilateral supporting plane calculation, corner coordinates are:\n");
				}
				ErrMsg("         (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n",
				       qvertex[0].x, qvertex[0].y, qvertex[0].z,
				       qvertex[1].x, qvertex[1].y, qvertex[1].z,
				       qvertex[2].x, qvertex[2].y, qvertex[2].z,
				       qvertex[3].x, qvertex[3].y, qvertex[3].z);
				ErrMsg("         panel definition found at line %d of the input file \"%s\"\n", linenum, fileinname);
				ErrMsg("         Skipping the panel\n");
			}
			else {

				// concat name with group name
				strcpy(name, groupname);
				// if this is a dielectric interface, ignore specific conductor names
				if(isdiel == false) {
					// concat name with group name
					strcat(name, tmpname);
				}

				ret = (long)GetConductor(&itc, &dielIndex, name, isdiel, outpermRe, outpermIm, inpermRe, inpermIm, dielrefpoint);
				if(ret != FC_NORMAL_END) {
					break;
				}

                // only if dumping input geometry, must keep (valid) quadrilateral panels as they are
                if(globalVars->m_bDumpInputGeo == true) {
					// create (quadrilateral) panel
					ret = (long)CreateQPanel(qvertex, tmpname, dielIndex, &itc, fileinname, linenum, globalVars, uselocaldiel, localDielrefpoint, plane);
					if(ret != FC_NORMAL_END) {
						break;
					}
                }
                else {

                    // check for planarity
                    is_planar = true;
                    for(i=0; i<4; i++) {
                        dist = op.PointPlaneDist(qvertex[i], plane);
                        if(fabs(dist) > AUTOPANEL_EPS) {
                            if(globalVars->m_bWarnGivenSkew == false) {
                                // signal we already warned the user about skewed input quadrilaterals
                                globalVars->m_bWarnGivenSkew = true;
                                ErrMsg("Warning: non-planar quadrilateral panels found in the input file\n");
                                ErrMsg("         splitting into two triangles only (no quality triangulation) along the smaller diagonal\n");
                            }
                            if(globalVars->m_bVerboseOutput == true) {
                                ErrMsg("Warning: non-planar quadrilateral panel 'Q' found (skewed vertexes), corner coordinates are:\n");
                                ErrMsg("         (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n",
                                       qvertex[0].x, qvertex[0].y, qvertex[0].z,
                                       qvertex[1].x, qvertex[1].y, qvertex[1].z,
                                       qvertex[2].x, qvertex[2].y, qvertex[2].z,
                                       qvertex[3].x, qvertex[3].y, qvertex[3].z);
                                ErrMsg("         distance of vertex #%d from the supporting plane is %g against max tolerance of %g\n", i+1, fabs(dist), AUTOPANEL_EPS);
                                ErrMsg("         panel definition found at line %d of the input file \"%s\"\n", linenum, fileinname);
                                is_planar = false;
                                break;
                            }
                        }
                    }

                    // if check of planarity was not ok
                    if(is_planar == false) {

                        diag1 = Mod(qvertex[2] - qvertex[0]);
                        diag2 = Mod(qvertex[3] - qvertex[1]);

                        // create first triangle
                        //

                        if(diag1 <= diag2) {
                            vertex[0] = qvertex[0];
                            vertex[1] = qvertex[1];
                            vertex[2] = qvertex[2];
                        }
                        else {
                            vertex[0] = qvertex[0];
                            vertex[1] = qvertex[1];
                            vertex[2] = qvertex[3];
                        }

                        // create (triangular) panels

                        ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_TRIANGULATE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
                        if(ret != FC_NORMAL_END) {
                            break;
                        }


                        // create second triangle
                        //

                        if(diag1 <= diag2) {
                            vertex[0] = qvertex[2];
                            vertex[1] = qvertex[3];
                            vertex[2] = qvertex[0];
                        }
                        else {
                            vertex[0] = qvertex[1];
                            vertex[1] = qvertex[2];
                            vertex[2] = qvertex[3];
                        }

                        // create (triangular) panels

                        ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_TRIANGULATE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
                        if(ret != FC_NORMAL_END) {
                            break;
                        }

                    }
                    // if check of planarity was ok
                    else {
                        //
                        // quality-triangulate the quadrilateral panel
                        //

                        // two options here. If the quadrilateral panel is (almost) a rectangle, divide
                        // it into well formed rectangles and split them in triangles (so the minimum angle is
                        // not too small and we avoid the issues with the quality triangulation of long
                        // thin rectangles: in that case, in fact, if we create two thin triangles splitting
                        // the rectangle in two, then the triangulation recovers the quality of the resulting
                        // triangles but leaves some artifacts in the form of very small triangles)

                        // calculate sides and sides length
                        sides[0] = qvertex[1] - qvertex[0];
                        sides[1] = qvertex[2] - qvertex[1];
                        sides[2] = qvertex[3] - qvertex[2];
                        sides[3] = qvertex[0] - qvertex[3];
                        for(i=0; i<4; i++) {
                            dsides[i] = Mod(sides[i]);
                        }

                        cosangle = DotProd(sides[0], sides[1]) / (dsides[0] * dsides[1]);

                        // if the quadrilateral panel has almost equal length opposite sides
                        // and the angles are around 90 deg (i.e. the quadrilateral is not skewed)
                        if( fabs(dsides[2] - dsides[0]) < dsides[0] * AUTOREFINE_QUAD_SIDE_TOL &&
                                fabs(dsides[3] - dsides[1]) < dsides[1] * AUTOREFINE_QUAD_SIDE_TOL &&
                                fabs(cosangle) < AUTOREFINE_QUAD_ANGLE_TOL ) {

                            // if taller than thinner, rotate corners
                            if(dsides[0] < dsides[1]) {
                                // nodes[0] is used as temp var
                                nodes[0] = qvertex[0];
                                qvertex[0] = qvertex[1];
                                qvertex[1] = qvertex[2];
                                qvertex[2] = qvertex[3];
                                qvertex[3] = nodes[0];

                                nodes[0] = sides[0];
                                sides[0] = sides[1];
                                sides[1] = sides[2];
                                sides[2] = sides[3];
                                sides[3] = nodes[0];

                                // mod is used as temp var
                                mod = dsides[0];
                                dsides[0] = dsides[1];
                                dsides[1] = dsides[2];
                                dsides[2] = dsides[3];
                                dsides[3] = mod;
                            }

                            // must divide side 0 and 2. Remark: casting will always round down to the lower integer,
                            // so we add 0.5 as a trick to round to the closest integer
                            numsegs = (long)(dsides[0] / (dsides[1] * 5.0) + 0.5);

                            nodes[0] = qvertex[0];
                            nodes[3] = qvertex[3];
                            // numsegs = 1 means 'do not subdivide'
                            if(numsegs == 0) {
                                numsegs = 1;
                            }
                            span01 = sides[0] / (double)numsegs;
                            span32 = -sides[2] / (double)numsegs;

                            for(i=0; i<numsegs-1; i++) {

                                nodes[1] = nodes[0] + span01;
                                nodes[2] = nodes[3] + span32;

                                // create first triangle
                                vertex[0] = nodes[0];
                                vertex[1] = nodes[1];
                                vertex[2] = nodes[2];
                                // create (triangular) panels
                                ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_TRIANGULATE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
                                if(ret != FC_NORMAL_END) {
                                    break;
                                }
                                // create second triangle
                                vertex[0] = nodes[0];
                                vertex[1] = nodes[2];
                                vertex[2] = nodes[3];
                                // create (triangular) panels
                                ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_TRIANGULATE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
                                if(ret != FC_NORMAL_END) {
                                    break;
                                }

                                // and move to next one
                                nodes[0] = nodes[1];
                                nodes[3] = nodes[2];
                            }

                            if(ret != FC_NORMAL_END) {
                                break;
                            }

                            // last ones (uses real vertexes to eliminate roundoffs)
                            nodes[1] = qvertex[1];
                            nodes[2] = qvertex[2];

                            // create first triangle
                            vertex[0] = nodes[0];
                            vertex[1] = nodes[1];
                            vertex[2] = nodes[2];
                            // create (triangular) panels
                            ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_TRIANGULATE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
                            if(ret != FC_NORMAL_END) {
                                break;
                            }
                            // create second triangle
                            vertex[0] = nodes[0];
                            vertex[1] = nodes[2];
                            vertex[2] = nodes[3];
                            // create (triangular) panels
                            ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_TRIANGULATE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
                            if(ret != FC_NORMAL_END) {
                                break;
                            }

                        }
                        else {

                            // find a suitable projection plane
                            // (triangulation is 2D!)
                            c3 = op.FindProjPlane(&c1, &c2, plane.m_vecNormal);

                            // Note that the cartesian axes follow the right hand rule, so
                            // a different behaviour is observed for Y and X, Z axes. I.e:
                            //
                            //        Z
                            //        ^  _ Y
                            //        |  /|
                            //        | /
                            //         /----> X
                            //
                            // if normal is along X, projection plane x,y is Y,Z
                            // if normal is along Y, projection plane x,y is -X,Z.
                            // if normal is along Z, projection plane x,y is X,Y.
                            //
                            // in fact, for projection on X,Z (normal along Y), the direction of X
                            // must be reversed, to be coherent with the general rule.
                            // Otherwise, the face loops would be oriented CW in the projection plane,
                            // instead than CCW, so must re-orient them.
                            // To re-orient the loops, simply multiply by 'dir' one of the two coordinates.
                            if(c3 == C3D_DIR_Y) {
                                dir = plane.m_vecNormal[c3] < 0 ? 1.0 : -1.0;
                            }
                            else
                                dir = plane.m_vecNormal[c3] > 0 ? 1.0 : -1.0;

                            //
                            // create vertexes list for Delaunay triangulation
                            //

                            // clean up the previous list (if any)
                            tri.vertexes.resize(0);
                            // and insert the vertexes
                            for(i=0; i<4; i++) {

                                tri.vertexes.push_back( C2DVector(dir * qvertex[i][c1], qvertex[i][c2]) );
                            }

                            // Delaunay triangulate
                            tri.Triangulate();


                            //
                            // Now insert constraints (useful for concave quadrilaterals)
                            //

                            for(i=0; i<4; i++) {

                                if(i<3) {
                                    j = i+1;
                                }
                                else {
                                    j = 0;
                                }

                                tri.InsertConstrEdge(C2DVector(dir * qvertex[i][c1], qvertex[i][c2]), C2DVector(dir * qvertex[j][c1], qvertex[j][c2]));
                            }

                            //
                            // make holes / concavities
                            //

                            for(i=0; i<4; i++) {

                                if(i<3) {
                                    j = i+1;
                                }
                                else {
                                    j = 0;
                                }
                                tri.DeleteCWEdges(C2DVector(dir * qvertex[i][c1], qvertex[i][c2]), C2DVector(dir * qvertex[j][c1], qvertex[j][c2]));
                            }

                            // refine mesh triangles whose smallest angle is < minAngle degrees
                            tri.Refine(AUTOREFINE_MIN_ANGLE);

                            // then at last generate triangulation
                            tri.GenerateTriangles();

    // debug
    //FILE *fp;
    //fp = fopen("triangulation.txt", "w");
    //fprintf(fp, "0 triangulation debug\n\n");

                            // and insert shiny new triangles in front of the face list of
                            // this face's parent shell
                            for(i=0; i<tri.GetTriangleNum(); i++) {

                                tri.GetTriangleCoords(i, &vertex2d[0], &vertex2d[1], &vertex2d[2]);

                                // flip coordinates if required
                                vertex2d[0].x *= dir;
                                vertex2d[1].x *= dir;
                                vertex2d[2].x *= dir;

                                vertex[0] = plane.Point3Dfrom2D(vertex2d[0], c1, c2, c3);
                                vertex[1] = plane.Point3Dfrom2D(vertex2d[1], c1, c2, c3);
                                vertex[2] = plane.Point3Dfrom2D(vertex2d[2], c1, c2, c3);

    // debug
    //fprintf(fp, "T 1  %g %g %g  %g %g %g  %g %g %g\n",
    //        vertex[0].x, vertex[0].y, vertex[0].z,
    //        vertex[1].x, vertex[1].y, vertex[1].z,
    //        vertex[2].x, vertex[2].y, vertex[2].z);


                                // create (triangular) panels

                                ret = (long)CreatePanel(vertex, tmpname, dielIndex, &itc, fileinname, linenum, AUTOREFINE_TRIANGULATE_CREATE_PANEL, globalVars, uselocaldiel, localDielrefpoint);
                                if(ret != FC_NORMAL_END) {
                                    break;
                                }
                            }
    // debug
    //fclose(fp);
                            if(ret != FC_NORMAL_END) {
                                break;
                            }
                        }
                    }
                }
			}
			// counting 'raw' panels (i.e. the panel # in the input file, no input refinement, e.g. Q panels
			// split into triangles)
			if(isdiel == false) {
				m_ulInputRawCondPanelNum++;
			}
			else {
				m_ulInputRawDielPanelNum++;
			}

		}
		// if conductor rename instruction
		else if(line[0] == 'N') {

			// if this is a dielectric interface, ignore specific conductor names
			if(isdiel == false) {
				// read conductor name to rename
				point = line+1;
				sscanf(point, "%s%n", condname, &skip);
				point += skip;

				// read new conductor name
				sscanf(point, "%s%n", tmpname, &skip);
				point += skip;

				// concat name with group name
				strcpy(name, groupname);
				// concat name with group name
				strcat(name, condname);

				// concat newname with group name
				strcpy(newname, groupname);
				// concat newname with group name
				strcat(newname, tmpname);

				// if old and new names are different
				if(strcmp(name, newname) != 0) {
					// if a conductor to be renamed exists
					if( FindConductorFromName(name, itc_old)) {
						// if new conductor name is already existing as well, should merge the two
						if( FindConductorFromName(newname, itc)) {
							MergeConductors(itc, itc_old);
							// make new conductor current: 'itc' already points to the conductor,
							// 'condname' will contain its name
							strcpy(condname, (*itc)->m_sName);
						}
						else {
							// otherwise, simply rename it
							strcpy((*itc_old)->m_sName, newname);
							// make it current
							itc = itc_old;
							// and update also current name, since itc is the current conductor
							strcpy(condname, newname);
						}
					}
				}
			}
		}

		fgets(line, AUTOREFINE_MAX_LINE_LEN, fid);
		linenum++;
	}

	if( ferror(fid) ) {
		ErrMsg("Error %d while reading file '%s', aborting\n", ferror(fid), fileinname);
		ret = FC_FILE_ERROR;
	}

	// close the file ID only if we opened it; if it was inheredited from the caller,
	// the caller will be responsible to close it
	if(itp == parentFilePosMap->end()) {
		fclose(fid);
	}
	else {
		// restore the file position (as of before sub-file parsing)
		status = fsetpos(fid, &startPos);
		if(status != 0)  {
			ErrMsg("Error: cannot set current file '%s' position, aborting\n", fileinname);
			// return without closing the 'fid' since it was opened by the caller
			return FC_FILE_ERROR;
		}
	}

	// decrement the recursion level counter
	m_iParseLevel--;

	return ret;
}


int CAutoRefine::Parse2DInputFile(char *fileinname, FILE *parentFid, StlFilePosMap *parentFilePosMap, CAutoRefGlobalVars *globalVars, bool isdiel, C2DVector offset, double outpermRe, double outpermIm,
                                  const char *groupname, double inpermRe, double inpermIm, C3DVector dielrefpoint)
{
	FILE *fid;
	long linenum, i;
	int skip, ret, status;
	char line[AUTOREFINE_MAX_LINE_LEN], *point;
	char name[AUTOCONDUCTOR_MAX_NAME_LEN], tmpname[AUTOCONDUCTOR_MAX_NAME_LEN], newname[AUTOCONDUCTOR_MAX_NAME_LEN];
	char condname[AUTOCONDUCTOR_MAX_NAME_LEN], localGroupname[AUTOCONDUCTOR_MAX_NAME_LEN];
	float x, y, localOutperm[2], localInperm[2], tmpperm;
	StlAutoCondDeque::iterator itc, itc_old;
	C2DVector localOffset;
	C3DVector localDielrefpoint;
	C2DVector_float vertex[2];
	bool dummyDiel, uselocaldiel;
	unsigned char dielIndex;
	StlFilePosMap filePosMap, *filePosMapPtr;
	StlFilePosMap::iterator itp;
	fpos_t startPos;

	// increment the recursion level counter
	m_iParseLevel++;

	ret = FC_NORMAL_END;

	itc = m_stlConductors.begin();


	// check if the conductor file is a sub-file
	itp = parentFilePosMap->find(std::string(fileinname));
	if(itp != parentFilePosMap->end()) {
		// if it is a sub-file, copy parent data
		filePosMapPtr = parentFilePosMap;
		fid = parentFid;
		// store current file position to restore it at the end
		status = fgetpos(fid, &startPos);
		if(status != 0)  {
			ErrMsg("Error: cannot get current file '%s' position, aborting\n", fileinname);
			// return without closing the 'fid' since it was opened by the caller
			return FC_FILE_ERROR;
		}
		//  and get linenum and position
		linenum = (*itp).second.second;
		status = fsetpos(fid, &((*itp).second.first));
		if(status != 0)  {
			ErrMsg("Error: cannot set sub-file '%s' position, aborting\n", fileinname);
			// return without closing the 'fid' since it was opened by the caller
			return FC_FILE_ERROR;
		}
	}
	else {
		// open the sub file id
		fid = fopen(fileinname, "rb");
		linenum = 1;

		if(fid == NULL) {
			ErrMsg("Error: cannot open file '%s', aborting\n", fileinname);
			return FC_CANNOT_OPEN_FILE;
		}

		// build the file map (for single input file)
		ret = CreateFileMap(fileinname, fid, &filePosMap);
		if(ret != FC_NORMAL_END) {
			return ret;
		}
		filePosMapPtr = &filePosMap;
	}

	fgets(line, AUTOREFINE_MAX_LINE_LEN, fid);

	// read all fastcap file storing panel info
	// according to the their conductors
	while( !feof(fid) && line[0] != 'E' && line[0] != 'e' && !ferror(fid) &&
	        line[0] != 'F' && line[0] != 'f' && ret == FC_NORMAL_END) {

		if(g_bFCContinue == false) {
			ret = FC_USER_BREAK;
			break;
		}

		// if conductor file
		if(line[0] == 'C') {

			// read file name
			point = line+1;
			sscanf(point, "%s%n", name, &skip);
			point += skip;

			// read outer permittivity

			sscanf(point, "%f%n", &localOutperm[0], &skip);
			point += skip;
			// check if it is a complex value
			status = sscanf(point, " - j %f%n", &localOutperm[1], &skip);
			if(status == 1) {
				// must add the sign
				localOutperm[1] = -localOutperm[1];
				// signal that we have at least one complex permittivity value; must solve as complex
				globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
				point += skip;
			}
			else {
				status = sscanf(point, " + j %f%n", &localOutperm[1], &skip);
				if(status == 1) {
					globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
					point += skip;
				}
				else {
					// local permittivity imaginary part is zero,
					// in case are permittivities are real and not complex values
					localOutperm[1] = 0.0;
				}
			}

			// read offset coordinates
			sscanf(point, "%f%f%n", &x, &y, &skip);
			localOffset.pos(x, y);
			localOffset += offset;
			point += skip;

			// compute group name (to distinguish between panels with the same
			// conductor name because in the same file called more than once)
			if(m_iParseLevel == 0) {
				strcpy(localGroupname, "g");
			}
			else {
				strcpy(localGroupname, groupname);
			}
			sprintf(tmpname, "%d_", m_iGroupNum[m_iParseLevel]);
			strcat(localGroupname, tmpname);

			// read optional '+'
			tmpname[0] = '\0';
			sscanf(point, "%s%n", tmpname, &skip);
			point += skip;
			// if '+', do not increase group name
			// in this way, panels with the same conductor name in different files
			// will belong to the same conductor
			if(tmpname[0] != '+') {
				// increase group name
				m_iGroupNum[m_iParseLevel]++;
			}

			// recurse into new conductor file
			m_iGroupNum[m_iParseLevel+1] = 1;
			ret = Parse2DInputFile(name, fid, filePosMapPtr, globalVars, false, localOffset, localOutperm[0], localOutperm[1], localGroupname);

			// if any error in Parse3DInputFile() (e.g. out of memory, user break, etc.)
			if(ret != FC_NORMAL_END) {
				break;
			}

		}
		// if dielectric file
		else if(line[0] == 'D') {

			// read file name
			point = line+1;
			sscanf(point, "%s%n", name, &skip);
			point += skip;

			// read outer permittivity
			sscanf(point, "%f%n", &localOutperm[0], &skip);
			point += skip;
			// check if it is a complex value
			status = sscanf(point, " - j %f%n", &localOutperm[1], &skip);
			if(status == 1) {
				// must add the sign
				localOutperm[1] = -localOutperm[1];
				globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
				point += skip;
			}
			else {
				status = sscanf(point, " + j %f%n", &localOutperm[1], &skip);
				if(status == 1) {
					globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
					point += skip;
				}
				else {
					// local permittivity imaginary part is zero,
					// in case are permittivities are real and not complex values
					localOutperm[1] = 0.0;
				}
			}

			// read inner permittivity
			sscanf(point, "%f%n", &localInperm[0], &skip);
			point += skip;
			// check if it is a complex value
			status = sscanf(point, " - j %f%n", &localInperm[1], &skip);
			if(status == 1) {
				// must add the sign
				localInperm[1] = -localInperm[1];
				globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
				point += skip;
			}
			else {
				status = sscanf(point, " + j %f%n", &localInperm[1], &skip);
				if(status == 1) {
					globalVars->m_ucHasCmplxPerm = AUTOREFINE_CPLX_PERM;
					point += skip;
				}
				else {
					// local permittivity imaginary part is zero,
					// in case are permittivities are real and not complex values
					localInperm[1] = 0.0;
				}
			}

			// read offset coordinates
			sscanf(point, "%f%f%n", &x, &y, &skip);
			localOffset.pos(x, y);
			localOffset += offset;
			point += skip;

			// read dielectric reference point coordinates
			sscanf(point, "%f%f%n", &x, &y, &skip);
			// to have a single def of conductor, we used the 'trick' to store
			// the diel ref point as 2D inside a 3D point, ignoring 'z'
			localDielrefpoint.pos(x, y, 0.0);
			localDielrefpoint.x += offset.x;
			localDielrefpoint.y += offset.y;
			point += skip;

			// read optional '-'
			tmpname[0] = '\0';
			sscanf(point, "%s%n", tmpname, &skip);
			point += skip;
			// if '-', reverse outperm and inperm;
			// in this way, the reference point is always on the outperm side
			if(tmpname[0] == '-') {
				tmpperm = localOutperm[0];
				localOutperm[0] = localInperm[0];
				localInperm[0] = tmpperm;
				tmpperm = localOutperm[1];
				localOutperm[1] = localInperm[1];
				localInperm[1] = tmpperm;
			}

			dummyDiel = false;
			if(localInperm[0]<= 0.0 || localOutperm[0]<=0.0) {
				dummyDiel = true;
			}
			if(globalVars->m_ucHasCmplxPerm == AUTOREFINE_CPLX_PERM) {
				if(fabs(localInperm[0] - localOutperm[0])/localInperm[0] < AUTOREFINE_REL_DIEL_TOL &&
				        fabs(localInperm[1] - localOutperm[1])/localInperm[1] < AUTOREFINE_REL_DIEL_TOL) {
					dummyDiel = true;
				}
			}
			else {
				if(fabs(localInperm[0] - localOutperm[0])/localInperm[0] < AUTOREFINE_REL_DIEL_TOL) {
					dummyDiel = true;
				}
			}

			if(dummyDiel == true) {
				ErrMsg("Warning: dummy dielectric-dielectric interface found in the input file, skipping\n");
				if(globalVars->m_bVerboseOutput == true) {
					ErrMsg("         Dummy interface defined at line %d of the input file \"%s\"\n", linenum, fileinname);
					ErrMsg("         Dielectric constants are: inperm %f-j%f, outperm %f-j%f; differing less than %f%%\n", localInperm[0], localInperm[1], localOutperm[0], localOutperm[1], AUTOREFINE_REL_DIEL_TOL * 100.0f);
				}
			}
			else {
				// compute dielectric name (to distinguish between panels
				// in the same file called more than once)
				sprintf(localGroupname, "diel%ld", m_lGroupDielNum);
				// increase group name
				m_lGroupDielNum++;

				// recurse into new dielectric file
				ret = Parse2DInputFile(name, fid, filePosMapPtr, globalVars, true, localOffset,  localOutperm[0], localOutperm[1],
				                       localGroupname, localInperm[0], localInperm[1], localDielrefpoint);

				// if any error in Parse3DInputFile() (e.g. out of memory, user break, etc.)
				if(ret != FC_NORMAL_END) {
					break;
				}
			}

		}
		// if linear patch (segment)
		else if(line[0] == 'S') {

			// read conductor name to which the patch belongs
			point = line+1;
			sscanf(point, "%s%n", tmpname, &skip);
			point += skip;

			// read panel coordinates
			for(i=0; i<2; i++) {
				sscanf(point, "%f%f%n", &x, &y, &skip);
				vertex[i].pos(x,y);
				vertex[i] += offset;

				// enlarge the bbox
				m_clsGlobal2DBbox += vertex[i];

				point += skip;
			}

			// read optional reference point
			status = sscanf(point, "%f%f%n", &x, &y, &skip);
			if(status == 2) {
				// to have a single def of conductor, we used the 'trick' to store
				// the diel ref point as 2D inside a 3D point, ignoring 'z'
				localDielrefpoint.pos(x, y, 0.0);
				localDielrefpoint.x += offset.x;
				localDielrefpoint.y += offset.y;
				uselocaldiel = true;
			}
			else {
				uselocaldiel = false;
			}

			strcpy(name, groupname);
			// if this is a dielectric interface, ignore specific conductor names
			if(isdiel == false) {
				// concat name with group name
				strcat(name, tmpname);
			}

			ret = (long)GetConductor(&(itc), &dielIndex, name, isdiel, outpermRe, outpermIm, inpermRe, inpermIm, dielrefpoint);
			if(ret != FC_NORMAL_END) {
				break;
			}

			// create segment
			ret = (long)CreateSegment(vertex, tmpname, dielIndex, &itc, fileinname, linenum, globalVars, uselocaldiel, localDielrefpoint);
			if(ret != FC_NORMAL_END) {
				break;
			}

			// counting 'raw' panels (i.e. the panel # in the input file, no input refinement, e.g. Q panels
			// split into triangles)
			if(isdiel == false) {
				m_ulInputRawCondPanelNum++;
			}
			else {
				m_ulInputRawDielPanelNum++;
			}
		}
		// if conductor rename instruction
		else if(line[0] == 'N') {

			// if this is a dielectric interface, ignore specific conductor names
			if(isdiel == false) {
				// read conductor name to rename
				point = line+1;
				sscanf(point, "%s%n", condname, &skip);
				point += skip;

				// read new conductor name
				sscanf(point, "%s%n", tmpname, &skip);
				point += skip;

				// concat name with group name
				strcpy(name, groupname);
				// concat name with group name
				strcat(name, condname);

				// concat newname with group name
				strcpy(newname, groupname);
				// concat newname with group name
				strcat(newname, tmpname);

				// if old and new names are different
				if(strcmp(name, newname) != 0) {
					// if a conductor to be renamed exists
					if( FindConductorFromName(name, itc_old)) {
						// if new conductor name is already existing as well, should merge the two
						if( FindConductorFromName(newname, itc)) {
							MergeConductors(itc, itc_old);
							// make new conductor current: 'itc' already points to the conductor,
							// 'condname' will contain its name
							strcpy(condname, (*itc)->m_sName);
						}
						else {
							// otherwise, simply rename it
							strcpy((*itc_old)->m_sName, newname);
							// make it current
							itc = itc_old;
							// and update also current name, since itc is the current conductor
							strcpy(condname, newname);
						}
					}
				}
			}
		}

		fgets(line, AUTOREFINE_MAX_LINE_LEN, fid);
		linenum++;
	}

	if( ferror(fid) ) {
		ErrMsg("Error %d while reading file '%s', aborting\n", ferror(fid), fileinname);
		ret = FC_FILE_ERROR;
	}

	// close the file ID only if we opened it; if it was inheredited from the caller,
	// the caller will be responsible to close it
	if(itp == parentFilePosMap->end()) {
		fclose(fid);
	}
	else {
		// restore the file position (as of before sub-file parsing)
		status = fsetpos(fid, &startPos);
		if(status != 0)  {
			ErrMsg("Error: cannot set current file '%s' position, aborting\n", fileinname);
			// return without closing the 'fid' since it was opened by the caller
			return FC_FILE_ERROR;
		}
	}

	// decrement the recursion level counter
	m_iParseLevel--;

	return ret;
}

int CAutoRefine::CreatePanel(C3DVector_float vertex[3], char *nakedname, unsigned char dielIndex,
                             StlAutoCondDeque::iterator *itc, char *fileinname, long linenum, char opType, CAutoRefGlobalVars *globalVars,
                             bool uselocaldiel, C3DVector &localDielrefpoint)

{
	int i;
	float tmpref;
	double cosmin;
	CAutoPanel *newpanel;

	// create new panel
	// SAFENEW_RET(TYPE, VAR, MEM)
	SAFENEW_RET(CAutoPanel, newpanel, g_clsMemUsage.m_ulPanelsMem)

	// store panel coordinates
	for(i=0; i<3; i++) {
		newpanel->m_clsVertex[i] = vertex[i];
	}

	// compute panel geometrical parameters
	cosmin = newpanel->CalcPanelGeomPar();

	// if dielectric panel, calculate outperm side with respect to normal
	if((**itc)->m_bIsDiel == true) {
		if(uselocaldiel == true) {
			tmpref = DotProd( ( localDielrefpoint - newpanel->m_clsVertex[0]), newpanel->m_clsNormal );
			// mark the panel to remember that the outperm direction was calculated w.r.t. a local panel indication,
			// and not the global conductor refpoint definition
			newpanel->m_ucType |= AUTOPANEL_OUTPERM_ELEMENT_LEVEL;
		}
		else {
			tmpref = DotProd( ( (**itc)->m_clsDielRef3DPoint - newpanel->m_clsVertex[0]), newpanel->m_clsNormal );
		}

		if( tmpref > 0 ) {
			newpanel->m_ucType |= (AUTOPANEL_IS_DIEL | AUTOPANEL_OUTPERM_NORMAL_DIR);
		}
		else {
			newpanel->m_ucType |= AUTOPANEL_IS_DIEL;
		}
	}
	// otherwise, conductor panel; store a reference to the dielectric constant index
	else {
		newpanel->m_ucDielIndex = dielIndex;
	}

	// cosmin is the cosinus of the minimum angle of the triangle; check if it is too thin
	if(cosmin >= AUTOREFINE_COS_MIN) {
		if(globalVars->m_bWarnGivenThin == false) {
			// signal we already warned the user about thin triangles
			globalVars->m_bWarnGivenThin = true;
			ErrMsg("Warning: thin triangles (min angle less than %g degrees) ", acos(AUTOREFINE_COS_MIN) * 180.0 / PI);
			if(opType == AUTOREFINE_TRIANGULATE_CREATE_PANEL) {
				ErrMsg("created in the triangulation of a quadrilateral panel\n");
			}
			else {
				ErrMsg("found in the input file\n");
			}
			ErrMsg("         This may impact the precision of the result due to numerical rounding errors\n");
		}
		if(globalVars->m_bVerboseOutput == true) {
			// signal we already warned the user about thin triangles
			ErrMsg("Warning: thin triangle found, corner coordinates are:\n");
			ErrMsg("         (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n",
			       newpanel->m_clsVertex[0].x, newpanel->m_clsVertex[0].y, newpanel->m_clsVertex[0].z,
			       newpanel->m_clsVertex[1].x, newpanel->m_clsVertex[1].y, newpanel->m_clsVertex[1].z,
			       newpanel->m_clsVertex[2].x, newpanel->m_clsVertex[2].y, newpanel->m_clsVertex[2].z );
			ErrMsg("         min angle is: %g\n", acos(cosmin) * 180.0 / PI);
			ErrMsg("         conductor name in the input file is %s, decorated conductor name is %s\n", nakedname, (**itc)->m_sName);
			if(opType == AUTOREFINE_TRIANGULATE_CREATE_PANEL) {
				ErrMsg("         the thin triangle was created in the triangulation of the quadrilateral panel found");
			}
			else {
				ErrMsg("         thin triangle panel definition found");
			}
			ErrMsg(" at line %d of the input file \"%s\"\n", linenum, fileinname);

		}
	}

	// increase panel counter
	m_ulInputPanelNum++;

	// now itc points to the correct conductor,
	// can add the panel to panel list
	(**itc)->m_stlPanels.push_back(newpanel);
	(**itc)->m_ulInputPanelNum++;
#ifdef DEBUG_DUMP_BASIC
	// store reference to conductor
	newpanel->m_pCond = **itc;
#endif

	return FC_NORMAL_END;
}

// create quadrilateral panel (only for input dumping)
int CAutoRefine::CreateQPanel(C3DVector qvertex[4], char *nakedname, unsigned char dielIndex,
                             StlAutoCondDeque::iterator *itc, char *fileinname, long linenum, CAutoRefGlobalVars *globalVars,
                             bool uselocaldiel, C3DVector &localDielrefpoint, C3DPlane plane)

{
    int i;
	float tmpref;
	double cosmin, dist;
	CAutoQPanel *newpanel;
    C3DOperation op;

	// create new panel
	// SAFENEW_RET(TYPE, VAR, MEM)
	SAFENEW_RET(CAutoQPanel, newpanel, g_clsMemUsage.m_ulPanelsMem)

	// store panel coordinates
	for(i=0; i<4; i++) {
		newpanel->m_clsQVertex[i] = qvertex[i];
	}

    // check for planarity
    for(i=0; i<4; i++) {
        dist = op.PointPlaneDist(qvertex[i], plane);
        if(fabs(dist) > AUTOPANEL_EPS) {
            if(globalVars->m_bWarnGivenSkew == false) {
                // signal we already warned the user about skewed input quadrilaterals
                globalVars->m_bWarnGivenSkew = true;
                ErrMsg("Warning: non-planar quadrilateral panels found in the input file\n");
                ErrMsg("         Panel geometrical parameters (e.g. normal) will be incorrect\n");
            }
            if(globalVars->m_bVerboseOutput == true) {
                ErrMsg("Warning: non-planar quadrilateral panel 'Q' found (skewed vertexes), corner coordinates are:\n");
                ErrMsg("         (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n",
                       qvertex[0].x, qvertex[0].y, qvertex[0].z,
                       qvertex[1].x, qvertex[1].y, qvertex[1].z,
                       qvertex[2].x, qvertex[2].y, qvertex[2].z,
                       qvertex[3].x, qvertex[3].y, qvertex[3].z);
                ErrMsg("         distance of vertex #%d from the supporting plane is %g against max tolerance of %g\n", i+1, fabs(dist), AUTOPANEL_EPS);
                ErrMsg("         panel definition found at line %d of the input file \"%s\"\n", linenum, fileinname);
                break;
            }
        }
    }

	// compute panel geometrical parameters
	cosmin = newpanel->CalcPanelGeomPar();

	// if dielectric panel, calculate outperm side with respect to normal
	if((**itc)->m_bIsDiel == true) {
		if(uselocaldiel == true) {
			tmpref = DotProd( ( localDielrefpoint - newpanel->m_clsQVertex[0]), newpanel->m_clsNormal );
			// mark the panel to remember that the outperm direction was calculated w.r.t. a local panel indication,
			// and not the global conductor refpoint definition
			newpanel->m_ucType |= AUTOPANEL_OUTPERM_ELEMENT_LEVEL;
		}
		else {
			tmpref = DotProd( ( (**itc)->m_clsDielRef3DPoint - newpanel->m_clsQVertex[0]), newpanel->m_clsNormal );
		}

		if( tmpref > 0 ) {
			newpanel->m_ucType |= (AUTOPANEL_IS_DIEL | AUTOPANEL_OUTPERM_NORMAL_DIR);
		}
		else {
			newpanel->m_ucType |= AUTOPANEL_IS_DIEL;
		}
	}
	// otherwise, conductor panel; store a reference to the dielectric constant index
	else {
		newpanel->m_ucDielIndex = dielIndex;
	}

	// cosmin is the cosinus of the minimum angle of the quadrilateral; check if it is too thin
	if(cosmin >= AUTOREFINE_COS_MIN) {
		if(globalVars->m_bWarnGivenThin == false) {
			// signal we already warned the user about thin triangles
			globalVars->m_bWarnGivenThin = true;
			ErrMsg("Warning: thin quadrilateral (min angle less than %g degrees) found in the input file\n", acos(AUTOREFINE_COS_MIN) * 180.0 / PI);
		}
		if(globalVars->m_bVerboseOutput == true) {
			// signal we already warned the user about thin triangles
			ErrMsg("Warning: thin quadrilateral found, corner coordinates are:\n");
            ErrMsg("         (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n",
                   qvertex[0].x, qvertex[0].y, qvertex[0].z,
                   qvertex[1].x, qvertex[1].y, qvertex[1].z,
                   qvertex[2].x, qvertex[2].y, qvertex[2].z,
                   qvertex[3].x, qvertex[3].y, qvertex[3].z);
			ErrMsg("         min angle is: %g\n", acos(cosmin) * 180.0 / PI);
			ErrMsg("         conductor name in the input file is %s, decorated conductor name is %s\n", nakedname, (**itc)->m_sName);
			ErrMsg("         thin triangle panel definition found at line %d of the input file \"%s\"\n", linenum, fileinname);
		}
	}

	// increase panel counter
	m_ulInputPanelNum++;

	// now itc points to the correct conductor,
	// can add the panel to panel list
	(**itc)->m_stlPanels.push_back(newpanel);
	(**itc)->m_ulInputPanelNum++;
#ifdef DEBUG_DUMP_BASIC
	// store reference to conductor
	newpanel->m_pCond = **itc;
#endif

	return FC_NORMAL_END;
}

int CAutoRefine::CreateSegment(C2DVector_float vertex[2], char *nakedname, unsigned char dielIndex,
                               StlAutoCondDeque::iterator *itc, char *fileinname, long linenum, CAutoRefGlobalVars *globalVars,
                               bool uselocaldiel, C3DVector &localDielrefpoint)

{
	int i;
	float tmpref;
	double length;
	CAutoSegment *newpanel;
	C2DVector dielrefpoint2D;

	// create new panel
	// SAFENEW_RET(TYPE, VAR, MEM)
	SAFENEW_RET(CAutoSegment, newpanel, g_clsMemUsage.m_ulPanelsMem)

	// store panel coordinates
	for(i=0; i<2; i++) {
		newpanel->m_clsVertex[i] = vertex[i];
	}

	// compute panel geometrical parameters
	newpanel->CalcSegmentGeomPar();

	// if dielectric panel, calculate outperm side with respect to normal
	if((**itc)->m_bIsDiel == true) {
		// copy in 2D (to have a single def of conductor, we used the 'trick' to store
		// the diel ref point as 2D inside a 3D point, ignoring 'z')
		if(uselocaldiel == true) {
			dielrefpoint2D.x = localDielrefpoint.x;
			dielrefpoint2D.y = localDielrefpoint.y;
		}
		else {
			dielrefpoint2D.x = (**itc)->m_clsDielRef3DPoint.x;
			dielrefpoint2D.y = (**itc)->m_clsDielRef3DPoint.y;
		}
		tmpref = DotProd( (dielrefpoint2D - newpanel->m_clsVertex[0]), newpanel->m_clsNormal );
		if( tmpref > 0 ) {
			newpanel->m_ucType |= (AUTOPANEL_IS_DIEL | AUTOPANEL_OUTPERM_NORMAL_DIR);
		}
		else {
			newpanel->m_ucType |= AUTOPANEL_IS_DIEL;
		}
	}
	// otherwise, conductor panel; store a reference to the dielectric constant index
	else {
		newpanel->m_ucDielIndex = dielIndex;
	}

	// check if dummy segment
	length = Mod(newpanel->m_clsVertex[1] - newpanel->m_clsVertex[0]);
	if( length < AUTOREFINE_MIN_LEN) {
		if(globalVars->m_bWarnGivenThin == false) {
			// signal we already warned the user about dummy segments
			globalVars->m_bWarnGivenThin = true;
			ErrMsg("Warning: very short segment found in the input file (length less than %g) ", AUTOREFINE_MIN_LEN);
			ErrMsg("         This may impact the precision of the result due to numerical errors\n");
		}
		if(globalVars->m_bVerboseOutput == true) {
			// signal we already warned the user about thin triangles
			ErrMsg("Warning: very short segment found, end point coordinates are:\n");
			ErrMsg("         (%g,%g) (%g,%g)\n",
			       newpanel->m_clsVertex[0].x, newpanel->m_clsVertex[0].y,
			       newpanel->m_clsVertex[1].x, newpanel->m_clsVertex[1].y);
			ErrMsg("         length is: %g\n", length);
			ErrMsg("         conductor name in the input file is %s, decorated conductor name is %s\n", nakedname, (**itc)->m_sName);
			ErrMsg("         very short segment definition found at line %d of the input file \"%s\"\n", linenum, fileinname);
		}
	}

	// increase panel counter
	m_ulInputPanelNum++;

	// now itc points to the correct conductor,
	// can add the panel to panel list
	(**itc)->m_stlSegments.push_back(newpanel);
	(**itc)->m_ulInputPanelNum++;
#ifdef DEBUG_DUMP_BASIC
	// store reference to conductor
	newpanel->m_pCond = **itc;
#endif

	return FC_NORMAL_END;
}

// in case of conductor (and not dielectric) will return in 'dielIndex' the index to the element of the permittivity array,
// stored in the target CAutoConductor object, that contains the permittivity value for the considered surface
int CAutoRefine::GetConductor(StlAutoCondDeque::iterator *itc, unsigned char *dielIndex, char *name,
                              bool isdiel, double outpermRe, double outpermIm, double inpermRe, double inpermIm, C3DVector &dielrefpoint)

{
	CAutoConductor *newconductor;
	bool createnewcond;

	createnewcond = false;

	// if no conductors in the list yet
	if( (*itc) == m_stlConductors.end()) {
		createnewcond = true;
	}
	// else if current conductor name is not matching the 'name'
	else if( strcmp( ((**itc)->m_sName), name) != 0) {
		// check if not already encountered; if yes, 'itc' contains the reference
		if(!FindConductorFromName(name, (*itc))) {
			// otherwise must create
			createnewcond = true;
		}
	}

	if(	createnewcond == true) {
		// create and add to the list
		if(isdiel == false) {
			// SAFENEW_CTOR_RET(TYPE, PARAM, VAR, MEM)
			SAFENEW_CTOR_RET(CAutoConductor, name DEF_COMMA false DEF_COMMA outpermRe DEF_COMMA outpermIm, newconductor, g_clsMemUsage.m_ulCondMem)
			m_lCondNum++;

			// conductors are pushed to the back of the deque.
			// This is needed for the algorithm in 2D capacitance extraction
			// that enforces the zero total charge condition
			m_stlConductors.push_back(newconductor);
			// make it current
			(*itc) = m_stlConductors.end();
			(*itc)--;
		}
		else {
			// SAFENEW_CTOR_RET(TYPE, PARAM, VAR, MEM)
			SAFENEW_CTOR_RET(CAutoConductor, name DEF_COMMA true DEF_COMMA outpermRe DEF_COMMA outpermIm DEF_COMMA inpermRe DEF_COMMA inpermIm DEF_COMMA dielrefpoint, newconductor, g_clsMemUsage.m_ulCondMem)
			m_lDielNum++;

			// dielectrics are pushed to the front of the deque.
			// This is needed for the algorithm in 2D capacitance extraction
			// that enforces the zero total charge condition
			m_stlConductors.push_front(newconductor);
			// make it current
			(*itc) = m_stlConductors.begin();
		}

	}

	// now the conductor is existing, if it is a real conductor and not a dielectric i/f,
	// assign the correct dielectric constant index
	if(isdiel == false) {
		for(*dielIndex = 0; *dielIndex < (**itc)->m_ucMaxSurfOutperm; (*dielIndex)++) {
			if( (**itc)->m_dSurfOutperm[*dielIndex][0] == outpermRe && (**itc)->m_dSurfOutperm[*dielIndex][1] == outpermIm) {
				// found
				break;
			}
		}
		// if not found
		if(*dielIndex == (**itc)->m_ucMaxSurfOutperm) {
			// if there is still space in the array
			if((**itc)->m_ucMaxSurfOutperm < AUTOPANEL_MAX_DIEL_NUM - 1) {
				(**itc)->m_dSurfOutperm[*dielIndex][0] = outpermRe;
				(**itc)->m_dSurfOutperm[*dielIndex][1] = outpermIm;
				(**itc)->m_ucMaxSurfOutperm++;
			}
			// otherwise raise error
			else {
				(*dielIndex)--;
				ErrMsg("Error: maximum number (%d) of mediums with different dielectric constants\n", AUTOPANEL_MAX_DIEL_NUM-1);
				ErrMsg("       surrounding the single conductor %s reached\n", (**itc)->m_sName);
				// if real value
				if( (**itc)->m_dSurfOutperm[*dielIndex][1] == 0.0 && outpermIm == 0.0) {
					ErrMsg("       Continuing using permittivity value %g instead of %g\n", (**itc)->m_dSurfOutperm[*dielIndex][0], outpermRe);
				}
				else {
					ErrMsg("       Continuing using permittivity value %g+j%g instead of %g+j%g\n", (**itc)->m_dSurfOutperm[*dielIndex][0], (**itc)->m_dSurfOutperm[*dielIndex][1], outpermRe, outpermIm);
				}
			}
		}
	}

	return FC_NORMAL_END;
}

bool CAutoRefine::FindConductorFromName(char *name, StlAutoCondDeque::iterator &itc)
{
	// scan all conductors
	for(itc = m_stlConductors.begin(); itc != m_stlConductors.end(); itc++) {
		// if conductor name matches, conductor found
		if(strcmp(name, (*itc)->m_sName) == 0)
			return true;
	}
	// no conductor with name 'name'
	return false;
}

void CAutoRefine::MergeConductors(StlAutoCondDeque::iterator itc, StlAutoCondDeque::iterator itc_old)
{
	CAutoConductor *tmpcond;

#ifdef DEBUG_DUMP_OTHER
	// must change the reference to the parent conductor to all panels going to be merged
	StlAutoPanelDeque::iterator itp;
	// scan all panels
	for(itp = (*itc_old)->m_stlPanels.begin(); itp != (*itc_old)->m_stlPanels.end(); itp++) {
		(*itp)->m_pCond = itc;
	}
#endif

	// copy panel list from old conductor to new conductor
	// TBC warning: not sure if the 'insert' deque operation is very fast for large number of panels
	(*itc)->m_stlPanels.insert((*itc)->m_stlPanels.end(),(*itc_old)->m_stlPanels.begin(),(*itc_old)->m_stlPanels.end());
	// sum the panel number
	(*itc)->m_ulInputPanelNum += (*itc_old)->m_ulInputPanelNum;

	// save a temporary reference to the old conductor
	tmpcond = *itc_old;
	// remove it from the conductors deque
	m_stlConductors.erase(itc_old);
	// delete the conductor itself
	delete tmpcond;
	g_clsMemUsage.m_ulCondMem -= sizeof(CAutoConductor);
	// and decrement the number of conductors accordingly
	m_lCondNum--;

}

int CAutoRefine::BuildSuperHierarchy()
{
	StlAutoCondDeque::iterator itc1;
	StlAutoPanelDeque::iterator itp1;
	StlAutoSegmentDeque::iterator its1;
	unsigned long i;
	clock_t start, finish;
	double scale;

	// start timer to time superhierarchy building
	start = clock();

	m_dMaxArea = 0.0;
	m_dMaxSide = 0.0;
	m_dMaxLength = 0.0;
	m_dTotalArea = 0.0;

	// scan all conductor groups
	for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {

		// build differently for segments or panels
		if(g_ucSolverType == SOLVERGLOBAL_3DSOLVER) {

			// declare the array of references to panels used to divide
			// and sort panels into tree branches
			// remark: this is only tested for out of memory and not used in memory
			// count, since the array is temporary only, local to this function
			// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
			SAFENEW_ARRAY_NOMEM_RET(CAutoPanel *, m_pPanelArray1, (*itc1)->m_ulInputPanelNum)
			SAFENEW_ARRAY_NOMEM_RET(CAutoPanel *, m_pPanelArray2, (*itc1)->m_ulInputPanelNum)

			// scan every panel inside the given conductor group
			for(itp1=(*itc1)->m_stlPanels.begin(), i=0; itp1!=(*itc1)->m_stlPanels.end(); itp1++, i++) {

				// update conductor bbox
				(*itc1)->m_cls3DBbox += (*itp1)->GetCentroid();

				// fill panel array
				m_pPanelArray1[i] = *itp1;
			}

			// now build super hierarchy
			m_iLevel = -1;
			(*itc1)->m_uTopElement.m_pTopPanel = RecurBuild3DSuperHier(0, (*itc1)->m_ulInputPanelNum, &((*itc1)->m_cls3DBbox));
			// if out of memory
			if((*itc1)->m_uTopElement.m_pTopPanel == NULL) {
				return FC_OUT_OF_MEMORY;
			}

			// clear structures of linked list of panels (as we have now a hierarchical structure
			// with a single top panel)
			(*itc1)->m_stlPanels.clear();

			delete [] m_pPanelArray1;
			delete [] m_pPanelArray2;

			// compute max area and diameter (used to normalize panel areas)
			// TBC warning: for optimization, should normalize in the eps constant, not in panel areas!

			// keep track of maximum side of any panel in the model
			if( (*itc1)->m_uTopElement.m_pTopPanel->GetMaxSideLen() > m_dMaxSide) {
				m_dMaxSide = (*itc1)->m_uTopElement.m_pTopPanel->GetMaxSideLen();
			}

			// keep track of maximum area of any panel in the model
			if( (*itc1)->m_uTopElement.m_pTopPanel->GetDimension() > m_dMaxArea) {
				m_dMaxArea = (*itc1)->m_uTopElement.m_pTopPanel->GetDimension();
			}

			m_dTotalArea += (*itc1)->m_uTopElement.m_pTopPanel->GetDimension();
		}
		else if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {

			// declare the array of references to segments used to divide
			// and sort segments into tree branches
			// remark: this is only tested for out of memory and not used in memory
			// count, since the array is temporary only, local to this function
			// SAFENEW_ARRAY_NOMEM_RET(TYPE, VAR, LEN)
			SAFENEW_ARRAY_NOMEM_RET(CAutoSegment *, m_pSegmentArray1, (*itc1)->m_ulInputPanelNum)
			SAFENEW_ARRAY_NOMEM_RET(CAutoSegment *, m_pSegmentArray2, (*itc1)->m_ulInputPanelNum)

			// define the scale parameter, such that the max distance of any point within the bbox is AUTOREFINE_MAX_2D_DIM
			scale = AUTOREFINE_MAX_2D_DIM / Mod(m_clsGlobal2DBbox.max_point - m_clsGlobal2DBbox.min_point);

			// scale conductor
			(*itc1)->Scale(scale);

			// scan every panel inside the given conductor group
			for(its1=(*itc1)->m_stlSegments.begin(), i=0; its1!=(*itc1)->m_stlSegments.end(); its1++, i++) {

				// scale segment
				(*its1)->Scale(scale);

				// update conductor bbox
				(*itc1)->m_cls3DBbox += (*its1)->GetCentroid();

				// fill panel array
				m_pSegmentArray1[i] = *its1;
			}

			// now build super hierarchy
			m_iLevel = -1;
			(*itc1)->m_uTopElement.m_pTopSegment = RecurBuild2DSuperHier(0, (*itc1)->m_ulInputPanelNum, &((*itc1)->m_cls3DBbox));
			// if out of memory
			if((*itc1)->m_uTopElement.m_pTopSegment == NULL) {
				return FC_OUT_OF_MEMORY;
			}

			// clear structures of linked list of panels (as we have now a hierarchical structure
			// with a single top panel)
			(*itc1)->m_stlSegments.clear();

			delete [] m_pSegmentArray1;
			delete [] m_pSegmentArray2;

			// compute max area and diameter (used to normalize panel areas)
			// TBC warning: for optimization, should normalize in the eps constant, not in panel areas!

			// keep track of maximum length of any segment in the model
			if( (*itc1)->m_uTopElement.m_pTopSegment->GetLength() > m_dMaxLength) {
				m_dMaxLength = (*itc1)->m_uTopElement.m_pTopSegment->GetLength();
			}

		}
		else {
			ASSERT(false);
		}

	}

	finish = clock();
	m_fDurationSuperH = (float)(finish - start) / CLOCKS_PER_SEC;

	return FC_NORMAL_END;
}


CAutoPanel *CAutoRefine::RecurBuild3DSuperHier(unsigned long firstPanel, unsigned long panelNum, C3DBBox *bbox)
{
	unsigned long i, low_i, hi_i;
	unsigned char maxcoord;
	double midPlane;
	C3DBBox lBbox, rBbox;
	CAutoPanel *leftSubPanel, *rightSubPanel, *currPanel;

	// could assert if there are multiple panels with the same centroid
	// TBC warning: should issue an error message to the user
	ASSERT(panelNum >= 1);

	m_iLevel++;

	// if only one panel left, we have reached the bottom; so return
	if(panelNum == 1) {
		currPanel = m_pPanelArray1[firstPanel];
	}
	else {

		// find division coordinate
		maxcoord = bbox->MaxSide();
		midPlane = (bbox->max_point[maxcoord] + bbox->min_point[maxcoord]) / 2.0;

		// divide panels in two groups according to their position relative to 'midPlane'
		for(i=firstPanel, low_i=firstPanel, hi_i = firstPanel+panelNum-1; i<firstPanel+panelNum; i++) {

			if(m_pPanelArray1[i]->GetCentroid()[maxcoord] <= midPlane) {
				m_pPanelArray2[low_i] = m_pPanelArray1[i];
				low_i++;
			}
			else {
				m_pPanelArray2[hi_i] = m_pPanelArray1[i];
				hi_i--;
			}
		}

		// compute bboxes for left and right groups
		// and copy sorted pointers back to the original array
		for(i=firstPanel; i<low_i; i++) {
			// copy data
			m_pPanelArray1[i] = m_pPanelArray2[i];
			// compute bbox
			lBbox += m_pPanelArray1[i]->GetCentroid();
		}
		for(i=low_i; i<firstPanel+panelNum; i++) {
			// copy data
			m_pPanelArray1[i] = m_pPanelArray2[i];
			// compute bbox
			rBbox += m_pPanelArray1[i]->GetCentroid();
		}

		// recursively build left and right subtrees
		leftSubPanel = RecurBuild3DSuperHier(firstPanel, low_i - firstPanel, &lBbox);
		// if out of memory
		if(leftSubPanel == NULL) {
			return NULL;
		}
		rightSubPanel = RecurBuild3DSuperHier(low_i, firstPanel + panelNum - low_i, &rBbox);
		// if out of memory
		if(rightSubPanel == NULL) {
			return NULL;
		}

		// build new super panel
		//

		// allcate
		// SAFENEW_RET_NULL(TYPE, VAR, MEM)
		SAFENEW_RET_NULL(CAutoPanel, currPanel, g_clsMemUsage.m_ulPanelsMem)
		// and make it
		currPanel->MakeSuperPanel(leftSubPanel, rightSubPanel);

	}

#ifdef DEBUG_DUMP_BASIC
	currPanel->m_iLevel = m_iLevel;
#endif

	m_iLevel--;

	return currPanel;
}

CAutoSegment *CAutoRefine::RecurBuild2DSuperHier(unsigned long firstPanel, unsigned long panelNum, C3DBBox *bbox)
{
	unsigned long i, low_i, hi_i;
	unsigned char maxcoord;
	double midPlane;
	C3DBBox lBbox, rBbox;
	CAutoSegment *leftSubPanel, *rightSubPanel, *currPanel;

	// could assert if there are multiple panels with the same centroid
	// TBC warning: should issue an error message to the user
	ASSERT(panelNum >= 1);

	m_iLevel++;

	// if only one panel left, we have reached the bottom; so return
	if(panelNum == 1) {
		currPanel = m_pSegmentArray1[firstPanel];
	}
	else {

		// find division coordinate
		maxcoord = bbox->MaxSide();
		ASSERT(maxcoord != 2);
		midPlane = (bbox->max_point[maxcoord] + bbox->min_point[maxcoord]) / 2.0;

		// divide panels in two groups according to their position relative to 'midPlane'
		for(i=firstPanel, low_i=firstPanel, hi_i = firstPanel+panelNum-1; i<firstPanel+panelNum; i++) {

			if(m_pSegmentArray1[i]->GetCentroid()[maxcoord] <= midPlane) {
				m_pSegmentArray2[low_i] = m_pSegmentArray1[i];
				low_i++;
			}
			else {
				m_pSegmentArray2[hi_i] = m_pSegmentArray1[i];
				hi_i--;
			}
		}

		// compute bboxes for left and right groups
		// and copy sorted pointers back to the original array
		for(i=firstPanel; i<low_i; i++) {
			// copy data
			m_pSegmentArray1[i] = m_pSegmentArray2[i];
			// compute bbox
			lBbox += m_pSegmentArray1[i]->GetCentroid();
		}
		for(i=low_i; i<firstPanel+panelNum; i++) {
			// copy data
			m_pSegmentArray1[i] = m_pSegmentArray2[i];
			// compute bbox
			rBbox += m_pSegmentArray1[i]->GetCentroid();
		}

		// recursively build left and right subtrees
		leftSubPanel = RecurBuild2DSuperHier(firstPanel, low_i - firstPanel, &lBbox);
		// if out of memory
		if(leftSubPanel == NULL) {
			return NULL;
		}
		rightSubPanel = RecurBuild2DSuperHier(low_i, firstPanel + panelNum - low_i, &rBbox);
		// if out of memory
		if(rightSubPanel == NULL) {
			return NULL;
		}

		// build new super panel
		//

		// allcate
		// SAFENEW_RET_NULL(TYPE, VAR, MEM)
		SAFENEW_RET_NULL(CAutoSegment, currPanel, g_clsMemUsage.m_ulPanelsMem)
		// and make it
		currPanel->MakeSuperSegment(leftSubPanel, rightSubPanel);

	}

#ifdef DEBUG_DUMP_BASIC
	currPanel->m_iLevel = m_iLevel;
#endif

	m_iLevel--;

	return currPanel;
}

/*

// memory info workaround
//

void wxGetMemoryInfo(long *memTotal, long *memUsed, long *memFree)
{
   // zero everything to start with
   *(memTotal) = *(memUsed) = *(memFree) = 0;
#if defined(__WINDOWS__)
    MEMORYSTATUS memStatus;
    memStatus.dwLength = sizeof(MEMORYSTATUS);
    GlobalMemoryStatus(&memStatus);
    *(memTotal) = memStatus.dwTotalPhys;
    *(memUsed) = memStatus.dwMemoryLoad;
    *(memFree) = memStatus.dwAvailPhys;
#elif defined(__UNIX__) || defined(__UNIX_LIKE__)
#  if !defined(__BSD__) && !defined(__DARWIN__)

   // this should be fine on most unix-like platforms
   // note: it doesn't report memory dedicated to the kernel
   //       but is otherwise, reasonably accurate.
   //       e.g. 256MiB installed, but reports 248MiB
    long page_size = sysconf(_SC_PAGESIZE);
    *(memTotal) = sysconf(_SC_PHYS_PAGES)*page_size;
    *(memFree)  = sysconf(_SC_AVPHYS_PAGES)*page_size;
    *(memUsed)  = *(memTotal) - *(memFree);

#   elif defined(__BSD__) || defined(__DARWIN__)

    long total_pages, free_pages, page_size;

    size_t size = sizeof(long);
    if(sysctlbyname("hw.pagesize", &page_size, &size, 0, 0) != 0 ||
       size != sizeof(long)) return;

    size = sizeof(long);
    if(sysctlbyname("vm.stats.vm.v_page_count", &total_pages, &size, 0, 0) != 0 ||
       size != sizeof(long)) return;

    // this may not work under FreeBSD?
    size = sizeof(long);
    if(sysctlbyname("vm.stats.vm.v_free_count", &free_pages, &size, 0, 0) != 0 ||
       size != sizeof(long)) return;

    *(memTotal) = page_size * total_pages;
    *(memFree) = page_size * free_pages;
    *(memUsed) = page_size * (total_pages - free_pages);
#   endif
#else
    // not implemented but defined for
    // portability - values were nulled above
#endif
}

*/


