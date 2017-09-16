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


// AutoRefine.h : autorefine class header file
//

#if !defined(AFX_AUTOREFINE_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
#define AFX_AUTOREFINE_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_

#include "SolverGlobal.h"

#include <string>
#include <map>

#ifdef MS_VS
// for memory state and debug macros (e.g. _ASSERT), when using MS VisualC++
#include <crtdbg.h>
#else
#define _ASSERT wxASSERT
#endif

// for wxFileName
#include <wx/filename.h>

// link with FasterCap main frame
#include "../FasterCapGlobal.h"

// includes for CAutoPanel, CAutoSegment, CAutoConductor class
#include "AutoPanel.h"
#include "AutoSegment.h"
#include "AutoConductor.h"

// includes for potential calculations
#include "Potential.h"

// includes for LinALg
#include "LinAlgebra/Vect.h"

using namespace std;

// Size of the single chunk in which the link array will be split,
// to optimize memory allocations
//
// debug values
//#define AUTOREFINE_LINK_CHUNK_SIZE		100
//#define AUTOREFINE_LINK_CHUNK_SIZE		131072
// actual value
#define AUTOREFINE_LINK_CHUNK_SIZE		1048576


// PotEstimateOpt() return error codes
#define AUTOREFINE_NO_ERROR				0
#define AUTOREFINE_ERROR_AUTOCAP		-1
#define AUTOREFINE_ERROR_ZERO_DIST		-2
#define AUTOREFINE_ERROR_SMALL_DIST		-4
#define AUTOREFINE_ERROR_NAN_OR_INF		-8

// constant used to define the min relative distance at which:
// 1. two super panels can never interact directly
#define AUTOREFINE_MIN_SUPERPANEL_DIST			5.0
// 2. two super panels can interact directly but only
//    if their children can interact directly as well
#define AUTOREFINE_REF_SUPERPANEL_DIST			10.0

// hierarchical levels (only two so far, waiting to possibly extend to n,
// with proper usage of variables, and not static #define)
#define AUTOREFINE_HIER_PRE_0_LEVEL			0
#define AUTOREFINE_HIER_PRE_1_LEVEL			1

// deallocate memory commands
#define AUTOREFINE_DEALLMEM_AT_END			1
#define AUTOREFINE_DEALLMEM_AT_START		2
#define AUTOREFINE_DEALLMEM_CLEAN_CHARGES	4
#define AUTOREFINE_DEALLMEM_ALL				128

// max number of levels in the input file recursion
#define AUTOREFINE_MAX_PARSE_LEVEL				128

// flag meaning that the panel does not belong to a conductor,
// i.e. no reference to the vector of the permittivities stored in the conductor
#define AUTOREFINE_NO_DIEL_INDEX				-1

class CAutoRefine
{

public:
	CAutoRefine();
	~CAutoRefine();
	int AutoRefinePanels(CAutoRefGlobalVars globalVars, unsigned char interactLevel = 0);
	int AutoRefineLinks(CAutoRefGlobalVars globalVars);
	int ReadFastCapFile(CAutoRefGlobalVars *globalVars);
	void Clean(int command, CAutoRefGlobalVars globalVars);
	void OutputFastCapFile(std::string fileinname, std::string suffix, CLin_Vector *condCharges = NULL);
//	void CleanOnlyUpToSupHie();
	void SelfPotential(CAutoPanel *panel, double *potestRe, double *potestIm);
    void SelfPotential(CAutoSegment *panel, double *potestRe, double *potestIm);
	int BuildSuperHierarchy();
	int PotEstimateOpt(CAutoElement *element1, CAutoElement *element2, double &potestim1);
	int PotEstimateOpt(CAutoPanel *panel1, CAutoPanel *panel2, double &potestim1, unsigned char computePrecond = AUTOREFINE_PRECOND_NONE);
	int PotEstimateOpt(CAutoPanel *panel1, CAutoPanel *panel2, double &potestim1, double &potestim2, unsigned char computePrecond = AUTOREFINE_PRECOND_NONE);
    int PotEstimateOpt(CAutoSegment *panel1, CAutoSegment *panel2, double &potestim1, unsigned char computePrecond = AUTOREFINE_PRECOND_NONE);
    int PotEstimateOpt(CAutoSegment *panel1, CAutoSegment *panel2, double &potestim1, double &potestim2, unsigned char computePrecond = AUTOREFINE_PRECOND_NONE);
    void SetCurrentConductor(CAutoConductor *m_pCurrCond);
	inline void SetInteractionLevel(unsigned char interactionLevel)
	{
		m_ucInteractionLevel = interactionLevel;
	}
	inline unsigned char GetInteractionLevel()
	{
		return m_ucInteractionLevel;
	}
	inline unsigned long GetPanelNum()
	{
		return m_ulPanelNum[m_ucInteractionLevel];
	}
	inline unsigned long GetPanelNum(unsigned char interactLevel)
	{
		return m_ulPanelNum[interactLevel];
	}
	inline unsigned long GetNodeNum()
	{
		return m_ulNodeNum[m_ucInteractionLevel];
	}
	inline unsigned long GetNodeNum(unsigned char interactLevel)
	{
		return m_ulNodeNum[interactLevel];
	}
	inline unsigned long GetLinksNum()
	{
		return m_ulLinksNum[m_ucInteractionLevel];
	}
	inline unsigned long GetLinksNum(unsigned char interactLevel)
	{
		return m_ulLinksNum[interactLevel];
	}
	// number of links + number of autolinks (self potentials, one for each panel)
	inline unsigned long GetTotalLinksNum(unsigned char interactLevel)
	{
		return (m_ulLinksNum[interactLevel] + m_ulPanelNum[interactLevel]);
	}
	inline double GetTotalArea()
	{
		return m_dTotalArea;
	}
	inline double GetMaxArea()
	{
		return m_dMaxArea;
	}
 	inline unsigned long GetInputPanelNum()
	{
		return m_ulInputPanelNum;
	}

	CAutoRefGlobalVars m_clsGlobalVars;
	unsigned char m_ucInteractionLevel;
	float *m_fGlobalCharges;
	double m_dMaxMeshEps;
	CLin_Vector m_clsSelfPotCoeff, m_clsImgSelfPotCoeff;
	unsigned char *m_pucDielIndex;

	// statistics vars
	int m_iLevel;
	int m_iMaxLevel;
	unsigned long m_ulNumofpotest;
//	unsigned long m_ulNumofFastPotest;
	unsigned long m_ulInputRawCondPanelNum, m_ulInputRawDielPanelNum;
	unsigned long m_ulPanelNum[AUTOPANEL_MAX_NUM_OF_HIERARCHIES], m_ulInputPanelNum;
	unsigned long m_ulNodeNum[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	unsigned long m_ulLinksNum[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	float m_fDurationDiscretize, m_fDurationRefine, m_fDurationSuperH, m_fDurationReadFile;
	StlAutoCondDeque m_stlConductors;
	long m_lCondNum, m_lDielNum, m_lGroupDielNum;

#ifdef DEBUG_DUMP_BASIC
	int m_iaLinksBtwLevels[64][64];
	int m_iaPotestBtwLevels[64][64];
	long m_laLinkNumStatistics[64];
	long m_laPanelsPerLevel[64];
	long m_laLinksPerLevel[64];
#endif //DEBUG_DUMP_BASIC


protected:
    typedef std::pair<fpos_t, long > StlPosLinenumPair;
    typedef std::map<std::string, StlPosLinenumPair, less<std::string> > StlFilePosMap;

	int ComputeLinks();
	int SaveLinks(bool saveAlsoPot = true);
	int LoadLinks(unsigned long block, bool loadAlsoPot = true);
	void DumpMemoryInfo();
	void CopyCharges(CAutoElement *panel);
	void RecurseIndex(CAutoElement *panel);
	void DeletePanels();
	void DeleteLinkArray(unsigned int level);
	CAutoPanel *RecurBuild3DSuperHier(unsigned long firstPanel, unsigned long panelNum, C3DBBox *bbox);
	CAutoSegment *RecurBuild2DSuperHier(unsigned long firstPanel, unsigned long panelNum, C3DBBox *bbox);
    int CreateFileMap(char *fileinname, FILE *fid, StlFilePosMap *filePosMap);
	int Parse3DInputFile(char *fileinname, FILE *parentFid, StlFilePosMap *parentFilePosMap, CAutoRefGlobalVars *globalVars, bool isdiel = false, C3DVector offset = C3DVector(0,0,0),
						double outpermRe = 1.0, double outpermIm = 0.0, const char *groupname = "",
						double inpermRe = 0.0, double inpermIm = 0.0, C3DVector dielrefpoint = C3DVector(0,0,0));
	int Parse2DInputFile(char *fileinname, FILE *parentFid, StlFilePosMap *parentFilePosMap, CAutoRefGlobalVars *globalVars, bool isdiel = false, C2DVector offset = C2DVector(0,0),
						double outpermRe = 1.0, double outpermIm = 0.0, const char *groupname = "",
						double inpermRe = 0.0, double inpermIm = 0.0, C3DVector dielrefpoint = C3DVector(0,0,0));
	int CreatePanel(C3DVector_float vertex[3], char *nakedname, unsigned char dielIndex,
							  StlAutoCondDeque::iterator *itc, char *fileinname, long linenum, char opType, CAutoRefGlobalVars *globalVars,
							  bool uselocaldiel, C3DVector &localDielrefpoint);
	int CreateSegment(C2DVector_float vertex[2], char *nakedname, unsigned char dielIndex,
							  StlAutoCondDeque::iterator *itc, char *fileinname, long linenum, CAutoRefGlobalVars *globalVars,
							  bool uselocaldiel, C3DVector &localDielrefpoint);
	void OutputPanelTree(char *condname, CAutoPanel *panel, FILE *fout, int dielIndex = AUTOREFINE_NO_DIEL_INDEX);
	void OutputPanelTree(char *condname, CAutoSegment *panel, FILE *fout, int dielIndex = AUTOREFINE_NO_DIEL_INDEX);
	int Discretize(CAutoPanel *panel);
	int DiscretizeSelf(CAutoPanel *panel);
	int DiscretizeSelf(CAutoSegment *panel);
	int DiscretizeMutual(CAutoPanel *panel1, CAutoPanel *panel2, bool selfCond);
	int DiscretizeMutual(CAutoSegment *panel1, CAutoSegment *panel2, bool selfCond);
	int RefineSelf(CAutoPanel *panel);
	int RefineSelf(CAutoSegment *panel);
	int RefineMutual(CAutoPanel *panel1, CAutoPanel *panel2);
    int RefineMutual(CAutoSegment *panel1, CAutoSegment *panel2);
	bool RefineCriteria(double *panel1crit, double *panel2crit, CAutoPanel *panel1, CAutoPanel *panel2, C3DVector *dist,
                        double rdist, double rmax, double eps, double ccoeff = 1.0f);
	bool RefineCriteria(double *panel1crit, double *panel2crit, CAutoSegment *panel1, CAutoSegment *panel2, C2DVector *dist,
                        double rdist, double rmax, double eps, double ccoeff = 1.0f);
//	double RecurseSelfPotential(CAutoPanel *panel);
//	double RecurseMutualPotential(CAutoPanel *panel);
	double MaxSideQ(CAutoPanel *panel);
	int GetConductor(StlAutoCondDeque::iterator *itc, unsigned char *dielIndex, char *name,
						bool isdiel, double outpermRe, double outpermIm, double inpermRe, double inpermIm, C3DVector &dielrefpoint);
	bool FindConductorFromName(char *name, StlAutoCondDeque::iterator &itc);
	void MergeConductors(StlAutoCondDeque::iterator itc, StlAutoCondDeque::iterator itc_old);
    unsigned long PortableGetTempFileName(const wxString &dirName, const wxString &prefix, unsigned long id, wxFileName &filenameobj);


	CPotential m_clsPotential;
	double m_dMaxSide, m_dMaxArea, m_dTotalArea, m_dMaxLength;
	unsigned long m_ulCountPanelNum, m_ulBasePanelNum, m_ulBaseLinksNum, m_ulBlocksNum, m_ulCurrBlock, m_ulCountNodeNum;
	CAutoPanel **m_pPanelArray1, **m_pPanelArray2;
	CAutoSegment **m_pSegmentArray1, **m_pSegmentArray2;
	C3DVector *m_pCentroid;
	bool m_bComputeLinks;
	double **m_dPotCoeffLinks[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	CAutoElement ***m_pdPanelPtrLinks[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	unsigned long m_ulLinkChunkNum[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	unsigned long *m_ulUniqueLinkChunkIDs, *m_ulUniquePotChunkIDs;
	CAutoConductor *m_pCurrentConductor;
	double m_dMaxSigma, m_dMidSigma, m_dMinSigma;
	bool m_bInitCharges, m_bPopulateNodeArray;
	int m_iParseLevel;
	int m_iGroupNum[AUTOREFINE_MAX_PARSE_LEVEL];
	CLin_Vector *m_pLocalCondCharge;
    static unsigned long m_ulTempFileID;
    C2DBBox m_clsGlobal2DBbox;
    CAutoElement **m_pNodes;

#ifdef DEBUG_DUMP_OTHER
public:
	long **m_pBlockMtx;
	double **m_pElemMtx;

	void DebugOutputTree();
	void DebugOutputTree2();
	void DebugOutputHierarchy(CAutoPanel *panel);
protected:
	void DebugOutputPanels2(CAutoPanel *panel, long &blockNo);
	void DebugMarkBlock(CAutoPanel *panel1, CAutoPanel *panel2, long blockNo, double potCoeff);
	void DebugOutputPanels(CAutoPanel *panel, FILE *fout);
	void DebugDepth(CAutoPanel *panel, char depth[]);
	void DebugOutputInteractions(CAutoPanel *panel);
	void DebugOutputInteracting(CAutoPanel *panel1, CAutoPanel *panel2);
	void DebugOutputPanelTree(char *condname, CAutoPanel *panel, FILE *fout);
#endif  //DEBUG_DUMP_OTHER
public:
	void DebugDumpInteractions();
};



#endif //!defined(AFX_AUTOREFINE_H__E89AAF21_5486_11D5_9282_04F014C10000__INCLUDED_)
