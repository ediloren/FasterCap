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


// SolveCapacitance.h : capacitance solve class header file
//

#if !defined(AFX_SOLVECAPACITANCE_H__ENRI1975_5487_55D2_9284_04F014CF5600__INCLUDED_)
#define AFX_SOLVECAPACITANCE_H__ENRI1975_5487_55D2_9284_04F014CF5600__INCLUDED_

#include "SolverGlobal.h"

#include "MultiplyHierarchical.h"
// includes for LinALg
#include "LinAlgebra/Vect.h"
#include "LinAlgebra/Mtx.h"

// block preconditioner
// remark: SOLVE_MAX_BLOCK_PRECOND_NUM cannot be bigger than SOLVE_MAX_SUPER_PRECOND_NUM
// ('m_pclsSupPotMtx' matrix is shared between super and block precond and defined of SOLVE_MAX_SUPER_PRECOND_NUM dimension)
// and, due to 'unsigned char' limit of 'm_pucBlockPrecondDim' dimension, cannot exceed 255
#define SOLVE_MAX_BLOCK_PRECOND_NUM		128
#define SOLVE_MIN_BLOCK_PRECOND_NUM		2

// max number of allowed GMRES iterations
// (must be less than the max number fitting in an unsigned int)
#define SOLVE_GMRES_ITER_MAX 1000

// test
#define SOLVE_TEST_ELEM_MAX 1216


class CSolveCap
{

public:
	CSolveCap();
	~CSolveCap();
	int Run(CAutoRefGlobalVars globalVars);
	int SolveMain(CAutoRefGlobalVars globalVars, CLin_Matrix &capRe, CLin_Matrix &capIm);
	void CopyCondNames(StlStringList &stringList);
	void DeallocateMemory(int command, CAutoRefGlobalVars globalVars);
	void DeallocatePrecond();
    void PrintRetError(int retErr);
	    
	CMultHier m_clsMulthier;
	// only for debug
	float m_fDurationSolve;
	float m_fDurationPrecond;

protected:
	int InputFile(CAutoRefGlobalVars *globalVars);
	void AutoSetPrecondType(unsigned long numOfLinks, long condNum, CAutoRefGlobalVars &globalVars);
	void OutputSolveParams(CAutoRefGlobalVars globalVars);
	void OutputMeshParams(CAutoRefGlobalVars globalVars);
	void OutputSolvePrecondType(CAutoRefGlobalVars globalVars);
	void OutputSolveStats(CAutoRefGlobalVars globalVars);
	int OutputCapMtx(const CLin_Matrix &matrixRe, const CLin_Matrix &matrixIm, CAutoRefGlobalVars globalVars);
	int OutputCapMtxToFile(const CLin_Matrix &matrixRe, const CLin_Matrix &matrixIm, CAutoRefGlobalVars globalVars);
    void PrintTime(float solveTime);
	int RefineGeoAndCountLinks();
	int SolveComputeLinks();
	int SolveForCapacitance(CLin_Matrix *cRe, CLin_Matrix *cIm);
	void MakePreconditioner();
	int Solve(CLin_Matrix *cRe, CLin_Matrix *cIm);
	int AllocateMemory();
	void RecurseHierSuperPre(CAutoPanel* panel);
	void IncrementSupPreNum();
	void RecurseComputePrecond(CAutoElement* element);
	void ComputePrecond(CAutoPanel* panel);
	void ComputeSuperPrecond();
	void ComputeBlockPrecond();
	void InvertMatrix(double (*matrix)[SOLVE_MAX_SUPER_PRECOND_NUM], double (*invMatrix)[SOLVE_MAX_SUPER_PRECOND_NUM], unsigned long size);
	int gmresPrecondSFast_test(CLin_Vector *b, CLin_Vector *x, double gmresTol);
	int gmresPrecondSFastAll(CLin_Vector *b, CLin_Vector *x, double gmresTol);
	int gmresPrecondSFastAllUpper(CLin_Vector *b, CLin_Vector *x, double gmresTol, unsigned char precondType);
	int gmresFlexPrecondSFastAll(CLin_Vector *b, CLin_Vector *x, double gmresTol);
	int gmresPrecondSFastAllX0(CLin_Vector *b, CLin_Vector *x, CLin_Vector *x0);
	void ComputePrecondVectFast(CLin_Vector *Pq, CLin_Vector *q, unsigned char precondType);

#ifdef DEBUG_DUMP_POT
    void DebugDumpPotMtxAndIndex();
	void DebugDumpPotMtx();
	void DumpCondPanelPerm(CAutoConductor *cond, CAutoElement *panel, FILE *fp, long condindex);
	void DumpDielPanelPerm(CAutoConductor *cond, CAutoElement *panel, FILE *fp, long condindex);
#endif
#ifdef DEBUG_DUMP_UNCOMP_POT
	void DebugDumpUncomprPotMtx();
	void DebugRecurseCopyPanels(CAutoPanel* panel);
	void DebugRecurseCopyPanels(CAutoSegment* panel);
#endif
#ifdef DEBUG_DUMP_BASIC
#  ifdef DEBUG_DUMP_OTHER
	void DebugDumpPrecondMtx();
#  endif
#endif

	// debug
//	void qsortStub(void *base, size_t num, size_t width, int (__cdecl *compare )(const void *elem1, const void *elem2 ) );


	CAutoRefGlobalVars m_clsGlobalVars;
	unsigned long m_ulPanelNum;
	unsigned long m_ulNumOfLeaves;
	CAutoPanel **m_clsBlockPrecondElements;
	double (*m_pdBlockPrecond)[SOLVE_MAX_BLOCK_PRECOND_NUM];
	unsigned char *m_pucBlockPrecondDim;
	double (*m_pclsTmpCapMtx)[SOLVE_MAX_SUPER_PRECOND_NUM];
	double (*m_pclsSupPrecondMtx)[SOLVE_MAX_SUPER_PRECOND_NUM];
	double (*m_pclsSupPotMtx)[SOLVE_MAX_SUPER_PRECOND_NUM];
	unsigned int m_uiBlockPreNum;
	unsigned long m_ulBlockPreBaseNum;
	bool m_bIsComputingBlock;
	unsigned int *m_puiSupPrecondIndex;
	float *m_pfSupPrecondAreae;
	unsigned long m_ulTmpSupPreIndex;
	unsigned int m_uiSupPreNum, m_uiSupPreLevel;
	unsigned long m_ulHierSupPreNum;
	int m_iLevel;
	unsigned long m_ulSupPreNumOfLeaves;
	CAutoConductor *m_pCurrCond;
	class CSuperPrecondElement
	{
		public:
		CAutoElement *m_pclsHiLevLeaf;
//		long m_lLoLevLeavesNum;
		CAutoConductor *m_pCond;
	};
	CSuperPrecondElement *m_clsSupPrecondElements;
	unsigned int m_uiGmresPrealloc[AUTOPANEL_MAX_NUM_OF_HIERARCHIES];
	unsigned char m_ucAlternatePrecond;
	CLin_Vector *m_pCondCharges;

	// gmres vars
	CLin_Vector *m_pclsGmres_q, *m_pclsGmres_h, *m_pclsGmres_zf;
	CLin_Vector m_clsGmres_g, m_clsGmres_v, m_clsGmres_w, m_clsGmres_z, m_clsGmres_y;
	// c,s vectors
	CLin_Vector m_clsGmres_c, m_clsGmres_s;
	// r vector
	CLin_Vector m_clsGmres_r;
	// Pq vector (preconditioned vector): must declare as visible to whole function, even in case is not used
	// (but does not allocate memory)
	CLin_Vector m_clsGmres_Pq;
	// x0 vector (initial solution vector): must declare as visible to whole function, even in case is not used
	// (but does not allocate memory)
	CLin_Vector m_clsGmres_x0;

	CLin_Vector *m_pclsGmres1_q, *m_pclsGmres1_h;
	CLin_Vector m_clsGmres1_g, m_clsGmres1_v, m_clsGmres1_w, m_clsGmres1_z, m_clsGmres1_y;
	// c,s vectors
	CLin_Vector m_clsGmres1_c, m_clsGmres1_s;
	// r vector
	CLin_Vector m_clsGmres1_r;
	// Pq vector (preconditioned vector): must declare as visible to whole function, even in case is not used
	// (but does not allocate memory)
	CLin_Vector m_clsGmres1_Pq;
	// x0 vector (initial solution vector): must declare as visible to whole function, even in case is not used
	// (but does not allocate memory)
	CLin_Vector m_clsGmres1_x0;

#ifdef DEBUG_DUMP_UNCOMP_POT
	// debug vars
	CAutoPanel **m_pLeafPanels;
	CAutoSegment **m_pLeafSegments;
	unsigned long m_pLeafPanIndex;
#endif //DEBUG_DUMP_UNCOMP_POT
#ifdef DEBUG_DUMP_OTHER
	unsigned long m_ulSkipdPanels;
#endif //DEBUG_DUMP_OTHER

};


#endif //!defined(AFX_SOLVECAPACITANCE_H__ENRI1975_5487_55D2_9284_04F014CF5600__INCLUDED_)
