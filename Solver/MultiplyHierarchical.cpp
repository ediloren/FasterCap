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


// MultiplyHierarchical.cpp : hierarchical multiplication class
//
// Herarchical multiplication algorithm
// (uses structures built by Autorefine)
//
// Enrico Di Lorenzo, 2005/10/04

#include "../stdafx.h"

#include <time.h>
// for multi-thread timing
#include "omp.h"

#include <cmath>

#include "MultiplyHierarchical.h"

// scale factor for the last term in 2D capacitance extraction
// (sum of charges to be enforced = 0), to improve the conditioning of the matrix
// This value has been empirically found
#define MULTIPLYHIER_SCALE_FACTOR   1.0

CMultHier::CMultHier()
{

}

CMultHier::~CMultHier()
{
	// garbage collection
	DeallocateMemory();
}

int CMultHier::AllocateMemory()
{
	return FC_NORMAL_END;
}

void CMultHier::DeallocateMemory()
{
	// garbage collection not strictly needed (setting to NULL and zeroing memory counter)
	// but in case this is moved out of destructor, it saves time and errors
	g_clsMemUsage.m_ulHierMem = 0;
}

// Compute the charge of all nodes
// Unrolling the recursion, for speed
void CMultHier::ComputePanelCharges_fast()
{
	int i;
	StlAutoCondDeque::iterator itc1;

	// init leaf panel counter
	m_dIndex = 0;

	// zeroes recursion vector
	for(i=0; i<MULTHIER_MAX_RECURS_DEPTH; i++) {
		m_clsRecursVec[i] = NULL;
	}

	// scan every conductor
	for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
		// and compute the charge of its panels
		i = 0;
		m_clsRecursVec[i] = (*itc1)->m_uTopElement.m_pTopPanel;

		do {
			if(m_clsRecursVec[i]->IsLeaf()) {
				// if leaf panel, simply copy charge value from input array
				m_clsRecursVec[i]->m_dCharge = m_clsChargeVect[m_dIndex];
				// pre-calculate leaf panel potential based on self-coefficient of potential
				// (for next ComputePanelPotentials() routine), instead of zeroing it and calculating it later,
				// since the array of self-potential is separated from the mutual coefficients array
				m_clsRecursVec[i]->m_dPotential = m_clsSelfPotCoeff[m_dIndex] * m_clsRecursVec[i]->m_dCharge;
				m_dIndex++;
				// and go up in the stack
				i--;
			}
			else {
				// if already visited right branch, we have finished here
				// and charge is the sum of left and right sub-tree charges
				if(m_clsRecursVec[i+1] == m_clsRecursVec[i]->m_pRight) {
					m_clsRecursVec[i]->m_dCharge = m_clsRecursVec[i]->m_pLeft->m_dCharge + m_clsRecursVec[i]->m_pRight->m_dCharge;
					ASSERT(fabs(m_clsRecursVec[i]->m_dCharge) < 1E20);
					// clear potentials (for next ComputePanelPotentials() routine)
					m_clsRecursVec[i]->m_dPotential = 0.0;
					i--;
				}
				else if(m_clsRecursVec[i+1] == m_clsRecursVec[i]->m_pLeft) {
					// finished with left branch, must go visit the right one
					m_clsRecursVec[i+1] = m_clsRecursVec[i]->m_pRight;
					i++;
				}
				else {
					// if none of the two, must go visit the left branch
					m_clsRecursVec[i+1] = m_clsRecursVec[i]->m_pLeft;
					i++;
				}
			}
		}
		while (i>=0);
	}
}

// Compute the potential of each panel due to the directly
// interacting panels (bottom-up or top-down is the same, no distribution
// or gathering between different levels is performed)
//
// Unrolling the recursion, for speed
//
// This is the most costly operation in terms of time spent, since
// it must cycle through all links
//
// Further accelerate, using links stored in array instead
// of linked list, mimicking what is done by direct inversion of the
// matrix, that in spite of n^3 complexity, is very fast for small n
int CMultHier::ComputePanelPotentials_2fast()
{
	int ret;
	unsigned long linkIndex, chunk, block, nodeIndex, nodeBlockEnd, linksPerBlock;
	// 'i' must have signed integral type due to MS OpenMP limitation (does not accept unsigned)
	long i;
	double **localPotCoeffLinks, **localBottomPotCoeffLinks;
	CAutoElement ***localPanelPtrLinks, ***localBottomPnPtrLinks;

	// init pointers to link arrays
	localPotCoeffLinks = m_dPotCoeffLinks[m_ucInteractionLevel];
	localPanelPtrLinks = m_pdPanelPtrLinks[m_ucInteractionLevel];
	// needed only for hierarchical precond upper level
	localBottomPotCoeffLinks = m_dPotCoeffLinks[AUTOREFINE_HIER_PRE_0_LEVEL];
	localBottomPnPtrLinks = m_pdPanelPtrLinks[AUTOREFINE_HIER_PRE_0_LEVEL];

	m_ulCurrBlock = 0;
	// only if we went out-of-core, pre-load first set of chunks
	if(m_ulBlocksNum > 1) {
		ret = LoadLinks(m_ulCurrBlock);
		if(ret != FC_NORMAL_END) {
			return ret;
		}
	}

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
			ret = LoadLinks(block);
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

        // now calculate potentials based on charges, performing multiplication and accumulation

#pragma omp parallel for
        for(i=nodeIndex; i<(long)nodeBlockEnd; i++) {
            unsigned long localLinkIndex, localChunk, localPosInChunk;
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
                        // node potentials have already been zeroed in ComputePanelCharges_fast()
                        m_pNodes[i]->m_dPotential += (localPanelPtrLinks[localChunk][localPosInChunk])->m_dCharge * localPotCoeffLinks[localChunk][localPosInChunk];
                        ASSERT(fabs(m_pNodes[i]->m_dPotential) < 1E20);
                    }
                }
            }
        }
        // and move to next nodeIndex position
        nodeIndex = nodeBlockEnd - 1;
    }

	return FC_NORMAL_END;
}

// Compute the potential of each leaf panels from potential
// on all panels, summing it down
// Unrolling the recursion, for speed
void CMultHier::ComputeLeafPotentials_fast()
{
	int i;
	StlAutoCondDeque::iterator itc1;

	// init leaf panel counter
	m_dIndex = 0;

	// zeroes recursion vector
	for(i=0; i<MULTHIER_MAX_RECURS_DEPTH; i++) {
		m_clsRecursVec[i] = NULL;
	}

	// scan every conductor
	for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
		// and compute the charge of its panels
		i = 0;
		m_clsRecursVec[i] = (*itc1)->m_uTopElement.m_pTopPanel;

		do {
			if(m_clsRecursVec[i]->IsLeaf()) {
				// if leaf, store resulting potential in the array
				// to be returned to the caller as mult result
				m_clsPotVect[m_dIndex] = m_clsRecursVec[i]->m_dPotential;
				m_dIndex++;
				// and go up in the stack
				i--;
			}
			else {
				// if already visited right branch, we have finished here
				// and charge is the sum of left and right sub-tree charges
				if(m_clsRecursVec[i+1] == m_clsRecursVec[i]->m_pRight) {
					i--;
				}
				else if(m_clsRecursVec[i+1] == m_clsRecursVec[i]->m_pLeft) {
					// finished with left branch, must go visit the right one
					m_clsRecursVec[i+1] = m_clsRecursVec[i]->m_pRight;
					i++;
				}
				else {
					// if none of the two
					// add panel potential to children potentials
					m_clsRecursVec[i]->m_pLeft->m_dPotential += m_clsRecursVec[i]->m_dPotential;
					m_clsRecursVec[i]->m_pRight->m_dPotential += m_clsRecursVec[i]->m_dPotential;
					// and go visit the left branch
					m_clsRecursVec[i+1] = m_clsRecursVec[i]->m_pLeft;
					i++;
				}
			}
		}
		while (i>=0);
	}
}

int CMultHier::MultiplyMatByVec_fast(CLin_Vector *v, CLin_Vector *q)
{
	int ret;
	unsigned long size, halfsize, i, j;
	double lastrow;

	// allocate vectors
	if( m_clsGlobalVars.m_ucHasCmplxPerm == AUTOREFINE_REAL_PERM ) {
		// store locally the reference to
		// the charge and the potential vectors
		m_clsChargeVect = *q;
		m_clsPotVect = *v;

		// gather charges, calculate potentials and
		// sum contributions into the leaves
		ComputePanelCharges_fast();
		//	ComputePanelPotentials_fast();
		ret = ComputePanelPotentials_2fast();
		ComputeLeafPotentials_fast();

		if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
			// if 2D, must correct for the 'k' integration constant;
			// we use the work-around from A. R. Djordjevic, R. F. Harrington, T. K. Sarkar,
			// "Evaluation of quasi-static matrix parameters for multiconductor transmission
			// lines using Galerkin's method", IEEE Transactions on Microwave Theory and Techniques,
			// Vol. 42, No. 7, 1994

			size = (*q).size();
			ASSERT(size == (*v).size());

			// subtract last row from all the previous conductor rows, except the last one
			lastrow = (*v)[size-1];
			for(i=m_ulFirstCondElemIndex; i<size-1; i++) {
				(*v)[i] -= lastrow;
			}
			// and 'create' last row as sum of all charges
			// scaling by MULTIPLYHIER_SCALE_FACTOR and in case by TWO_PI_TIMES_E0
			// is to have a better matrix conditioning
			// having these values of more or less the same dimension of the other ones
			// (anyway the sum should be zero)
			lastrow = 0.0;
			for(i=0; i<size; i++) {
				lastrow += (*q)[i];
			}
			(*v)[size-1] = MULTIPLYHIER_SCALE_FACTOR * lastrow;
		}
	}
	else {

		// the complex matrix is represented in real format as the block matrix [R -C; C R]
		// so since the multiplication of R by a vector is accelerated (R is hierarchically represented),
		// we can perform the multiplication of the complex matrix in two steps
		// remark: also the vector represents complex numbers, where the first half of the vector contains
		// the real parts and the second one the complex parts

		size = (*q).size();
		ASSERT(size == (*v).size());
		halfsize = size / 2;
		ASSERT(halfsize * 2 == size);

		//
		// first step, multiply the first n rows by '*q': vRe = [R -C] * (*q)
		//

		// start mutliplying the R block by the first half of the 'q' vector

		// store locally the reference to
		// the charge and the potential vectors
		m_clsChargeVect = CLin_Range(*q, 0, halfsize-1);
		m_clsPotVect = CLin_Range(*v, 0, halfsize-1);

		// gather charges, calculate potentials and
		// sum contributions into the leaves
		ComputePanelCharges_fast();
		//	ComputePanelPotentials_fast();
		ret = ComputePanelPotentials_2fast();
		ComputeLeafPotentials_fast();

		// then add to the potential vector the contribution of the multiplication of the C block
		// by the second half of the 'q' vector. This is easier since the C block is diagonal
		// (btw many elements on the diagonal will be zeros, i.e. the ones not corresponding to a diel. interface)

		for(i=0, j=halfsize; i<halfsize; i++, j++) {
			(*v)[i] -= m_clsImgSelfPotCoeff[i] * (*q)[j];
		}

		if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
			// if 2D, must correct for the 'k' integration constant;
			// we use the work-around from A. R. Djordjevic, R. F. Harrington, T. K. Sarkar,
			// "Evaluation of quasi-static matrix parameters for multiconductor transmission
			// lines using Galerkin's method", IEEE Transactions on Microwave Theory and Techniques,
			// Vol. 42, No. 7, 1994

			// this is the correction for the real part

			// subtract last row from all the previous conductor rows, except the last one
			lastrow = (*v)[halfsize-1];
			for(i=m_ulFirstCondElemIndex; i<halfsize-1; i++) {
				(*v)[i] -= lastrow;
			}
			// and 'create' last row as sum of all charges
			// scaling by MULTIPLYHIER_SCALE_FACTOR and in case by TWO_PI_TIMES_E0
			// is to have a better matrix conditioning
			// having these values of more or less the same dimension of the other ones
			// (anyway the sum should be zero)

			// this is the sum of the real parts

			lastrow = 0.0;
			for(i=0; i<halfsize; i++) {
				lastrow += (*q)[i];
			}
			(*v)[halfsize-1] = MULTIPLYHIER_SCALE_FACTOR * lastrow;
		}

		//
		// second step, multiply the remaining n rows by '*q': vIm = [C R] * (*q)
		//

		// start mutliplying the R block by the second half of the 'q' vector

		// store locally the reference to
		// the charge and the potential vectors
		m_clsChargeVect = CLin_Range(*q, halfsize, size);
		m_clsPotVect = CLin_Range(*v, halfsize, size);

		// gather charges, calculate potentials and
		// sum contributions into the leaves
		ComputePanelCharges_fast();
		//	ComputePanelPotentials_fast();
		ret = ComputePanelPotentials_2fast();
		ComputeLeafPotentials_fast();

		// then add to the potential vector the contribution of the multiplication of the C block
		// by the first half of the 'q' vector. This is easier since the C block is diagonal
		// (btw many elements on the diagonal will be zeros, i.e. the ones not corresponding to a diel. interface)

		for(i=halfsize, j=0; i<size; i++, j++) {
			(*v)[i] += m_clsImgSelfPotCoeff[j] * (*q)[j];
		}

		if(g_ucSolverType == SOLVERGLOBAL_2DSOLVER) {
			// if 2D, must correct for the 'k' integration constant;
			// we use the work-around from A. R. Djordjevic, R. F. Harrington, T. K. Sarkar,
			// "Evaluation of quasi-static matrix parameters for multiconductor transmission
			// lines using Galerkin's method", IEEE Transactions on Microwave Theory and Techniques,
			// Vol. 42, No. 7, 1994

			// this is the correction for the imaginary part

			// subtract last row from all the previous conductor rows, except the last one
			lastrow = (*v)[size-1];
			for(i=m_ulFirstCondElemIndex + halfsize; i<size-1; i++) {
				(*v)[i] -= lastrow;
			}
			// and 'create' last row as sum of all charges
			// scaling by MULTIPLYHIER_SCALE_FACTOR and in case by TWO_PI_TIMES_E0
			// is to have a better matrix conditioning
			// having these values of more or less the same dimension of the other ones
			// (anyway the sum should be zero)

			// this is the sum of the imaginary parts

			lastrow = 0.0;
			for(i=halfsize; i<size; i++) {
				lastrow += (*q)[i];
			}
			(*v)[size-1] = MULTIPLYHIER_SCALE_FACTOR * lastrow;
		}

	}

	return ret;
}

// Recursively copy the charges from panels to vector
void CMultHier::CopyCharges(CAutoElement* panel)
{

	// if panel is not leaf, go on
	if (panel->IsLeaf() == false) {
		CopyCharges(panel->m_pLeft);
		CopyCharges(panel->m_pRight);
	}
	else {
		// if leaf, store charge in the array
		// to be returned to the caller
		m_clsChargeVect[m_dIndex] = panel->m_dCharge;
		m_dIndex++;
	}
}

// Copy the charges from panels to vector
void CMultHier::CopyPanelCharges()
{
	StlAutoCondDeque::iterator itc1;

	// init leaf panel counter
	m_dIndex = 0;

	// scan every conductor
	for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
		// and compute the total potential of its leaf panels
		CopyCharges((*itc1)->m_uTopElement.m_pTopPanel);
	}
}

void CMultHier::CopyChargesToVec(CLin_Vector *q)
{

	// store locally the reference to
	// the charge vector
	m_clsChargeVect = *q;

	// copy charges from leaf panels to charge vector
	CopyPanelCharges();
}

// Recursively copy the charges from vector to panels
void CMultHier::CopyVec(CAutoElement* panel)
{

	// if panel is not leaf, go on
	if (panel->IsLeaf() == false) {
		CopyVec(panel->m_pLeft);
		CopyVec(panel->m_pRight);
	}
	else {
		// if leaf, store charge in the array
		// to be returned to the caller
		panel->m_dCharge = m_clsChargeVect[m_dIndex];
		m_dIndex++;
	}
}

// Copy the charges from vector to panels
void CMultHier::CopyVecCharges()
{
	StlAutoCondDeque::iterator itc1;

	// init leaf panel counter
	m_dIndex = 0;

	// scan every conductor
	for(itc1 = m_stlConductors.begin(); itc1 != m_stlConductors.end(); itc1++) {
		// and compute the total potential of its leaf panels
		CopyVec((*itc1)->m_uTopElement.m_pTopPanel);
	}
}

void CMultHier::CopyVecToCharges(CLin_Vector *q)
{

	// store locally the reference to
	// the charge vector
	m_clsChargeVect = *q;

	// copy charge vector to leaf panels
	CopyVecCharges();
}
