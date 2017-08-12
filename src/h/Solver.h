/*****************************************************************************/
/*  FEM++ : A C++ FEM code sharing the same input data file with STAP90      */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Domain.h"

using namespace std;

//	Base class for a solver
//	New solver should be derived from this base class, and match the storage scheme 
//	of the global stiffness matrix employed in Domain class.
class Solver
{
protected:

	Domain* FEMData;

public:

	Solver(Domain* FEMData);

	virtual void Solve() = 0;
};

//	LDLT solver for skyline storage scheme
class LDLTSolver : public Solver
{
public:

//	Constructor
	LDLTSolver(Domain* FEMData) :Solver(FEMData) {};

//	Perform LDLT factorization of the global stiffness matrix
	void LDLT();

//	Calculate global nodal displacement vector
	void ComputeDisplacement(); 

//	Solve
	virtual void Solve();       // Ω‚FEMŒ Ã‚
};