/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "SkylineMatrix.h"

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver
{
private:
    
    CSkylineMatrix<double>& K;

public:

//!	Constructor
	CLDLTSolver(CSkylineMatrix<double>* K): K(*K) {};

//!	Perform L*D*L(T) factorization of the stiffness matrix
	void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	void BackSubstitution(double* Force); 
};
