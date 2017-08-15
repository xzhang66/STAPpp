/*****************************************************************************/
/*  FEM++ : A C++ FEM code sharing the same input data file with STAP90      */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"

#include <cmath>
#include <cfloat>
#include <iostream>

using namespace std;

Solver::Solver(Domain* FEMData) : FEMData(FEMData) {};

void LDLTSolver::Solve()
{ 
	Outputter* Output = Outputter::Instance();

//	Perform L*D*L(T) factorization of stiffness matrix
	LDLT();
	
#ifdef _DEBUG_
	Output->PrintStiffnessMatrix();
#endif

//	Loop over for all load cases
	for (int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
	{
//		Assemble righ-hand-side vector (force vector)
		FEMData->AssembleForce(lcase + 1);

//		Reduce right-hand-side force vector and back substitute 
		BackSubstitution();
	
#ifdef _DEBUG_
		Output->PrintDisplacement(lcase);
#endif

		Output->OutputNodalDisplacement(lcase);
	}

	return; 
};

// LDLT facterization
void LDLTSolver::LDLT()
{
	double* K = FEMData->GetStiffnessMatrix();
	unsigned int* Address = FEMData->GetDiagonalAddress();	// Numbering starting from 1
	unsigned int N = FEMData->GetNEQ();

	for (int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
        // Row number of the first non-zero element in column j (Numbering starting from 1)
		int mj = j - (Address[j] - Address[j-1]) + 1;

		// Address of K_jj and K_mj,j in banded matrix (Numbering starting from 0);
		int Address_jj = Address[j-1] - 1;
		int Address_mjj = Address[j-1] + j - mj - 1;
        
		for (int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
            // Row number of the first nonzero element in column i (Numbering starting from 1)
			int mi = i - (Address[i] - Address[i-1]) + 1;
            
			int mm = max(mi, mj);

			// Address of L_mm,i and U_mm,j in banded matrix (Numbering starting from 0)
			int Address_mmi = Address[i-1] + i - mm - 1;
			int Address_mmj = Address[j-1] + j - mm - 1;
            
			double C = 0.0;
			for (int r = mm; r <= i-1; r++)	// Loop for max(mi,mj):i-1 (Numbering starting from 1)
				C += K[Address_mmi++] * K[Address_mmj++];		// C += L_ri * U_rj

			K[++Address_mjj] = K[Address_mjj] - C;	// U_ij = K_ij - C
		}

		// Address of K_mj,j in banded matrix (Numbering starting from 0);
		Address_mjj = Address[j-1] + j - mj - 1;

		for (int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = K[Address_mjj] / K[Address[r-1] - 1];	// L_rj = U_rj / D_rr
			K[Address_jj] = K[Address_jj]  - Lrj * K[Address_mjj];	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			K[Address_mjj++] = Lrj;

			if (fabs(K[Address_jj]) <= FLT_MIN)
			{
				cout << "*** Error *** Stiffness matrix is not positive definite !" << endl
                     << "    Euqation no = " << r + 1 << endl
					 << "    Pivot = " << K[Address_jj] << endl;

                    exit(4);
				}
		}
	}
};

// Solve displacement by back substitution
void LDLTSolver::BackSubstitution()
{
	double* Force = FEMData->GetForce();        //  Force vector (Numering starting from 1)
	double* K = FEMData->GetStiffnessMatrix();  //  Factorized stiffness matrix

	unsigned int* Address = FEMData->GetDiagonalAddress(); // Numering starting from 1
	unsigned int N = FEMData->GetNEQ();

//	Reduce right-hand-side load vector (LV = R)
	for (int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
		int mi = i - (Address[i] - Address[i-1]) + 1;

		for (int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] = Force[i-1] - K[Address[i-1] + i - j - 1] * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] = Force[i-1] / K[Address[i-1] - 1];	// Vbar = D^(-1) V

	for (int j = N; j >= 2; j--)	// Loop for j=N:2
	{
		int mj = j - (Address[j] - Address[j-1]) + 1;

		for (int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= K[Address[j-1] + j - i - 1] * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
};
