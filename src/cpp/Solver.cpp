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
		// Array containing elements in column j
		double* Columnj = &K[Address[j-1] - 1];

		// Address of K_jj in banded matrix (Numbering starting from 0);
		int Address_jj = Address[j-1] - 1;

        // Row number of the first non-zero element in column j (Numbering starting from 1)
		int mj = j - (Address[j] - Address[j-1]) + 1;
        
		for (int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
            // Row number of the first nonzero element in column i (Numbering starting from 1)
			int mi = i - (Address[i] - Address[i-1]) + 1;
            
			int mij = max(mi, mj);

			// Array containing elements in column i
			double* Columni = &K[Address[i-1] - 1];
            
			double C = 0.0;
			for (int r = mij; r <= i-1; r++)	// Loop for max(mi,mj):i-1 (Numbering starting from 1)
			{
				// Address of L_ri and U_rj in banded matrix (Numbering starting from 0)
				int Address_ri = Address[i-1] + i - r - 1;
				int Address_rj = Address[j-1] + j - r - 1;

				C += K[Address_ri] * K[Address_rj];		// C += L_ri * U_rj
			}

			// Address of K[i,j] in banded matrix (Numbering starting from 0)
			int Address_ij = Address[j-1] + j - i - 1;

			K[Address_ij] = K[Address_ij] - C;
		}

		for (int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			// Address of U_rj and D_rr in banded matrix (Numbering starting from 0)
			int Address_rj = Address[j-1] + j - r - 1;	
			int Address_rr = Address[r-1] - 1;

			double Lrj = K[Address_rj] / K[Address_rr];	// L_rj = U_rj / D_rr
			K[Address_jj] = K[Address_jj]  - Lrj * K[Address_rj];	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			K[Address_rj] = Lrj;

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

//	Reduce right-hand-side load vector
	for (int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
		int mi = i - (Address[i] - Address[i-1]) + 1;

		for (int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
		{
			int Address_ji = Address[i-1] + i - j - 1;
			Force[i-1] = Force[i-1] - K[Address_ji] * Force[j-1];
		}
	}

//	Back substitute
	for (int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] = Force[i-1] / K[Address[i-1] - 1];

	for (int j = N; j >= 2; j--)	// Loop for j=N:2
	{
		int mj = j - (Address[j] - Address[j-1]) + 1;

		for (int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
		{
			int Address_ij = Address[j-1] + j - i - 1;
			Force[i-1] -= K[Address_ij] * Force[j-1];
		}
	}
};
