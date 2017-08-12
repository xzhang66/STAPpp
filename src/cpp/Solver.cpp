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

//	Loop over for all load cases
	for (int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
	{
//		Assemble righ-hand-side vector (force vector)
		FEMData->AssembleForce(lcase + 1);

//		Reduce right-hand-side force vector and back substitute 
		BackSubstitution();
	
#ifdef _DEBUG_
		Outputter* Output = Outputter::Instance();
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

	unsigned int* Address = FEMData->GetDiagonalAddress();

	unsigned int N = FEMData->GetNEQ();

	for (int j = 0; j < N; j++)      // Loop for column
	{
		double* Columnj = &K[Address[j] - 1];

        // Total number of non-zero elements in column j (column height + 1)
		int nEleColumnj = Address[j + 1] - Address[j];

        // Row number of the first non-zero element in column j
		int mj = j - nEleColumnj + 1;
        
//		Loop for all nonzero elements in column j (upper triangular)
		for (int i = mj; i <= j; i++) 
		{
            // Total number of nonzero elements in column i
			int nEleColumni = Address[i + 1] - Address[i];
			double* Columni = &K[Address[i] - 1];
			int Address_ij = Address[j] + j - i - 1;

			double C = 0;
            
            // Row number of the first nonzero element in column i
			int mi = i - nEleColumni + 1;
            
			int mij = max(mi, mj);
			for (int r = mij; r < i; r++)
			{
				int AddressI = Address[i] + i - r - 1;
				int AddressJ = Address[j] + j - r - 1;

				C += K[AddressI] * K[AddressJ] * K[Address[r] - 1];
			}

			if (i == j)
			{
				K[Address_ij] = K[Address_ij] - C;
				if (fabs(K[Address_ij]) < FLT_MIN)
				{
					cout << "*** Error *** Stiffness matrix is not positive definite !" << endl
                         << "    Euqation no = " << i + 1 << endl;

                    exit(4);
				}
			}
			else 
				K[Address_ij] = (K[Address_ij] - C) / K[Address[i] - 1];
		}
	}
};

// Solve displacement by back substitution
void LDLTSolver::BackSubstitution()
{
	double* Force = FEMData->GetForce();        //  Force vector
	double* K = FEMData->GetStiffnessMatrix();  //  Factorized stiffness matrix
	double* U = FEMData->GetDisplacement();     //  Displacement vector

	unsigned int* Address = FEMData->GetDiagonalAddress();
	unsigned int NEQ = FEMData->GetNEQ();

	// L * V = F , V and F share same storage
	for (int i = 0; i < NEQ; i++)
	{
		int Height = Address[i + 1] - Address[i];
		int CurPos = Address[i + 1] - 2;
		for (int M = i - Height + 1; M < i; M++)
		{
			Force[i] -= K[CurPos] * Force[M];
			CurPos--;
		}
	}

	// D * S = V,  S, V and F share same storage
	for (int i = 0; i < NEQ; i++)
	{
		Force[i] = Force[i] / K[Address[i] - 1];
	}

	// LT * U = V
	for (int i = NEQ - 1; i >= 0; i--)
	{
		double C = 0;
		for (int M = NEQ - 1; M > i; M--)
		{
			int Height = Address[M + 1] - Address[M];
			if (M - Height + 1 <= i) 
				C += K[Address[M] - 1 + M - i] * Force[M];
		}
		Force[i] = Force[i] - C;
	}

	for (int i = 0; i < NEQ; i++) 
		U[i] = Force[i];
};
