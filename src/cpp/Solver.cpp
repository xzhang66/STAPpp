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

int MIN(int I, int J)
{
	if (I < J) 
		return I;
	else 
		return J;
};

int MAX(int I, int J)
{
	if (I > J) 
		return I;
	else 
		return J;
};

Solver::Solver(Domain* FEMData) : FEMData(FEMData) {};

// LDLT facterization
void LDLTSolver::LDLT()
{
	double* K = FEMData->GetStiffnessMatrix();

	unsigned int* Address = FEMData->GetDiagonalAddress();

	unsigned int N = FEMData->GetNEQ();

	for (int j = 0; j < N; j++)      // Loop for column
	{
		double* Columnj = &K[Address[j] - 1];

        // Total number of non-zero elements in column j
		int ColumnNumberj = Address[j + 1] - Address[j];

        // Row number of the first non-zero element in column j
		int Heightj = j - ColumnNumberj + 1;
        
		for (int i = Heightj; i <= j; i++)  // Loop for all nonzero elements in column j
		{
            // Total number of nonzero elements in column i
			int ColumnNumberi = Address[i + 1] - Address[i];
			double* Columni = &K[Address[i] - 1];
			int CurPostion = Address[j] + j - i - 1;

			double C = 0;
            
            // Row number of the first nonzero element in column i
			int Heighti = i - ColumnNumberi + 1;
            
			int Height = MAX(Heighti, Heightj);
			for (int M = Height; M < i; M++)
			{
				int AddressI = Address[i] + i - M - 1;
				int AddressJ = Address[j] + j - M - 1;
				C += K[AddressI] * K[AddressJ] * K[Address[M] - 1];
			}

			if (i == j)
			{
				K[CurPostion] = K[CurPostion] - C;
				if (fabs(K[CurPostion]) < FLT_MIN)
				{
					cout << "*** Error *** Stiffness matrix is not positive definite !" << endl
                         << "    Euqation no = " << i + 1 << endl;

                    exit(4);
				}
			}
			else 
				K[CurPostion] = (K[CurPostion] - C) / K[Address[i] - 1];
		}
	}
};

// Solve displacement by back substitution
void LDLTSolver::ComputeDisplacement()
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

void LDLTSolver::Solve()
{ 
	Outputter* Output = Outputter::Instance();

	LDLT();

	for (int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
	{
		FEMData->AssembleForce(lcase + 1);

		ComputeDisplacement();
	
#ifdef _DEBUG_
		Outputter* Output = Outputter::Instance();
		Output->PrintDisplacement(lcase);
#endif

		Output->OutputNodalDisplacement(lcase);
	}

	return; 
};