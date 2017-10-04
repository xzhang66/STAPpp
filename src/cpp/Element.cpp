/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <algorithm>

#include "Element.h"

//  Calculate the column height, used with the skyline storage scheme
void CElement::CalculateColumnHeight(unsigned int* ColumnHeight)
{
    
//  Generate location matrix
    GenerateLocationMatrix();

//  Look for the row number of the first non-zero element
    int nfirstrow = INT_MAX;
    for (int i = 0; i < ND; i++)
        if (LocationMatrix[i] && LocationMatrix[i] < nfirstrow)
            nfirstrow = LocationMatrix[i];

//	Calculate the column height contributed by this element
	for (int i = 0; i < ND; i++)
	{
		int column = LocationMatrix[i];
		if (!column)
			continue;

		unsigned int Height = column - nfirstrow;
		if (ColumnHeight[column-1] < Height) ColumnHeight[column-1] = Height;
	}
}

//	Assemble the banded global stiffness matrix (skyline storage scheme)
void CElement::assembly(double* Matrix, double* StiffnessMatrix, unsigned int* DiagonalAddress)
{
//	Calculate element stiffness matrix
	ElementStiffness(Matrix);
	
//	Assemble global stiffness matrix
	for (int j = 0; j < ND; j++)
	{
		int Lj = LocationMatrix[j];	// Global equation number corresponding to jth DOF of the element
		if (!Lj) 
			continue;

//		Address of diagonal element of column j in the one dimensional element stiffness matrix
		int DiagjElement = (j+1)*j/2 + 1;

		for (int i = 0; i <= j; i++)
		{
			int Li = LocationMatrix[i];	// Global equation number corresponding to ith DOF of the element
			if (!Li) 
				continue;

            if (Lj>=Li)
                StiffnessMatrix[DiagonalAddress[Lj - 1] + Lj - Li - 1] += Matrix[DiagjElement + j - i - 1];
            else
                StiffnessMatrix[DiagonalAddress[Li - 1] + Li - Lj - 1] += Matrix[DiagjElement + j - i - 1];
		}
	}

	return;
}
