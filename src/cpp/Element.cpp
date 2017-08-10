/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "Element.h"

//  Calculate the column height, used with the skyline storage scheme
void Element::ColumnHeight(unsigned int* ColumnHeight)
{
//	Obtain the location matrix: the global equation number that corresponding to each DOF of the element
	vector<int> LocationMatrix;
	for (int N = 0; N < NEN; N++)
		for (int D = 0; D < Node::NDF; D++)
			LocationMatrix.push_back(nodes[N]->bcode[D]);

//	Calculate the column height contributed by this element
	for (int i = 0; i < LocationMatrix.size(); i++)
	{
		int Ilocation = LocationMatrix[i];
		if (!Ilocation)
			continue;

		for (int j = i + 1; j < LocationMatrix.size(); j++)
		{
			int Jlocation = LocationMatrix[j];
			if (!Jlocation)
				continue;

			// Upper triangular part (row number <= column number)
			if (Ilocation > Jlocation)
			{
				int temp = Ilocation;
				Ilocation = Jlocation;
				Jlocation = temp;
			}

			int Height = Jlocation - Ilocation;
			if (ColumnHeight[Jlocation] < Height) ColumnHeight[Jlocation] = Height;
		}
	}
}

//	Assemble global stiffness matrix (this should be same for all element type ? to be moved to approriate posion)
void Element::assembly(double* Matrix, double* StiffnessMatrix, unsigned int* DiagonalAddress)
{
//	Calculate element stiffness matrix
	ElementStiffness(Matrix);

//	Obtain the location matrix
	vector<int> LocationMatrix;
	for (int i=0; i<NEN; i++)
		for (int j=0; j<Node::NDF; j++)
			LocationMatrix.push_back(nodes[i]->bcode[j]);

//	Assemble global stiffness matrix
	for (int j = 0; j < LocationMatrix.size(); j++)
	{
		int Lj = LocationMatrix[j];	// Global equation number corresponding to jth DOF of the element
		if (!Lj) 
			continue;

//		Address of diagonal element of column j in the one dimensional element stiffness matrix
		int DiagjElement = (j+1)*j/2 + 1;

//		Address of diagonal element of column j in the banded global stiffness matrix
		int DiagjGlobal = DiagonalAddress[Lj - 1];

		for (int i = 0; i <= j; i++)
		{
			int Li = LocationMatrix[i];	// Global equation number corresponding to ith DOF of the element
			if (!Li) 
				continue;

			StiffnessMatrix[DiagjGlobal + Lj - Li - 1] += Matrix[DiagjElement + j - i - 1];
		}
	}

	return;
}
