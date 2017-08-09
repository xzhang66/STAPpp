/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "Truss.h"

#include <math.h>
#include <iostream>

using namespace std;

//	Constructor
Bar::Bar()
{
	NEN = 2;	// Each element has 2 nodes
	nodes = new Node*[NEN];

	ElementMaterial = NULL;
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int Bar::SizeOfStiffnessMatrix() { return 21; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void Bar::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	int DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (int i = 0; i < 3; i++)
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];

	int DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

//	Calculate element stiffness matrix

	BarMaterial* material = (BarMaterial*)ElementMaterial;	// Pointer to material of the element

	double k = material->E * material->Area / L / L2;

	Matrix[0] = k*DX2[0];
	Matrix[1] = k*DX2[1];
	Matrix[2] = k*DX2[3];
	Matrix[3] = k*DX2[2];
	Matrix[4] = k*DX2[4];
	Matrix[5] = k*DX2[5];
	Matrix[6] = k*DX2[0];
	Matrix[7] = -k*DX2[5];
	Matrix[8] = -k*DX2[3];
	Matrix[9] = -k*DX2[0];
	Matrix[10] = k*DX2[1];
	Matrix[11] = k*DX2[3];
	Matrix[12] = -k*DX2[4];
	Matrix[13] = -k*DX2[1];
	Matrix[14] = -k*DX2[3];
	Matrix[15] = k*DX2[2];
	Matrix[16] = k*DX2[4];
	Matrix[17] = k*DX2[5];
	Matrix[18] = -k*DX2[2];
	Matrix[19] = -k*DX2[4];
	Matrix[20] = -k*DX2[5];
}
