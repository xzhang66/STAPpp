/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Bar.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CBar::CBar()
{
	NEN = 2;	// Each element has 2 nodes
	nodes = new CNode*[NEN];
    
    ND = 6;
    LocationMatrix = new unsigned int[ND];

	ElementMaterial = NULL;
}

//	Desconstructor
CBar::~CBar()
{
}

//	Read element data from stream Input
bool CBar::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cout << "*** Error *** Elements must be inputted in order !" << endl 
			 << "   Expected element : " << Ele + 1 << endl
			 << "   Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2;	// Left node number and right node number

	Input >> N1 >> N2 >> MSet;
	ElementMaterial = &((CBarMaterial*)MaterialSets)[MSet - 1];
	nodes[0] = &NodeList[N1 - 1];
	nodes[1] = &NodeList[N2 - 1];

	return true;
}

//	Write element data to stream OutputFile
void CBar::Write(ofstream& OutputFile, unsigned int Ele)
{
	cout << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber 
		 << setw(9) << nodes[1]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
	OutputFile << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber 
			   << setw(9) << nodes[1]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CBar::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//	has 21 elements
unsigned int CBar::SizeOfStiffnessMatrix() { return 21; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBar::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

//	Calculate bar length
	double DX[3];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];

	double DX2[6];	//  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
	DX2[0] = DX[0] * DX[0];
	DX2[1] = DX[1] * DX[1];
	DX2[2] = DX[2] * DX[2];
	DX2[3] = DX[0] * DX[1];
	DX2[4] = DX[1] * DX[2];
	DX2[5] = DX[0] * DX[2];

	double L2 = DX2[0] + DX2[1] + DX2[2];
	double L = sqrt(L2);

//	Calculate element stiffness matrix

	CBarMaterial* material = (CBarMaterial*)ElementMaterial;	// Pointer to material of the element

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

//	Calculate element stress 
void CBar::ElementStress(double* stress, double* Displacement)
{
	CBarMaterial* material = (CBarMaterial*)ElementMaterial;	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
		L2 = L2 + DX[i]*DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material->E / L2;
		S[i+3] = -S[i];
	}
	
	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix[i])
			*stress += S[i] * Displacement[LocationMatrix[i]-1];
	}
}
