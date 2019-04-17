/*****************************************************************************/
//Auther:Fan Yang
//Date: 14th April 2019
/*****************************************************************************/

#include "4Q.h"

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

//	Constructor
C4Q::C4Q()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];


	ND_ = 8;// each node has 2 dimensions. It's a planar model
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Desconstructor
C4Q::~C4Q()
{
}

//	Read element data from stream Input
bool C4Q::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;

	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "    Expected element : " << Ele + 1 << endl
			<< "    Provided element : " << N << endl;

		return false;
	}

	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// node number of 4 nodes of 4Q element

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CBarMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

//	Write element data to stream
void C4Q::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void C4Q::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 2; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];// store the location matrix into an array
	//let's assume each node has 2 DoFs, then the from i*2 till (i+1)*2-1 means the ith element's location
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4 nodes 4Q element, element stiffness is a 8x8 matrix, whose upper triangular part
//	has 36 elements
unsigned int C4Q::SizeOfStiffnessMatrix() { return 36; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void C4Q::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//	Calculate bar length
	double DX[2];		//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	for (unsigned int i = 0; i < 3; i++)
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

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

	C4QMaterial* material_ = dynamic_cast<C4QMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
	double v = material_->Thickness;

	Matrix[0] = k * DX2[0];
	Matrix[1] = k * DX2[1];
	Matrix[2] = k * DX2[3];
	Matrix[3] = k * DX2[2];
	Matrix[4] = k * DX2[4];
	Matrix[5] = k * DX2[5];
	Matrix[6] = k * DX2[0];
	Matrix[7] = -k * DX2[5];
	Matrix[8] = -k * DX2[3];
	Matrix[9] = -k * DX2[0];
	Matrix[10] = k * DX2[1];
	Matrix[11] = k * DX2[3];
	Matrix[12] = -k * DX2[4];
	Matrix[13] = -k * DX2[1];
	Matrix[14] = -k * DX2[3];
	Matrix[15] = k * DX2[2];
	Matrix[16] = k * DX2[4];
	Matrix[17] = k * DX2[5];
	Matrix[18] = -k * DX2[2];
	Matrix[19] = -k * DX2[4];
	Matrix[20] = -k * DX2[5];
}

//	Calculate element stress 
void C4Q::ElementStress(double* stress, double* Displacement)
{
	CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double DX[3];	//	dx = x2-x1, dy = y2-y1, dz = z2-z1
	double L2 = 0;	//	Square of bar length (L^2)

	for (unsigned int i = 0; i < 3; i++)
	{
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 = L2 + DX[i] * DX[i];
	}

	double S[6];
	for (unsigned int i = 0; i < 3; i++)
	{
		S[i] = -DX[i] * material_->E / L2;
		S[i + 3] = -S[i];
	}

	*stress = 0.0;
	for (unsigned int i = 0; i < 6; i++)
	{
		if (LocationMatrix_[i])
			*stress += S[i] * Displacement[LocationMatrix_[i] - 1];
	}
}
