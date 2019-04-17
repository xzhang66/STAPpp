/*****************************************************************************/
//Auther:Fan Yang
//Date: 14th April 2019
/*****************************************************************************/

#include "4Q.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Dense> 
using namespace Eigen;
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
	C4QMaterial* material_ = dynamic_cast<C4QMaterial*>(ElementMaterial_);	// Pointer to material of the element

	double E = material_->E;
	double v = material_->Thickness;
	MatrixXd D(3, 3);
	D(0, 0) = 1;
	D(0, 1) = v;
	D(0, 2) = 0;
	D(1, 0) = v;
	D(1, 1) = 1;
	D(1, 2) = 0;
	D(2, 0) = 0;
	D(2, 1) = 0;
	D(2, 2) = (1 - v) / 2;
	D = D * E / (1 - pow(v,2));
	double x1 = nodes_[0]->XYZ[0];
	double y1 = nodes_[0]->XYZ[1];
	double x2 = nodes_[1]->XYZ[0];
	double y2 = nodes_[1]->XYZ[1];
	double x3 = nodes_[2]->XYZ[0];
	double y3 = nodes_[2]->XYZ[1];
	double x4 = nodes_[3]->XYZ[0];
	double y4 = nodes_[3]->XYZ[1];
	// use 4 point gaussian intergration, while 0 denotes -0.577,1 represents 0.577
	//the first one is x the second one is y,eg,J00 denotes the determinant of Jacobian at -0.577,-0.577
	double psi[2];
	double eta[2];
	psi[0] = -0.57735;
	psi[1] = 0.57735;
	eta[0] = -0.57735;
	eta[1] = 0.57735;
	MatrixXd C(4, 2);
	C(0, 0) = x1;
	C(0, 1) = y1;
	C(1, 0) = x2;
	C(1, 1) = y2;
	C(2, 0) = x3;
	C(2, 1) = y3;
	C(3, 0) = x4;
	C(3, 1) = y4;
	MatrixXd GN(2, 4);
	MatrixXd J(2, 2);
	MatrixXd Bpart(2, 4);
	MatrixXd B(3, 8);
	MatrixXd K=MatrixXd::Zero(8, 8);
	
	for (int i = 0; i < 2; i++)//psi
	{
		for (int j = 0; j < 2; j++)//eta
		{
			GN(0, 0) = eta[j] - 1;
			GN(0, 1) = 1 - eta[j];
			GN(0, 2) = 1 + eta[j];
			GN(0, 3) = -eta[j] - 1;
			GN(1, 0) = psi[i] - 1;
			GN(1, 1) = -psi[i] - 1;
			GN(1, 2) = 1 + psi[i];
			GN(1, 3) = 1 - psi[i];
			GN = GN / 4;
			J = GN * C;
			Bpart = J.inverse()*GN;
			B(0, 0) = Bpart(0, 0);
			B(0, 1) = 0;
			B(0, 2) = Bpart(0, 1);
			B(0, 3) = 0;
			B(0, 4) = Bpart(0, 2);
			B(0, 5) = 0;
			B(0, 6) = Bpart(0, 3);
			B(0, 7) = 0;

			B(1, 0) = 0;
			B(1, 1) = Bpart(1, 0);
			B(1, 2) = 0;
			B(1, 3) = Bpart(1, 1);
			B(1, 4) = 0;
			B(1, 5) = Bpart(1, 2);
			B(1, 6) = 0;
			B(1, 7) = Bpart(1, 3);

			B(2, 0) = Bpart(1, 0);
			B(2, 1) = Bpart(0, 0);
			B(2, 2) = Bpart(1, 1);
			B(2, 3) = Bpart(0, 1);
			B(2, 4) = Bpart(1, 2);
			B(2, 5) = Bpart(0, 2);
			B(2, 6) = Bpart(1, 3);
			B(2, 7) = Bpart(0, 3);
			K = K + B.transpose()*D*B*J.determinant();
		}
	}
	//double J00 = 0.197*x1*y2 - 0.197*x2*y2 - 0.197*x1*y4 + 0.0528*x2*y3 - 0.0528*x3*y2 + 0.144*x2*y4 + 0.0528*x4*y2 + 0.0528*x3*y4 - 0.0528*x4*y3;
	//double J01 = 0.0528*x1*y2 + 0.144*x1*y3 - 0.0528*x2*y2 - 0.197*x1*y4 + 0.0528*x2*y3 - 0.197*x3*y2 + 0.197*x4*y2 + 0.197*x3*y4 - 0.197*x4*y3;
	//double J10 = 0.197*x1*y2 - 0.144*x1*y3 - 0.197*x2*y2 - 0.0528*x1*y4 + 0.197*x2*y3 - 0.0528*x3*y2 + 0.0528*x4*y2 + 0.0528*x3*y4 - 0.0528*x4*y3;
	//double J11 = 0.0528*x1*y2 - 0.0528*x2*y2 - 0.0528*x1*y4 + 0.197*x2*y3 - 0.197*x3*y2 - 0.144*x2*y4 + 0.197*x4*y2 + 0.197*x3*y4 - 0.197*x4*y3;
	Matrix[0] = K(0, 0);
	Matrix[1] = K(1, 1);
	Matrix[2] = K(1, 0);
	Matrix[3] = K(2, 2);
	Matrix[4] = K(2, 1);
	Matrix[5] = K(2, 0);
	Matrix[6] = K(3, 3);
	Matrix[7] = K(3, 2);
	Matrix[8] = K(3, 1);
	Matrix[9] = K(3, 0);
	Matrix[10] = K(4, 4);
	Matrix[11] = K(4, 3);
	Matrix[12] = K(4, 2);
	Matrix[13] = K(4, 1);
	Matrix[14] = K(4, 0);
	Matrix[15] = K(5, 5);
	Matrix[16] = K(5, 4);
	Matrix[17] = K(5, 3);
	Matrix[18] = K(5, 2);
	Matrix[19] = K(5, 1);
	Matrix[20] = K(5, 0);
	Matrix[21] = K(6, 6);
	Matrix[22] = K(6, 5);
	Matrix[23] = K(6, 4);
	Matrix[24] = K(6, 3);
	Matrix[25] = K(6, 2);
	Matrix[26] = K(6, 1);
	Matrix[27] = K(6, 0);
	Matrix[28] = K(7, 7);
	Matrix[29] = K(7, 6);
	Matrix[30] = K(7, 5);
	Matrix[31] = K(7, 4);
	Matrix[32] = K(7, 3);
	Matrix[33] = K(7, 2);
	Matrix[34] = K(7, 1);
	Matrix[35] = K(7, 0);
	
}

//	Calculate element stress 
void C4Q::ElementStress(double* stress, double* Displacement)
{
	double x1 = nodes_[0]->XYZ[0];
	double y1 = nodes_[0]->XYZ[1];
	double x2 = nodes_[1]->XYZ[0];
	double y2 = nodes_[1]->XYZ[1];
	double x3 = nodes_[2]->XYZ[0];
	double y3 = nodes_[2]->XYZ[1];
	double x4 = nodes_[3]->XYZ[0];
	double y4 = nodes_[3]->XYZ[1];
	MatrixXd C(4, 2);
	C(0, 0) = x1;
	C(0, 1) = y1;
	C(1, 0) = x2;
	C(1, 1) = y2;
	C(2, 0) = x3;
	C(2, 1) = y3;
	C(3, 0) = x4;
	C(3, 1) = y4;

	MatrixXd de(8, 1);
	de(0, 0) = x1;
	de(1, 0) = y1;
	de(2, 0) = x2;
	de(3, 0) = y2;
	de(4, 0) = x3;
	de(5, 0) = y3;
	de(6, 0) = x4;
	de(7, 0) = y4;
	C4QMaterial* material_ = dynamic_cast<C4QMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double E = material_->E;
	double v = material_->Thickness;
	MatrixXd D(3, 3);
	D(0, 0) = 1;
	D(0, 1) = v;
	D(0, 2) = 0;
	D(1, 0) = v;
	D(1, 1) = 1;
	D(1, 2) = 0;
	D(2, 0) = 0;
	D(2, 1) = 0;
	D(2, 2) = (1 - v) / 2;
	D = D * E / (1 - pow(v, 2));

	double psi[2];
	double eta[2];
	psi[0] = -0.57735;
	psi[1] = 0.57735;
	eta[0] = -0.57735;
	eta[1] = 0.57735;
	MatrixXd B;
	MatrixXd sigma(3,3);//use to store the stress matrix
	int count = 01;
	*stress = 0.0;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			count++;
			B = getB(psi[i], eta[j], C);
			sigma = D * B*de;
			stress[count * 3] = sigma(0, 0);
			stress[count * 3 + 1] = sigma(1, 1);
			stress[count * 3 + 2] = sigma(0, 1);
		}
	}
}
MatrixXd C4Q::getB(double psi, double eta, MatrixXd C)
{
	MatrixXd GN(2, 4);
	MatrixXd J(2, 2);
	MatrixXd Bpart(2, 4);
	MatrixXd B(3, 8);
	GN(0, 0) = eta - 1;
	GN(0, 1) = 1 - eta;
	GN(0, 2) = 1 + eta;
	GN(0, 3) = -eta - 1;
	GN(1, 0) = psi - 1;
	GN(1, 1) = -psi - 1;
	GN(1, 2) = 1 + psi;
	GN(1, 3) = 1 - psi;
	GN = GN / 4;
	J = GN * C;
	Bpart = J.inverse()*GN;
	B(0, 0) = Bpart(0, 0);
	B(0, 1) = 0;
	B(0, 2) = Bpart(0, 1);
	B(0, 3) = 0;
	B(0, 4) = Bpart(0, 2);
	B(0, 5) = 0;
	B(0, 6) = Bpart(0, 3);
	B(0, 7) = 0;

	B(1, 0) = 0;
	B(1, 1) = Bpart(1, 0);
	B(1, 2) = 0;
	B(1, 3) = Bpart(1, 1);
	B(1, 4) = 0;
	B(1, 5) = Bpart(1, 2);
	B(1, 6) = 0;
	B(1, 7) = Bpart(1, 3);

	B(2, 0) = Bpart(1, 0);
	B(2, 1) = Bpart(0, 0);
	B(2, 2) = Bpart(1, 1);
	B(2, 3) = Bpart(0, 1);
	B(2, 4) = Bpart(1, 2);
	B(2, 5) = Bpart(0, 2);
	B(2, 6) = Bpart(1, 3);
	B(2, 7) = Bpart(0, 3);
	return B;
}
