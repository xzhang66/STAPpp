#include "Q4.h"

#include <cmath>
#include <iomanip>

using namespace std;

CQ4::CQ4()
{
	NEN_ = 4;
	nodes_ = new CNode*[NEN_];

	ND_ = 8;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

CQ4::~CQ4()
{
}

bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2, N3, N4;

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	return true;
}

void CQ4::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(9) << nodes_[3]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CQ4::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		LocationMatrix_[i++] = nodes_[N]->bcode[0];
		LocationMatrix_[i++] = nodes_[N]->bcode[1];
	}
}

void shapeFunctionDerivatives(double xi, double eta, double dN_dxi[4], double dN_deta[4])
{
	dN_dxi[0] = -0.25 * (1.0 - eta);
	dN_dxi[1] =  0.25 * (1.0 - eta);
	dN_dxi[2] =  0.25 * (1.0 + eta);
	dN_dxi[3] = -0.25 * (1.0 + eta);

	dN_deta[0] = -0.25 * (1.0 - xi);
	dN_deta[1] = -0.25 * (1.0 + xi);
	dN_deta[2] =  0.25 * (1.0 + xi);
	dN_deta[3] =  0.25 * (1.0 - xi);
}

void assembleB(double B[3][8], const double dN_dx[4], const double dN_dy[4])
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 8; j++)
			B[i][j] = 0.0;

	for (int node = 0; node < 4; node++)
	{
		int col = 2 * node;
		B[0][col] = dN_dx[node];
		B[1][col + 1] = dN_dy[node];
		B[2][col] = dN_dy[node];
		B[2][col + 1] = dN_dx[node];
	}
}

void CQ4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);
	double E = material->E;
	double nu = material->Nu;
	double thickness = material->Thickness;

	double D[3][3] = {
		{1.0, nu, 0.0},
		{nu, 1.0, 0.0},
		{0.0, 0.0, (1.0 - nu) / 2.0}
	};
	double factor = E / (1.0 - nu * nu) * thickness;

	double x[4], y[4];
	for (int i = 0; i < 4; i++)
	{
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
	}

	const double g = 1.0 / sqrt(3.0);
	const double gauss[4][2] = {
		{-g, -g}, {g, -g}, {g, g}, {-g, g}
	};

	double K[8][8];
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			K[i][j] = 0.0;

	for (int gp = 0; gp < 4; gp++)
	{
		double xi = gauss[gp][0];
		double eta = gauss[gp][1];
		double dN_dxi[4], dN_deta[4];
		shapeFunctionDerivatives(xi, eta, dN_dxi, dN_deta);

		double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
		for (int i = 0; i < 4; i++)
		{
			J11 += dN_dxi[i] * x[i];
			J12 += dN_dxi[i] * y[i];
			J21 += dN_deta[i] * x[i];
			J22 += dN_deta[i] * y[i];
		}

		double detJ = J11 * J22 - J12 * J21;
		double invJ11 =  J22 / detJ;
		double invJ12 = -J12 / detJ;
		double invJ21 = -J21 / detJ;
		double invJ22 =  J11 / detJ;

		double dN_dx[4], dN_dy[4];
		for (int i = 0; i < 4; i++)
		{
			dN_dx[i] = invJ11 * dN_dxi[i] + invJ12 * dN_deta[i];
			dN_dy[i] = invJ21 * dN_dxi[i] + invJ22 * dN_deta[i];
		}

		double B[3][8];
		assembleB(B, dN_dx, dN_dy);

		double DB[3][8];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 8; j++)
			{
				DB[i][j] = 0.0;
				for (int k = 0; k < 3; k++)
					DB[i][j] += D[i][k] * B[k][j];
			}

		for (int i = 0; i < 8; i++)
			for (int j = i; j < 8; j++)
			{
				double value = 0.0;
				for (int k = 0; k < 3; k++)
					value += B[k][i] * DB[k][j];
				K[i][j] += value * factor * detJ;
			}
	}

	for (int j = 0; j < 8; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = K[i][j];
}

void CQ4::ElementStress(double* stress, double* Displacement)
{
	CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);

	double xi = 0.0;
	double eta = 0.0;
	double dN_dxi[4], dN_deta[4];
	shapeFunctionDerivatives(xi, eta, dN_dxi, dN_deta);

	double x[4], y[4];
	for (int i = 0; i < 4; i++)
	{
		x[i] = nodes_[i]->XYZ[0];
		y[i] = nodes_[i]->XYZ[1];
	}

	double J11 = 0.0, J12 = 0.0, J21 = 0.0, J22 = 0.0;
	for (int i = 0; i < 4; i++)
	{
		J11 += dN_dxi[i] * x[i];
		J12 += dN_dxi[i] * y[i];
		J21 += dN_deta[i] * x[i];
		J22 += dN_deta[i] * y[i];
	}

	double detJ = J11 * J22 - J12 * J21;
	double invJ11 =  J22 / detJ;
	double invJ12 = -J12 / detJ;
	double invJ21 = -J21 / detJ;
	double invJ22 =  J11 / detJ;

	double dN_dx[4], dN_dy[4];
	for (int i = 0; i < 4; i++)
	{
		dN_dx[i] = invJ11 * dN_dxi[i] + invJ12 * dN_deta[i];
		dN_dy[i] = invJ21 * dN_dxi[i] + invJ22 * dN_deta[i];
	}

	double B[3][8];
	assembleB(B, dN_dx, dN_dy);

	double strain[3] = {0.0, 0.0, 0.0};
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 8; j++)
			if (LocationMatrix_[j])
				strain[i] += B[i][j] * Displacement[LocationMatrix_[j] - 1];

	double coeff = material->E / (1.0 - material->Nu * material->Nu);
	stress[0] = coeff * (strain[0] + material->Nu * strain[1]);
	stress[1] = coeff * (material->Nu * strain[0] + strain[1]);
	stress[2] = coeff * ((1.0 - material->Nu) / 2.0) * strain[2];
}