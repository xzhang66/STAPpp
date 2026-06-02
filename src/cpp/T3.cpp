#include "T3.h"

#include <cmath>
#include <iomanip>

using namespace std;

CT3::CT3()
{
	NEN_ = 3;
	nodes_ = new CNode*[NEN_];

	ND_ = 6;
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

CT3::~CT3()
{
}

bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int MSet;
	unsigned int N1, N2, N3;

	Input >> N1 >> N2 >> N3 >> MSet;
	ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];

	return true;
}

void CT3::Write(COutputter& output)
{
	output << setw(11) << nodes_[0]->NodeNumber
		   << setw(9) << nodes_[1]->NodeNumber
		   << setw(9) << nodes_[2]->NodeNumber
		   << setw(12) << ElementMaterial_->nset << endl;
}

void CT3::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
	{
		LocationMatrix_[i++] = nodes_[N]->bcode[0];
		LocationMatrix_[i++] = nodes_[N]->bcode[1];
	}
}

void CT3::ElementStiffness(double* Matrix)
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

	double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1];
	double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1];
	double x3 = nodes_[2]->XYZ[0], y3 = nodes_[2]->XYZ[1];

	double area = 0.5 * fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

	double b1 = y2 - y3, c1 = x3 - x2;
	double b2 = y3 - y1, c2 = x1 - x3;
	double b3 = y1 - y2, c3 = x2 - x1;

	double inv_2A = 1.0 / (2.0 * area);

	double B[3][6] = {
		{b1 * inv_2A, 0,           b2 * inv_2A, 0,           b3 * inv_2A, 0          },
		{0,           c1 * inv_2A, 0,           c2 * inv_2A, 0,           c3 * inv_2A},
		{c1 * inv_2A, b1 * inv_2A, c2 * inv_2A, b2 * inv_2A, c3 * inv_2A, b3 * inv_2A}
	};

	double DB[3][6];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 6; j++)
		{
			DB[i][j] = 0.0;
			for (int k = 0; k < 3; k++)
				DB[i][j] += D[i][k] * B[k][j];
		}

	double K[6][6];
	for (int i = 0; i < 6; i++)
		for (int j = 0; j < 6; j++)
			K[i][j] = 0.0;

	for (int i = 0; i < 6; i++)
		for (int j = i; j < 6; j++)
		{
			double val = 0.0;
			for (int k = 0; k < 3; k++)
				val += B[k][i] * DB[k][j];
			K[i][j] = val * factor * area;
		}

	// Copy upper-triangular K to column-by-column storage
	for (int j = 0; j < 6; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = K[i][j];
}

void CT3::ElementStress(double* stress, double* Displacement)
{
	CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);

	double x1 = nodes_[0]->XYZ[0], y1 = nodes_[0]->XYZ[1];
	double x2 = nodes_[1]->XYZ[0], y2 = nodes_[1]->XYZ[1];
	double x3 = nodes_[2]->XYZ[0], y3 = nodes_[2]->XYZ[1];

	double area = 0.5 * fabs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

	double b1 = y2 - y3, c1 = x3 - x2;
	double b2 = y3 - y1, c2 = x1 - x3;
	double b3 = y1 - y2, c3 = x2 - x1;

	double inv_2A = 1.0 / (2.0 * area);

	double B[3][6] = {
		{b1 * inv_2A, 0,           b2 * inv_2A, 0,           b3 * inv_2A, 0          },
		{0,           c1 * inv_2A, 0,           c2 * inv_2A, 0,           c3 * inv_2A},
		{c1 * inv_2A, b1 * inv_2A, c2 * inv_2A, b2 * inv_2A, c3 * inv_2A, b3 * inv_2A}
	};

	double strain[3] = {0.0, 0.0, 0.0};
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 6; j++)
			if (LocationMatrix_[j])
				strain[i] += B[i][j] * Displacement[LocationMatrix_[j] - 1];

	double coeff = material->E / (1.0 - material->Nu * material->Nu);
	stress[0] = coeff * (strain[0] + material->Nu * strain[1]);
	stress[1] = coeff * (material->Nu * strain[0] + strain[1]);
	stress[2] = coeff * ((1.0 - material->Nu) / 2.0) * strain[2];
}
