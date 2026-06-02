#include "Q4R.h"

#include <cmath>
#include <iomanip>

using namespace std;

CQ4R::CQ4R()
{
}

CQ4R::~CQ4R()
{
}

void CQ4R::ElementStiffness(double* Matrix)
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

	// === 1-point integration at element center ===
	double xi0 = 0.0, eta0 = 0.0;
	double dN_dxi_0[4], dN_deta_0[4];
	shapeFunctionDerivatives(xi0, eta0, dN_dxi_0, dN_deta_0);

	double J11_0 = 0.0, J12_0 = 0.0, J21_0 = 0.0, J22_0 = 0.0;
	for (int i = 0; i < 4; i++)
	{
		J11_0 += dN_dxi_0[i] * x[i];
		J12_0 += dN_dxi_0[i] * y[i];
		J21_0 += dN_deta_0[i] * x[i];
		J22_0 += dN_deta_0[i] * y[i];
	}
	double detJ_0 = J11_0 * J22_0 - J12_0 * J21_0;

	double invJ11_0 =  J22_0 / detJ_0;
	double invJ12_0 = -J12_0 / detJ_0;
	double invJ21_0 = -J21_0 / detJ_0;
	double invJ22_0 =  J11_0 / detJ_0;

	double dN_dx_0[4], dN_dy_0[4];
	for (int i = 0; i < 4; i++)
	{
		dN_dx_0[i] = invJ11_0 * dN_dxi_0[i] + invJ12_0 * dN_deta_0[i];
		dN_dy_0[i] = invJ21_0 * dN_dxi_0[i] + invJ22_0 * dN_deta_0[i];
	}

	double B_0[3][8];
	assembleB(B_0, dN_dx_0, dN_dy_0);

	// K_1pt (8x8)
	double K1[8][8];
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			K1[i][j] = 0.0;

	for (int i = 0; i < 8; i++)
	{
		for (int j = i; j < 8; j++)
		{
			double val = 0.0;
			for (int k = 0; k < 3; k++)
			{
				double DB_kj = 0.0;
				for (int m = 0; m < 3; m++)
					DB_kj += D[k][m] * B_0[m][j];
				val += B_0[k][i] * DB_kj;
			}
			K1[i][j] += val * factor * detJ_0 * 4.0;
		}
	}

	// === 2x2 Gauss integration (same as Q4) ===
	const double g = 1.0 / sqrt(3.0);
	const double gauss[4][2] = {
		{-g, -g}, {g, -g}, {g, g}, {-g, g}
	};

	double Kfull[8][8];
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			Kfull[i][j] = 0.0;

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
				Kfull[i][j] += value * factor * detJ;
			}
	}

	// === K_total = K_1pt + r_hg * (K_full - K_1pt) ===
	double r_hg = 0.05;  // 5% hourglass stabilization (softer than full integration)

	double K[8][8];
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			K[i][j] = K1[i][j] + r_hg * (Kfull[i][j] - K1[i][j]);

	// Force symmetry
	for (int i = 0; i < 8; i++)
		for (int j = 0; j < i; j++)
			K[i][j] = K[j][i];

	// Copy upper-triangular K to column-by-column storage
	for (int j = 0; j < 8; j++)
		for (int i = 0; i <= j; i++)
			Matrix[(j + 1) * j / 2 + j - i] = K[i][j];
}

void CQ4R::ElementStress(double* stress, double* Displacement)
{
	// Stress at the single integration point (element center)
	CQ4::ElementStress(stress, Displacement);
}
