/************************************************/
/*              FEMCPP 教学程序                 */
/*     清华大学航天航空学院计算动力学教研室     */
/*     Solver.cpp  求解器类源文件               */
/*     作者：宋言                               */
/*     最近修改：2017/06/12                     */
/************************************************/
#include "Solver.h"
#include <math.h>
#include <iostream>
#define Error 1e-24  // 用来判断刚度阵是否正定
using namespace std;
int MIN(int I, int J)
{
	if (I < J) return I;
	else return J;
};
int MAX(int I, int J)
{
	if (I > J) return I;
	else return J;
};
Solver::Solver(FEM* FEMData) : FEMData(FEMData) {};
// LDLT分解
void LDLTSolver::LDLT()
{
	double* K = FEMData->GetStiffnessMatrix();
	unsigned int* Address = FEMData->GetDiagonalAddress();
	unsigned int N = FEMData->GetFreedom();
	for (int j = 0; j < N; j++)      //对列循环
	{
		double* Columnj = &K[Address[j] - 1];
		int ColumnNumberj = Address[j + 1] - Address[j];   //此列的非零元素数量
		int Heightj = j - ColumnNumberj + 1;               //第j列的最高非零元素位置
		for (int i = Heightj; i <= j; i++)                 //对所有列的非零元素循环
		{
			int ColumnNumberi = Address[i + 1] - Address[i];  //第i列的非零元素数量
			double* Columni = &K[Address[i] - 1];
			int CurPostion = Address[j] + j - i - 1;
			double C = 0;
			int Heighti = i - ColumnNumberi + 1;           //第i列的最高非零元素位置
			int Height = MAX(Heighti, Heightj);            //Height为i,j两列最低的高度
			for (int M = Height; M < i; M++)
			{
				int AddressI = Address[i] + i - M - 1;
				int AddressJ = Address[j] + j - M - 1;
				C += K[AddressI] * K[AddressJ] * K[Address[M] - 1];
			}
			if (i == j)
			{
				K[CurPostion] = K[CurPostion] - C;
				if (abs(K[CurPostion]) < Error)
				{
					cout << "在组装第" << i + 1 << "个自由度时，刚度阵不正定" << endl;
					exit(4);
				}
			}
			else K[CurPostion] = (K[CurPostion] - C) / K[Address[i] - 1];
		}
	}
};
// 解方程，计算位移
void LDLTSolver::ComputeDisplacement()
{
	double* Force = FEMData->GetForce();        //力向量
	double* K = FEMData->GetStiffnessMatrix();  //已经做过LDLT分解的刚度阵
	double* U = FEMData->GetDisplacement();     //位移
	unsigned int* Address = FEMData->GetDiagonalAddress();  //对角元素位置
	unsigned int Freedom = FEMData->GetFreedom(); //自由度数
	// L * V = F , V与F占用同样的位置
	for (int i = 0; i < Freedom; i++)
	{
		int Height = Address[i + 1] - Address[i]; //第i行非零元素个数
		int CurPos = Address[i + 1] - 2;
		for (int M = i - Height + 1; M < i; M++)
		{
			Force[i] -= K[CurPos] * Force[M];
			CurPos--;
		}
	}
	// D * S = V,  S V F占用同样的位置
	for (int i = 0; i < Freedom; i++)
	{
		Force[i] = Force[i] / K[Address[i] - 1];
	}
	// LT * U = V
	for (int i = Freedom - 1; i >= 0; i--)
	{
		double C = 0;
		for (int M = Freedom - 1; M > i; M--)
		{
			int Height = Address[M + 1] - Address[M];
			if (M - Height + 1 <= i) 
				C += K[Address[M] - 1 + M - i] * Force[M];
		}
		Force[i] = Force[i] - C;
	}
	for (int i = 0; i < Freedom; i++) U[i] = Force[i];
};
void LDLTSolver::Solve()
{ 
	Outputter* Output = Outputter::Instance();
	LDLT();
	for (int i = 0; i < FEMData->GetLoadCaseNumber(); i++)
	{
		FEMData->AssemblyForce(i + 1);
		ComputeDisplacement();
		Output->OutputLoadInfo(i + 1);
		Output->OutputDisplacement();
	}
	return; 
};