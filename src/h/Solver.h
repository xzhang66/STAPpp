/***************************************************************/
/*  FEM++ ：A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once

#include "FEM.h"

using namespace std;

// 求解器类，提供求解器的唯一接口Solve
// 实现新的求解器需要继承此类
// 需要与FEM类的数据存储格式相匹配

class Domain;

class Solver
{
protected:
	Domain* FEMData;

public:
	Solver(Domain* FEMData);

	virtual void Solve() = 0;
};

// LDLT分解求解器
// 与skyline存储格式匹配
class LDLTSolver : public Solver
{
public:
	LDLTSolver(Domain* FEMData) :Solver(FEMData) {};

	void LDLT();                // 对总刚度阵执行LDLT分解，分解后存储于原位置 

	void ComputeDisplacement(); // 使用当前的力向量计算位移 

	virtual void Solve();       // 解FEM问题
};