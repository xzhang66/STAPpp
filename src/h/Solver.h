/************************************************/
/*              FEMCPP 教学程序                 */
/*     清华大学航天航空学院计算动力学教研室     */
/*     Solver.h  求解器类头文件                 */
/*     作者：宋言                               */
/*     最近修改：2017/06/05                     */
/************************************************/
#pragma once
#include "FEM.h"
using namespace std;
// 求解器类，提供求解器的唯一接口Solve
// 实现新的求解器需要继承此类
// 需要与FEM类的数据存储格式相匹配
class FEM;
class Solver
{
protected:
	FEM* FEMData;
public:
	Solver(FEM* FEMData);
	virtual void Solve() = 0;
};
// LDLT分解求解器
// 与skyline存储格式匹配
class LDLTSolver : public Solver
{
public:
	LDLTSolver(FEM* FEMData) :Solver(FEMData) {};
	void LDLT();                // 对总刚度阵执行LDLT分解，分解后存储于原位置 
	void ComputeDisplacement(); // 使用当前的力向量计算位移 
	virtual void Solve();       // 解FEM问题
};