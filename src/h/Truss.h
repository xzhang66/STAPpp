/************************************************/
/*              FEMCPP 教学程序                 */
/*     清华大学航天航空学院计算动力学教研室     */
/*     Truss.h  杆单元头文件                    */
/*     实现了三维空间杆单元类和材料类的定义     */
/*     作者：宋言                               */
/*     最近修改：2017/06/02                     */
/************************************************/
#pragma once
#include "FEM.h"
using namespace std;
// 杆单元类
class Bar : public Element
{
public:
	Bar();
	friend FileReader;
	virtual void LocalStiffness(double* Matrix);
	virtual void Assembly(double* Matrix);
	virtual void ComputeColumnHeight(unsigned int* ColumnHeight);
	virtual unsigned int LocalMatrixSpace();
};
// 杆的材料
// 增加了截面积变量
class BarMaterial : public Material
{
public:
	double Area;
};