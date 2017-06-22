#pragma once
/************************************************/
/*              FEMCPP 教学程序                 */
/*     清华大学航天航空学院计算动力学教研室     */
/*     FEM.h  有限元数据类型头文件              */
/*     作者：宋言                               */
/*     最近修改：2017/06/18                     */
/************************************************/
#include <string>
#include <fstream>
#include <vector>
#include "Outputter.h"
#include "Solver.h"
#define _DEBUG_
using namespace std;
class FileReader;
class FEM;
// 材料和截面类
/************************/
// 材料和截面性质的基类，只存储杨氏模量一个参数。新的材料通过继承实现。
/************************/
class Material
{
public:
	double E;  //杨氏模量
};
// 结点类
// 三维结点类，存储了XYZ三个平动自由度的信息
class Node
{
public:
	double XYZ[3];  // XYZ坐标
	int Fix[3];     // 约束条件，1：约束。0：自由。
	unsigned int Freedom[3];  // 三个分量的自由度编号
	static unsigned int Dimension;  // 静态：结点的自由度数
	Node(double X, double Y, double Z);
};
// 单元类
// 单元基类，提供了单元的基本信息和接口
// 实现具体单元应继承此类
class Element
{
protected:
	unsigned int NodeNumber;    //单元的节点数
	Node** NodeList;            //单元的结点指针数组
	Material* ElementMaterial;  //单元的材料与截面
public:
	Element() :NodeNumber(0), NodeList(NULL), ElementMaterial(NULL) {};
	virtual void LocalStiffness(double* Matrix) = 0;  //计算单元刚度阵，存储在Matrix中
	virtual void Assembly(double* Matrix) = 0;        //将单元刚度阵组装入总刚，Matrix的存储格式需要与LocalStiffness相配
	virtual void ComputeColumnHeight(unsigned int* ColumnHeight) = 0; //计算列高，在skyline存储方法中使用
	virtual unsigned int LocalMatrixSpace() = 0;     //返回单元刚度阵所占的空间大小，需要与Matrix的存储格式相配
	friend FileReader;
};
// LoadData 用来暂存力信息
struct LoadData
{
	unsigned NodeNumber;   //结点号
	unsigned Direction;    //力的方向
	double Force;          //力的大小
};
// FileReader 文件读入器
class FileReader
{
private:
	ifstream Input;
public:
	FileReader(string InputFile);
	virtual bool ReadFile(FEM* FEMData);
};
// FEM类，单例，存储有限元算法用到的数据
class FEM
{
private:
	static FEM* _instance;
	string Title;                     //标题
	int Mode;						  //求解模式
	unsigned int NodeNumber;          //结点数
	Node* NodeList;                   //结点
	unsigned int ElementGroupNumber;  //单元组数
	unsigned int* ElementNumber;      //每个单元组的单元数
	Element** ElementGroupList;       //单元
	unsigned int* MaterialNumber;     //每个材料组的材料数
	Material** MaterialGroupList;     //材料
	LoadData** LoadList;              //载荷
	unsigned int LoadCaseNumber;      //工况数
	unsigned int* LoadNumber;         //每个工况中载荷数量

	unsigned int Freedom;             //自由度数
	double* StiffnessMatrix;          //总刚度阵
	unsigned int* DiagonalAddress;    //刚度阵对角元素的地址
	double* Displacement;             //位移
	double* Force;                    //力向量
	FEM();
public:
	static FEM* Instance();
	friend FileReader;
	friend Outputter;
	inline double* GetStiffnessMatrix() { return StiffnessMatrix; }
	inline unsigned int* GetDiagonalAddress() { return DiagonalAddress; }
	inline unsigned int GetFreedom() { return Freedom; }
	inline double* GetForce() { return Force; }
	inline double* GetDisplacement() { return Displacement; }
	inline unsigned int GetLoadCaseNumber() { return LoadCaseNumber; }
	void GenerateFreedom();           //计算自由度数
	void AllocateStiffnessMatrix();   //分配刚度阵的存储空间
	void AssemblyStiffnessMatrix();   //组装刚度阵
	bool AssemblyForce(unsigned int LoadCase);   //组装第LoadCase个力向量，成功则返回true
	bool Initial(FileReader* Reader); //从文件读入器初始化，读取基本信息
#ifdef _DEBUG_
	void Info();  //debug时用到的数据输出函数
#endif
};