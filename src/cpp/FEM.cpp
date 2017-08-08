/***************************************************************/
/*  FEM++ ：A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "FEM.h"
#include "Truss.h"

#include <iomanip>
#include <iostream>

using namespace std;

Node::Node(double X = 0, double Y = 0, double Z = 0) 
{
	XYZ[0] = X;		// Coordinates of the node
	XYZ[1] = Y;
	XYZ[2] = Z;

	bcode[0] = 0;	// Boundary codes
	bcode[1] = 0;
	bcode[2] = 0;

	EquationNo[0] = 0;	// Equation numbers 
	EquationNo[1] = 0;
	EquationNo[2] = 0;
};


Domain* Domain::_instance = NULL;

Domain::Domain() : Title(""), MODEX(0), NUMNP(0), NodeList(NULL), NUMEG(0),
ElementList(NULL), NUME(NULL),  NUMMAT(NULL),
MaterialList(NULL), NLCASE(0), NLOAD(NULL), LoadList(NULL), NEQ(0), StiffnessMatrix(NULL), Displacement(NULL),
Force(NULL) {};

Domain* Domain::Instance()
{
	if (!_instance) _instance = new Domain;
	return _instance;
}

// Initial从文件初始化
// 调用Reader的读取数据函数
// 这样做的原因是可以通过继承的方式改写文件读取方法，从而适配不同的文件格式
bool Domain::Initial(FileReader* Reader)
{
	return Reader->ReadFile(this);
}

/***FileReader***/
FileReader::FileReader(string Filename)
{
	Input.open(Filename);
	if (!Input) exit(3);
}

// 读取数据文件
bool FileReader::ReadFile(Domain* FEMData)
{
	int N, LL;

	Input >> FEMData->Title;

	// 读入控制行
	Input >> FEMData->NUMNP >> FEMData->NUMEG >> FEMData->NLCASE
		>> FEMData->MODEX;

	// 读入结点
	FEMData->NodeList = new Node[FEMData->NUMNP];
	Node* NodeList = FEMData->NodeList;
	for (int i = 0; i < FEMData->NUMNP; i++)
	{
		Input >> N;
		if (N != i + 1) return false;
		Input >> NodeList[i].bcode[0] >> NodeList[i].bcode[1] >> NodeList[i].bcode[2]
			>> NodeList[i].XYZ[0] >> NodeList[i].XYZ[1] >> NodeList[i].XYZ[2];
	}

	// 读入工况
	FEMData->LoadList = new LoadData*[FEMData->NLCASE];

	FEMData->NLOAD = new unsigned int[FEMData->NLCASE];
	LoadData** LoadList = FEMData->LoadList;

	unsigned int* NLOAD = FEMData->NLOAD;
	
	for (int i = 0; i < FEMData->NLCASE; i++)
	{
		Input >> LL >> N;

		if (LL != i + 1) 
			return false;
		
		NLOAD[i] = N;
		LoadList[i] = new LoadData[NLOAD[i]];
		for (int j = 0; j < NLOAD[i]; j++)
		{
			Input >>LoadList[i][j].node >> LoadList[i][j].dof >> LoadList[i][j].load;
		}
	}

	// 读入单元组数据
	FEMData->NUME = new unsigned int[FEMData->NUMEG];
	unsigned int* NUME = FEMData->NUME;
	
	FEMData->ElementList = new Element*[FEMData->NUMEG];
	FEMData->NUMMAT = new unsigned int[FEMData->NUMEG];
	unsigned int* NUMMAT = FEMData->NUMMAT;
	FEMData->MaterialList = new Material*[FEMData->NUMEG];

	unsigned int ElementType;
	for (int i = 0; i < FEMData->NUMEG; i++)
	{
		Input >> ElementType >> NUME[i] >> NUMMAT[i];

		// 根据不同的单元类型读取数据 
		// 此处可以尝试改写得更加灵活
		switch (ElementType)
		{
		case 1:
		{
			BarMaterial* MaterialGroup = new BarMaterial[NUMMAT[i]];
			FEMData->MaterialList[i] = MaterialGroup;
			for (int j = 0; j < NUMMAT[i]; j++)
			{
				Input >> N;
				if (N != j + 1) return false;
				Input >> MaterialGroup[j].E >> MaterialGroup[j].Area;
			}
			Bar* ElementList = new Bar[NUME[i]];
			FEMData->ElementList[i] = ElementList;
			for (int j = 0; j < NUME[i]; j++)
			{
				Input >> N;
				if (N != j + 1) return false;
				int MNumber;
				int N1, N2;
				Input >> N1 >> N2 >> MNumber;
				ElementList[j].ElementMaterial = &MaterialGroup[MNumber - 1];
				ElementList[j].nodes[0] = &NodeList[N1 - 1];
				ElementList[j].nodes[1] = &NodeList[N2 - 1];
			}
			break;
		}

		default:
			return false;
		}
	}
	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void Domain::EquationNumber()
{
	NEQ = 0;
	for (int i = 0; i < NUMNP; i++)
	{
		for (int j = 0; j < Node::NDF; j++)
		{
			if (NodeList[i].bcode[j]) 
				NodeList[i].EquationNo[j] = 0;
			else
			{
				NEQ++;
				NodeList[i].EquationNo[j] = NEQ;
			}
		}
	}
}

//	Allocate storage for the one dimensional array storing the global stiffness matrix,
//	and generate the address of diagonal elements
void Domain::AllocateStiffnessMatrix()
{
	Displacement = new double[NEQ];
	Force = new double[NEQ];
	DiagonalAddress = new unsigned int[NEQ + 1];

	// 考虑每个单元的不同自由度I < J, 说明在第(I,J)位置有元素，那么第J列的列高至少为J - I
	for (int i = 0; i < NEQ + 1; i++) 
		DiagonalAddress[i] = 0;

	for (int EG = 0; EG < NUMEG; EG++)
	{
		for (int EN = 0; EN < NUME[EG]; EN++)
		{
			ElementList[EG][EN].ComputeColumnHeight(DiagonalAddress);
		}
	}

	// 根据计算得到的列高DiagonalAddress来分配空间，计算对角元地址
	DiagonalAddress[0] = 1;
	for (int C = 1; C <= NEQ; C++)
	{
		DiagonalAddress[C] = DiagonalAddress[C - 1] + DiagonalAddress[C] + 1;
	}

	StiffnessMatrix = new double[DiagonalAddress[NEQ] - 1];

	for (int i = 0; i < DiagonalAddress[NEQ] - 1; i++) 
		StiffnessMatrix[i] = 0;
}

//组装刚度阵
//对每个单元组（同种单元）分配一次Matrix空间
//调用每个单元的assembly函数
//组装方式在单元中实现
void Domain::AssembleStiffnessMatrix()
{
	for (int EG = 0; EG < NUMEG; EG++)
	{
		unsigned int Space = ElementList[EG][0].LocalMatrixSpace();
		double* Matrix = new double[Space];

		for (int E = 0; E < NUME[EG]; E++)
		{
			ElementList[EG][E].assembly(Matrix);
		}

		delete [] Matrix;
	}
}

//组装力向量
bool Domain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	LoadData* Load = LoadList[LoadCase - 1];

	//先将力置为零
	for (int i = 0; i < NEQ; i++) 
		Force[i] = 0;

	//组装力向量
	for (int FN = 0; FN < NLOAD[LoadCase - 1]; FN++)
	{
		int FreedomDegree = NodeList[Load[FN].node - 1].EquationNo[Load[FN].dof - 1];
		Force[FreedomDegree - 1] += Load[FN].load;
	}

	return true;
}

#ifdef _DEBUG_
//Info
//debug时用的数据输出函数
//当前可以输出刚度阵的一维存储方式，对角元素的位置以及二维的刚度阵
void Domain::Info()
{
	cout << "StiffnessMatrix : " << endl;
	for (int i = 0; i < DiagonalAddress[NEQ] - 1; i++) 
		cout << StiffnessMatrix[i] << " ";
	cout << endl;

	cout << "Address : " << endl;
	for (int i = 0; i < NEQ + 1; i++) 
		cout << DiagonalAddress[i] << " ";
	cout << endl;

	cout << "Matrix : " << endl;
	for (int I = 0; I < NEQ; I++)
	{
		for (int J = 0; J < NEQ; J++)
		{
			int i = I;
			int j = J;
			if (i > j)
			{
				int temp = i;
				i = j;
				j = temp;
			}

			cout << setiosflags(ios::scientific);

			int H = DiagonalAddress[j + 1] - DiagonalAddress[j];
			if (j - i - H >= 0) 
				cout << setw(15) << 0.0;
			else 
				cout << setw(15) << StiffnessMatrix[DiagonalAddress[j] + j - i - 1] << "  ";
		}

		cout << endl;
	}
	cout << endl;

	cout << "U : " << endl;
	for (int I = 0; I < NEQ; I++)
	{
		cout << Displacement[I] << " ";
	}
	cout << endl;
}
#endif