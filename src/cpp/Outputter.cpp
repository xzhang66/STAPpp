/***************************************************************/
/*  FEM++ ：A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "FEM.h"
#include "Outputter.h"

#include <iostream>
#include <iomanip>
#include <time.h>

using namespace std;

// 时间输出函数
// 可以给定输出流output
void dsptime(const struct tm * ptm, ostream& output)
{
	char *pxq[] = { "日","一","二","三","四","五","六" };
	output << ptm->tm_year + 1900 << "年" << ptm->tm_mon + 1 << "月" << ptm->tm_mday << "日 ";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " ";
	output << " 星期" << pxq[ptm->tm_wday] << endl;
}

Outputter* Outputter::_instance = NULL;

// 构造函数
Outputter::Outputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile) 
		exit(3);
}

// 单例函数
Outputter* Outputter::Instance(string FileName)
{
	if (!_instance) _instance = new Outputter(FileName);
	return _instance;
}

// 打印文件头
void Outputter::OutputLogo()
{
	cout << "***********************************************************" << endl;
	cout << "xxxxxx  xxxxxx  xxx       xx      xxxxx  xxxxxx   xxxxxx" << endl;
	cout << "xx      xx      xxxx     xxx    xxx      xx   xx  xx   xx" << endl;
	cout << "xx      xx      xxxx     x x    xx       xx    xx xx    xx" << endl;
	cout << "xx      xx      xx x    xx x   xx        xx    xx xx    xx" << endl;
	cout << "xx      xx      xx xx   xx x   xx        xx    xx xx    xx" << endl;
	cout << "xxxxxx  xxxxxx  xx xx   x  x   xx        xx   xx  xx   xx" << endl;
	cout << "xx      xx      xx  xx xx  x   xx        xxxxxx   xxxxxx" << endl;
	cout << "xx      xx      xx  xx x   x   xx        xx       xx" << endl;
	cout << "xx      xx      xx   xxx   x   xx        xx       xx" << endl;
	cout << "xx      xx      xx   xxx   x    xxx   x  xx       xx" << endl;
	cout << "xx      xxxxxxx xx    x    x      xxxxx  xx       xx" << endl;
	cout << "***********************************************************" << endl;
	cout << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;

	OutputFile << "***********************************************************" << endl;
	OutputFile << "xxxxxx  xxxxxx  xxx       xx      xxxxx  xxxxxx   xxxxxx" << endl;
	OutputFile << "xx      xx      xxxx     xxx    xxx      xx   xx  xx   xx" << endl;
	OutputFile << "xx      xx      xxxx     x x    xx       xx    xx xx    xx" << endl;
	OutputFile << "xx      xx      xx x    xx x   xx        xx    xx xx    xx" << endl;
	OutputFile << "xx      xx      xx xx   xx x   xx        xx    xx xx    xx" << endl;
	OutputFile << "xxxxxx  xxxxxx  xx xx   x  x   xx        xx   xx  xx   xx" << endl;
	OutputFile << "xx      xx      xx  xx xx  x   xx        xxxxxx   xxxxxx" << endl;
	OutputFile << "xx      xx      xx  xx x   x   xx        xx       xx" << endl;
	OutputFile << "xx      xx      xx   xxx   x   xx        xx       xx" << endl;
	OutputFile << "xx      xx      xx   xxx   x    xxx   x  xx       xx" << endl;
	OutputFile << "xx      xxxxxxx xx    x    x      xxxxx  xx       xx" << endl;
	OutputFile << "***********************************************************" << endl;
	OutputFile << "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;

	Domain* FEMData = Domain::Instance();

	cout << "TITLE : " << FEMData->Title << endl;
	OutputFile << "TITLE : " << FEMData->Title << endl;

	time_t nowtime;
	struct tm *local = new struct tm;
	nowtime = time(NULL);

	localtime_s(local, &nowtime);
	dsptime(local, cout);
	dsptime(local, OutputFile);
}

// 打印结点信息
void Outputter::OutputNodeInfo()
{
	Domain* FEMData = Domain::Instance();

	Node* NodeList = FEMData->NodeList;

	int Page = 30;      // 每页行数

	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);

	cout << "*********************  N O D E **************************" << endl;
	OutputFile << "*********************  N O D E **************************" << endl;

	for (int i = 0; i < FEMData->NUMNP; i++)
	{
		if (i % Page == 0)
		{
			cout << "No. ......... X ........... Y ........... Z ..........." << endl;
			OutputFile << "No. ......... X ........... Y ........... Z ..........." << endl;
		}

		cout << setw(12) << i + 1 << setw(14) << NodeList[i].XYZ[0] << setw(14) << NodeList[i].XYZ[1] 
		     << setw(14) << NodeList[i].XYZ[2] << endl;
		OutputFile << setw(12) << i + 1 << setw(14) << NodeList[i].XYZ[0] << setw(14) << NodeList[i].XYZ[1] 
		           << setw(14) << NodeList[i].XYZ[2] << endl;
	}
}

// 打印载荷信息
void Outputter::OutputLoadInfo(int LoadCase)
{
	Domain* FEMData = Domain::Instance();

	unsigned int* NLOAD = FEMData->NLOAD;

	if (LoadCase > FEMData->NLCASE) 
		return;

	LoadData* Load = FEMData->LoadList[LoadCase - 1];

	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);

	cout << "****************** L O A D C A S E " << LoadCase <<  " ********************" << endl;
	OutputFile << "****************** L O A D C A S E " << LoadCase <<  " ********************" << endl;

	cout << "No. ...... NUMNP .. DIR ....... ....Force .........." << endl;
	OutputFile << "No. ...... NUMNP .. DIR ....... ....Force .........." << endl;

	for (int i = 0; i < NLOAD[LoadCase - 1]; i++)
	{
		cout << setw(10) << i + 1 << setw(14) << Load[i].node << setw(12) << Load[i].dof 
			 << setw(17) << Load[i].load << endl;
		OutputFile << setw(10) << i + 1 << setw(14) << Load[i].node << setw(12) << Load[i].dof 
			       << setw(17) << Load[i].load << endl;
	}
}

// 输出位移
void Outputter::OutputDisplacement()
{
	Domain* FEMData = Domain::Instance();

	int Page = 30;

	Node* NodeList = FEMData->NodeList;

	double* Displacement = FEMData->Displacement;

	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);

	cout << "************* D I S P L A C E l M E N T *****************" << endl;
	OutputFile << "************* D I S P L A C E l M E N T *****************" << endl;

	for (int i = 0; i < FEMData->NUMNP; i++)
	{
		if (i % Page == 0)
		{
			cout << "No. ......... X ........... Y ........... Z ..........." << endl;
			OutputFile << "No. ......... X ........... Y ........... Z ..........." << endl;
		}

		cout << setw(12) << i + 1;
		OutputFile << setw(12) << i + 1;

		for (int j = 0; j < 3; j++)
		{
			if (NodeList[i].EquationNo[j] == 0)
			{
				cout << setw(14) << 0.0;
				OutputFile << setw(14) << 0.0;
			}
			else
			{
				cout << setw(14) << Displacement[NodeList[i].EquationNo[j] - 1];
				OutputFile << setw(14) << Displacement[NodeList[i].EquationNo[j] - 1];
			}
		}

		cout << endl;
		OutputFile << endl;
	}
}