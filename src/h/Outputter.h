/***************************************************************/
/*  FEM++ ：A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

// 输出器文件，使用单例模式，可以在各出进行调用
class Outputter
{
private:
	static Outputter* _instance;

	ofstream OutputFile;      //数据文件输出流
	Outputter(string FileName);

public:
	static Outputter* Instance(string FileName = " ");

	void OutputLogo();

	void OutputNodeInfo();

	void OutputLoadInfo(int LoadCase);  //输出第LoadCase个工况的信息

	void OutputDisplacement();          //输出位移
};