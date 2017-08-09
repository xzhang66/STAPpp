/***************************************************************/
/*  FEM++ £ºA C++ finite element method code for teaching      */
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

// Outputer class for outputing results
class Outputter
{
private:

//	Designed as a single instance class
	static Outputter* _instance;

//	Print program logo
void PrintLogo(ostream& output);

	ofstream OutputFile;		// File stream for output
	Outputter(string FileName);

public:

//	Constructor
	static Outputter* Instance(string FileName = " ");

//	Output logo and heading 
	void OutputHeading();

//	Output nodal point data
	void OutputNodeInfo();

//	Output load data for load case LoadCase
	void OutputLoadInfo(int LoadCase); 

//	Output displacement data
	void OutputDisplacement();
};