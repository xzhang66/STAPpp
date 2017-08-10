/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
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

//	File stream for output
	ofstream OutputFile;

protected:

//	Constructor
	Outputter(string FileName);

//	Designed as a single instance class
	static Outputter* _instance;

//	Print program logo
	void PrintLogo(ostream& output);

public:

//	Return the single instance of the class
	static Outputter* Instance(string FileName = " ");

//	Return pointer to the output file stream
	inline ofstream* GetOutputFile() { return &OutputFile; }

//	Output current time and date
	void PrintTime(const struct tm * ptm, ostream& output);

//	Output logo and heading 
	void OutputHeading();

//	Output nodal point data
	void OutputNodeInfo();

//	Output load data for load case LoadCase
	void OutputLoadInfo(int LoadCase); 

//	Output displacement data
	void OutputDisplacement();
};