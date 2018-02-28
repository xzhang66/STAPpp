/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

//! Outputer class is used to output results
class COutputter
{
private:

//!	File stream for output
	ofstream OutputFile;

//!	Designed as a single instance class
	static COutputter* _instance;

//! Constructor
    COutputter(string FileName);

public:

//!	Return pointer to the output file stream
	inline ofstream* GetOutputFile() { return &OutputFile; }

//!	Return the single instance of the class
	static COutputter* GetInstance(string FileName = " ");

//!	Output current time and date
	void PrintTime(const struct tm * ptm, COutputter& output);

//!	Output logo and heading 
	void OutputHeading();

//!	Output nodal point data
	void OutputNodeInfo();

//!	Output equation numbers
	void OutputEquationNumber();

//!	Output element data
	void OutputElementInfo();

//!	Output bar element data
	void OutputBarElements(unsigned int EleGrp);

//!	Output load data 
	void OutputLoadInfo(); 

//!	Output displacement data
	void OutputNodalDisplacement();

//!	Output element stresses 
	void OutputElementStress();

//!	Print total system data
	void OutputTotalSystemData();

//! Overload the operator <<
	template <typename T>
	COutputter& operator<<(const T& item) 
	{
		std::cout << item;
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	COutputter& operator<<(CharOstream& (*op)(CharOstream&)) 
	{
		op(std::cout);
		op(OutputFile);
		return *this;
	}

#ifdef _DEBUG_

//!	Print banded and full stiffness matrix for debuging
	void PrintStiffnessMatrix();

//!	Print address of diagonal elements for debuging
	void PrintDiagonalAddress();

//!	Print column heights for debuging
	void PrintColumnHeights();

//!	Print displacement vector for debuging
	void PrintDisplacement();

#endif

};
