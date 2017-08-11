/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once

#include <string>
#include <fstream>

using namespace std;

// Structure LoadData is used to store load data
class LoadCaseData
{
public:

	unsigned int nloads;	// Number of concentrated loads in this load case
	unsigned int* node;		// Node number to which this load is applied
	unsigned int* dof;		// Degree of freedom number for this load component
	double* load;			// Magnitude of load

public:

	LoadCaseData() : nloads(0), node(NULL), dof(NULL), load(NULL) {}; 
	~LoadCaseData();

//	Set nloads, and new array node, dof and load
	void Allocate(int num);

//	Read load case data from stream Input
	bool Read(ifstream& Input, int lcase);

//	Write load case data to stream OutputFile
	void Write(ofstream& OutputFile, int lcase);
};
