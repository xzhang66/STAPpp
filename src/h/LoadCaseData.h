/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <string>
#include <fstream>

using namespace std;

//! Class LoadData is used to store load data
class CLoadCaseData
{
public:

	unsigned int nloads;	//!< Number of concentrated loads in this load case
	unsigned int* node;		//!< Node number to which this load is applied
	unsigned int* dof;		//!< Degree of freedom number for this load component
	double* load;			//!< Magnitude of load

public:

	CLoadCaseData() : nloads(0), node(NULL), dof(NULL), load(NULL) {};
	~CLoadCaseData();

//!	Set nloads, and new array node, dof and load
	void Allocate(unsigned int num);

//!	Read load case data from stream Input
	bool Read(ifstream& Input, unsigned int lcase);

//!	Write load case data to stream OutputFile
	void Write(ofstream& OutputFile, unsigned int lcase);
};
