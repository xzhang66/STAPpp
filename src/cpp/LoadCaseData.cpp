/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "LoadCaseData.h"

#include <iomanip>
#include <iostream>

using namespace std;

CLoadCaseData :: ~CLoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
}

void CLoadCaseData :: Allocate(unsigned int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
}; 

//	Read load case data from stream Input
bool CLoadCaseData :: Read(ifstream& Input, unsigned int lcase)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL)
	
	unsigned int LL, NL;

	Input >> LL >> NL;	

	if (LL != lcase + 1) 
	{
		cout << "*** Error *** Load case must be inputted in order !" << endl 
			 << "   Expected load case : " << lcase + 1 << endl
			 << "   Provided load case : " << LL << endl;

		return false;
	}

	Allocate(NL);

	for (unsigned int i = 0; i < NL; i++)
		Input >> node[i] >> dof[i] >> load[i];

	return true;
}

//	Write load case data to stream OutputFile
void CLoadCaseData::Write(ofstream& OutputFile, unsigned int lcase)
{
	for (unsigned int i = 0; i < nloads; i++)
	{
		cout << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
		OutputFile << setw(7) << node[i] << setw(13) << dof[i]  << setw(19) << load[i] << endl;
	}
}
