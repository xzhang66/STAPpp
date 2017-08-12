/*****************************************************************************/
/*  FEM++ : A C++ FEM code sharing the same input data file with STAP90      */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool BarMaterial::Read(ifstream& Input, int mset)
{
	Input >> num;	// Number of property set

	if (num != mset + 1)
	{
		cout << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "   Expected set : " << mset + 1 << endl
			 << "   Provided set : " << num << endl;

		return false;
	}

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream OutputFile
void BarMaterial::Write(ofstream& OutputFile, int mset)
{
	cout << setw(5) << mset+1 << setw(16) << E << setw(16) << Area << endl;
	OutputFile << setw(5) << mset+1 << setw(16) << E  << setw(16) << Area << endl;
}
