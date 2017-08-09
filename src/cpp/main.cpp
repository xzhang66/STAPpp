/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "Domain.h"
#include "Truss.h"
#include "Outputter.h"

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: FEM++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
	string InFile = filename + ".dat";
	string OutFile = filename + ".out";

	Domain* FEMData = Domain::Instance();

	if (!FEMData->ReadData(InFile))
	{
		cout << "*** Error *** Data input failed!" << endl;
		exit(1);
	}

	Outputter* Output = Outputter::Instance(OutFile);
	Output->OutputHeading();
	Output->OutputNodeInfo();

	FEMData->EquationNumber();
	
	FEMData->AllocateStiffnessMatrix();
	
	FEMData->AssembleStiffnessMatrix();
	
	LDLTSolver* S = new LDLTSolver(FEMData);
	S->Solve();
	
#ifdef _DEBUG_
	FEMData->Info();
#endif

	return 0;
}