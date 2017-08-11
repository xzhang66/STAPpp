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

	Outputter* Output = Outputter::Instance(OutFile);
	Output->OutputHeading();

	Domain* FEMData = Domain::Instance();
	FEMData->SetOutputFile(Output->GetOutputFile());

	if (!FEMData->ReadData(InFile))
	{
		cout << "*** Error *** Data input failed!" << endl;
		exit(1);
	}

	Output->OutputNodeInfo();

	FEMData->CalculateEquationNumber();
	Output->OutputEquationNumber();

	Output->OutputElementInfo();

	FEMData->AllocateMatrices();
	
	FEMData->AssembleStiffnessMatrix();
	
	LDLTSolver* S = new LDLTSolver(FEMData);
	S->Solve();

	return 0;
}