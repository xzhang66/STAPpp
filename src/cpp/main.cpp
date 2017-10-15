/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <string>
#include <iostream>

#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
	string InFile = filename + ".dat";
	string OutFile = filename + ".out";

	CDomain* FEMData = CDomain::Instance();

    Clock timer;
    timer.Start();

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cout << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
	FEMData->AllocateMatrices();
    
//  Assemble the banded gloabl stiffness matrix
	FEMData->AssembleStiffnessMatrix();
    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements
	CLDLTSolver* S = new CLDLTSolver(FEMData->GetStiffnessMatrix());
//	S->Solve();
    
//  Perform L*D*L(T) factorization of stiffness matrix
    S->LDLT();

    COutputter* Output = COutputter::Instance();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif
        
//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
//      Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);
            
//      Reduce right-hand-side force vector and back substitute
        S->BackSubstitution(FEMData->GetForce());
            
#ifdef _DEBUG_
        Output->PrintDisplacement(lcase);
#endif
            
        Output->OutputNodalDisplacement(lcase);
    }

    double time_solution = timer.ElapsedTime();

//  Calculate and output stresses of all elements
	Output->OutputElementStress();
    
    double time_stress = timer.ElapsedTime();
    
    timer.Stop();
    
    cout << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
         << "     TIME FOR INPUT PHASE = " << time_input << endl
         << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
         << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
         << "     T O T A L   S O L U T I O N   T I M E = " << time_stress << endl;

    
    ofstream* OutputFile = Output->GetOutputFile();
    *OutputFile << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
                << "     TIME FOR INPUT PHASE = " << time_input << endl
                << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
                << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
                << "     T O T A L   S O L U T I O N   T I M E = " << time_stress << endl;

	return 0;
}
