/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include <iostream>
#include <iomanip>

#include "Node.h"


Node::Node(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
};


//	Read element data from stream Input
bool Node::Read(ifstream& Input, int np)
{
	int N;

	Input >> N;	// node number
	if (N != np + 1) 
	{
		cout << "*** Error *** Nodes must be inputted in order !" << endl 
			 << "   Expected node number : " << np + 1 << endl
			 << "   Provided node number : " << N << endl;

		return false;
	}

	num = N;

	Input >> bcode[0] >> bcode[1] >> bcode[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];

	return true;
}

//	Output nodal point data to stream OutputFile
void Node::Write(ofstream& OutputFile, int np)
{
	cout << setw(9) << np + 1 << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		 << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
	OutputFile << setw(9) << np + 1 << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
			   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream OutputFile
void Node::WriteEquationNo(ofstream& OutputFile, int np)
{
	cout << setw(9) << np+1 << "       ";
	OutputFile << setw(9) << np+1 << "       ";

	for (int dof = 0; dof < Node::NDF; dof++)	// Loop over for DOFs of node np
	{
		cout << setw(5) << bcode[dof];
		OutputFile << setw(5) << bcode[dof];
	}

	cout << endl;
	OutputFile << endl;
}

//	Write nodal displacement
void Node::WriteNodalDisplacement(ofstream& OutputFile, int np, double* Displacement)
{
	cout << setw(5) << np + 1 << "        ";
	OutputFile << setw(5) << np + 1 << "        ";

	for (int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			cout << setw(18) << 0.0;
			OutputFile << setw(18) << 0.0;
		}
		else
		{
			cout << setw(18) << Displacement[bcode[j] - 1];
			OutputFile << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	cout << endl;
	OutputFile << endl;
}
