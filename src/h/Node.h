/*****************************************************************************/
/*  FEM++ : A C++ FEM code sharing the same input data file with STAP90      */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <iostream>
#include <fstream>

using namespace std;

// #define _DEBUG_    /* Used for print results for debuging purpose */

//	Node class
class Node
{
public:

//	Maximum number of degrees of freedom per node
//	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements NDF = 5 or 6
	const static unsigned int NDF = 3; 

//	Node numer
	unsigned int NodeNumber;

//	x, y and z coordinates of the node
	double XYZ[NDF];

//	Boundary code of each degree of freedom of the node
//		0: The corresponding degree of freedom is active (defined in the global system)
//		1: The corresponding degree of freedom in nonactive (not defined)
//	After call Domain::EquationNumber(), bcode stores the global equation number 
//	corresponding to each degree of freedom of the node
	int bcode[NDF];

//	Constructor
	Node(double X = 0, double Y = 0, double Z = 0);

//	Read nodal point data from stream Input
	bool Read(ifstream& Input, int np);

//	Output nodal point data to stream OutputFile
	void Write(ofstream& OutputFile, int np);

//	Output equation numbers of nodal point to stream OutputFile
	void WriteEquationNo(ofstream& OutputFile, int np);

//	Write nodal displacement
	void WriteNodalDisplacement(ofstream& OutputFile, int np, double* Displacement);

};
