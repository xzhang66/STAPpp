/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once

using namespace std;

//	Node class
class Node
{
public:

//	Maximum number of degrees of freedom per node
//	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements NDF = 5 or 6
	const static unsigned int NDF = 3; 

//	x, y and z coordinates of the node
	double XYZ[NDF];

//	Boundary code of each degree of freedom of the node
//		0: The corresponding degree of freedom is active (defined in the global system)
//		1: The corresponding degree of freedom in nonactive (not defined)
	int bcode[NDF];

//	Global equation number corresponding to each degree of freedom
	unsigned int EquationNo[NDF]; 

//	Constructor
	Node(double X = 0, double Y = 0, double Z = 0);
};
