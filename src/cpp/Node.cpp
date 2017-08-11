/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

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
