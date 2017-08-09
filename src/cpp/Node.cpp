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
    
    EquationNo[0] = 0;	// Equation numbers
    EquationNo[1] = 0;
    EquationNo[2] = 0;
};
