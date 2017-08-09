/***************************************************************/
/*  FEM++ £ºA C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once

#include "FEM.h"

using namespace std;

// Bar element class
class Bar : public Element
{
public:

	Bar();

	friend FileReader;

//	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//	Assemble global stiffness matrix (this should be same for all element type ? to be moved to approriate posion)
	virtual void assembly(double* Matrix);

//	Calculate column height  (this should be same for all element type ? to be moved to approriate posion)
	virtual void ColumnHeight(unsigned int* ColumnHeight);

//	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
};

//	Material class for bar element
class BarMaterial : public Material
{
public:
	double Area;	// Sectional area of a bar element
};