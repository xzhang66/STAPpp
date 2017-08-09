/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once

#include "Element.h"

using namespace std;

// Bar element class
class Bar : public Element
{
public:

	Bar();

//	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
};

//	Material class for bar element
class BarMaterial : public Material
{
public:
	double Area;	// Sectional area of a bar element
};