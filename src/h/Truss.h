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

//	Read element data from stream Input
	virtual bool Read(ifstream& Input, int Ele, Material* MaterialSets, Node* NodeList);

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

public:
	
//	Read material data from stream Input
	virtual bool Read(ifstream& Input, int mset);

//	Write material data to Stream OutputFile
//	virtual void Output(ofstream* OutputFile);

};