/*****************************************************************************/
/*  FEM++ : A C++ FEM code sharing the same input data file with STAP90      */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

// Bar element class
class Bar : public Element
{
public:

//	Constructor
	Bar();

//	Desconstructor
	~Bar();

//	Read element data from stream Input
	virtual bool Read(ifstream& Input, int Ele, Material* MaterialSets, Node* NodeList);

//	Write element data to stream OutputFile
	virtual void Write(ofstream& OutputFile, int Ele);

//	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

//	Calculate element stress 
	virtual void ElementStress(double* stress, double* Displacement);

//	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
};
