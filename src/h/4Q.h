/*****************************************************************************/
// the following code is coded by Fan Yang 
// date:14th April, 2019
/*****************************************************************************/

#pragma once

#include "Element.h"
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;

//! 4Q element class
class C4Q : public CElement
{
public:

	//!	Constructor
	C4Q();

	//!	Desconstructor
	~C4Q();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

	//! Generate location matrix: the global equation number that corresponding to each DOF of the element
	//	Caution:  Equation number is numbered from 1 !
	virtual void GenerateLocationMatrix();

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

	//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();
	C4Q::MatrixXd getB(double psi, double eta, MatrixXd C);

};
