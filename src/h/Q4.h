#pragma once

#include "Element.h"

using namespace std;

//! Shape function derivatives for 4-node bilinear element
void shapeFunctionDerivatives(double xi, double eta, double dN_dxi[4], double dN_deta[4]);

//! Assemble strain-displacement matrix B from shape function derivatives
void assembleB(double B[3][8], const double dN_dx[4], const double dN_dy[4]);

//! Q4 element class: bilinear 4-node plane stress element
class CQ4 : public CElement
{
public:
	//! Constructor
	CQ4();

	//! Destructor
	~CQ4();

	//! Read element data from stream Input
	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

	//! Write element data to stream
	virtual void Write(COutputter& output);

	//! Generate location matrix for the two in-plane degrees of freedom
	virtual void GenerateLocationMatrix();

	//! Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//! Calculate element stress at the element center
	virtual void ElementStress(double* stress, double* Displacement);
};