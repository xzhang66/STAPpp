#pragma once

#include "Element.h"

//! T3 element class: 3-node constant strain triangle (CST) plane stress element
class CT3 : public CElement
{
public:
	CT3();
	~CT3();

	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
	virtual void Write(COutputter& output);
	virtual void GenerateLocationMatrix();
	virtual void ElementStiffness(double* Matrix);
	virtual void ElementStress(double* stress, double* Displacement);
};
