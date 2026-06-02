#pragma once

#include "Q4.h"

//! Q4R element: Q4 with reduced integration (1-point) + hourglass control
class CQ4R : public CQ4
{
public:
	CQ4R();
	~CQ4R();

	//! Calculate element stiffness matrix with 1-point integration + hourglass stabilization
	virtual void ElementStiffness(double* Matrix);

	//! Calculate element stress at the single integration point (element center)
	virtual void ElementStress(double* stress, double* Displacement);
};
