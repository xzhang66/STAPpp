/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.0, October 14, 2017                                         */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include <stddef.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "Node.h"
#include "Material.h"
#include "SkylineMatrix.h"

using namespace std;

class CDomain;

template <class type> void clear( type* a, unsigned int N );	// Clear an array

//!	Element base class
/*!	All type of element classes should be derived from this base class */
class CElement
{
protected:

//!	Number of nodes per element
	unsigned int NEN;

//!	Nodes of the element
	CNode** nodes;

//!	Material of the element
	CMaterial* ElementMaterial;	//!< Pointer to an element of MaterialSetList[][]
    
//! Location Matrix of the element
    unsigned int* LocationMatrix;

//! Dimension of the location matrix
    unsigned int ND;

public:

//!	Constructor
	CElement() : NEN(0), nodes(NULL), ElementMaterial(NULL) {};

//! Virtual deconstructor
    virtual ~CElement();

//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList) = 0;

//!	Write element data to stream OutputFile
	virtual void Write(ofstream& OutputFile, unsigned int Ele) = 0;

//! Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
    virtual void GenerateLocationMatrix() = 0;
    
//! Calculate the column height, used with the skyline storage scheme
	void CalculateColumnHeight(unsigned int* ColumnHeight); 

//!	Assemble the element stiffness matrix to the global stiffness matrix
	void assembly(double* Matrix, CSkylineMatrix<double>* StiffnessMatrix);

//!	Calculate element stiffness matrix (Upper triangular matrix, stored as an array column by colum)
	virtual void ElementStiffness(double* stiffness) = 0; 

//!	Calculate element stress 
	virtual void ElementStress(double* stress, double* Displacement) = 0;

//!	Return nodes of the element
	inline CNode** GetNodes() { return nodes; }

//!	Return material of the element
	inline CMaterial* GetElementMaterial() { return ElementMaterial; }

//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix() = 0;     

	friend class CDomain;	// Allow class Domain to access its protected member
};
