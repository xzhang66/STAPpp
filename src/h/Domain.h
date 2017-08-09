/***************************************************************/
/*  FEM++ £ºA C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#pragma once

#include <string>
#include <fstream>
#include <vector>

#include "Node.h"
#include "Truss.h"
#include "Outputter.h"
#include "Solver.h"

#define _DEBUG_

using namespace std;

//	Clear an array
template <class type> void clear( type* a, int N );

class Domain;

// Structure LoadData is used to store load data
struct LoadData
{
	unsigned node;	// Node number to which this load is applied
	unsigned dof;	// Degree of freedom number for this load component
	double load;	// Magnitude of load
};

//	Domain class : Define the problem domain
//				Only a single instance of Domain class can be created
class Domain
{
private:

//	The instance of the Domain class
	static Domain* _instance;

//	Input file stream for reading data from input data file
	ifstream Input;

//	Heading information for use in labeling the outpu
	string Title; 

//	Solution MODEX
//		0 : Data check only
//		1 : Execution
	int MODEX;

//	Total number of nodal points
	unsigned int NUMNP;

//	List of all nodes in the domain
	Node* NodeList; 

//	Total number of element groups. An element group consists of a convenient
//	collection of elements with same type
	unsigned int NUMEG;

//	Number of elements in each element group
	unsigned int* NUME;

//	Element Set List 
//		ElementSetList[i] - ith element set
//		ElementSetList[i][j] - jth element in ith set
	Element** ElementSetList;

//	Number of different sets of material/section properties in each element group
	unsigned int* NUMMAT;

//	Material set list
//		MaterialSetList[i] - ith material set
//		MaterialSetList[i][j] - jth material in ith set
	Material** MaterialSetList;

//	Number of load cases
	unsigned int NLCASE;

//	List of loads in each load case
	LoadData** LoadList;

//	Number of concentrated loads applied in each load case
	unsigned int* NLOAD;

//	Total number of equations in the system
	unsigned int NEQ;

//	Banded stiffness matrix : A one-dimensional array storing only the elements below the 
//	skyline of the global stiffness matrix. 
	double* StiffnessMatrix;

//	Address of diagonal elements in banded stiffness matrix
	unsigned int* DiagonalAddress;

//	Global nodal displacement vector
	double* Displacement;

//	Global nodal force vector
	double* Force;

//	Constructor
	Domain();

public:

//	Return pointer to the instance of the Domain class
	static Domain* Instance();

	friend Outputter;

//	Return pointer to the banded stiffness matrix
	inline double* GetStiffnessMatrix() { return StiffnessMatrix; }

//	Return pointer to the array storing the address of diagonal elements
	inline unsigned int* GetDiagonalAddress() { return DiagonalAddress; }

//	Return the total number of equations
	inline unsigned int GetNEQ() { return NEQ; }

//	Return pointer to the global nodal force vector
	inline double* GetForce() { return Force; }

//	Return pointer to the global nodal displacement vector
	inline double* GetDisplacement() { return Displacement; }

//	Return the total number of load cases
	inline unsigned int GetNLCASE() { return NLCASE; }

//	Calculate global equation numbers corresponding to every degree of freedom of each node
	void EquationNumber();

//	Allocate storage for the one dimensional array storing the global stiffness matrix,
//	and generate the address of diagonal elements
	void AllocateStiffnessMatrix();

//	Assemble the banded gloabl stiffness matrix
	void AssembleStiffnessMatrix();

//	Assemble the global nodal force vector for load case LoadCase
	bool AssembleForce(unsigned int LoadCase); 

//	Read domain data from the input data file
	bool ReadData(string FileName);

#ifdef _DEBUG_
//	Print debug information
	void Info(); 
#endif

};