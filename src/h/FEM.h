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

#include "Outputter.h"
#include "Solver.h"

#define _DEBUG_

using namespace std;

//	Clear an array
template <class type> void clear( type* a, int N );

class FileReader;

class Domain;

/*
   Material base class which only define one data member
   All type of material classes should be derived from this base class
*/
class Material
{
public:
	double E;  // Young's modulus
};

//	Node class
class Node
{
public:

//	Maximum number of degrees of freedom per node
//	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements NDF = 5 or 6
	const static unsigned int NDF = 3; 

//	x, y and z coordinates of the node
	double XYZ[NDF];

//	Boundary code of each degree of freedom of the node
//		0: The corresponding degree of freedom is active (defined in the global system)
//		1: The corresponding degree of freedom in nonactive (not defined)
	int bcode[NDF];

//	Global equation number corresponding to each degree of freedom
	unsigned int EquationNo[NDF]; 

//	Constructor
	Node(double X, double Y, double Z);
};

//	Element base class
//	All type of element classes should be derived from this base class
class Element
{
protected:

//	Number of nodes per element in this type of element
	int NEN;

//	Nodes of the element
	Node** nodes;

//	Material of the element
	Material* ElementMaterial;

public:

//	Constructor
	Element() : NEN(0), nodes(NULL), ElementMaterial(NULL) {};

//	Calculate element stiffness matrix (Upper triangular matrix, stored as an array column by colum)
	virtual void ElementStiffness(double* stiffness) = 0;  

//  Calculate the column height, used with the skyline storage scheme
	void ColumnHeight(unsigned int* ColumnHeight); 

//	Assemble the element stiffness matrix to the global stiffness matrix
	void assembly(double* stiffness);

//	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix() = 0;     

	friend FileReader;
};

// Structure LoadData is used to store load data
struct LoadData
{
	unsigned node;	// Node number to which this load is applied
	unsigned dof;	// Degree of freedom number for this load component
	double load;	// Magnitude of load
};

// FileReader : The reader of the input data file
class FileReader
{
private:
	ifstream Input;

public:

	FileReader(string InputFile);

	virtual bool ReadData(Domain* FEMData);
};

//	Domain class : Define the problem domain
//				Only a single instance of Domain class can be created
class Domain
{
private:

//	The instance of the Domain class
	static Domain* _instance;

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

//	List of all elements in each element group
	Element** ElementList;

//	Number of different sets of material/section properties in each element group
	unsigned int* NUMMAT;

//	List of all material sets
	Material** MaterialList;

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

	friend FileReader;
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

//	Initialize the class data member by reading input data file
	bool Initial(FileReader* Reader); 

#ifdef _DEBUG_
//	Print debug information
	void Info(); 
#endif

};