/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "Material.h"

#include <iomanip>
#include <iostream>

using namespace std;

//	Clear an array
template <class type> void clear( type* a, int N )
{
	for (int i = 0; i < N; i++)
		a[i] = 0;
}

CDomain* CDomain::_instance = NULL;

//	Constructor
CDomain::CDomain()
{
	Title[0] = '0';
	MODEX = 0;

	NUMNP = 0;
	NodeList = NULL;
	
	NUMEG = 0;
	ElementTypes = NULL;
	NUME = NULL;
	ElementSetList = NULL;
	
	NUMMAT = NULL;
	MaterialSetList = NULL;
	
	NLCASE = 0;
	NLOAD = NULL;
	LoadCases = NULL;
	
	NEQ = 0;
	NWK = 0;
	MK = 0;

	StiffnessMatrix = NULL;
	ColumnHeights = NULL;
	DiagonalAddress = NULL;

	Force = NULL; 
}

//	Desconstructor
CDomain::~Domain()
{
	delete [] NodeList;

	delete [] ElementTypes;
	delete [] NUME;
	delete [] ElementSetList;

	delete [] NUMMAT;
	delete [] MaterialSetList;

	delete [] NLOAD;
	delete [] LoadCases;

	delete [] StiffnessMatrix;
	delete [] ColumnHeights;
	delete [] DiagonalAddress;

	delete [] Force; 
}

//	Return pointer to the instance of the Domain class
Domain* CDomain::Instance()
{
	if (!_instance) 
		_instance = new Domain;
	
	return _instance;
}

//	Read domain data from the input data file
bool CDomain::ReadData(string FileName, string OutFile)
{
	Input.open(FileName);

	if (!Input) 
	{
		cout << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

	Outputter* Output = Outputter::Instance(OutFile);

//	Read the heading line
	Input.getline(Title, 256);
	Output->OutputHeading();

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data
	ReadNodalPoints();
	Output->OutputNodeInfo();

//	Update equation number
	CalculateEquationNumber();
	Output->OutputEquationNumber();

//	Read load data
	ReadLoadCases();
	Output->OutputLoadInfo();

//	Read element data
	ReadElements();
	Output->OutputElementInfo();

	return true;
}

//	Read nodal point data
bool CDomain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new Node[NUMNP];

//	Loop over for all nodal points
	for (unsigned int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void CDomain::CalculateEquationNumber()
{
	NEQ = 0;
	for (unsigned int np = 0; np < NUMNP; np++)	// Loop over for all node
	{
		for (int dof = 0; dof < Node::NDF; dof++)	// Loop over for DOFs of node np
		{
			if (NodeList[np].bcode[dof]) 
				NodeList[np].bcode[dof] = 0;
			else
			{
				NEQ++;
				NodeList[np].bcode[dof] = NEQ;
			}
		}
	}
}

//	Read load case data
bool CDomain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new LoadCaseData[NLCASE];	// List all load cases

//	Loop over for all load cases
	for (unsigned int lcase = 0; lcase < NLCASE; lcase++)
		if (!LoadCases[lcase].Read(Input, lcase))
			return false;

	return true;
}

// Read element data
bool CDomain::ReadElements()
{

//	Read element group control line
	ElementTypes = new unsigned int[NUMEG];	// Element type of each group
	NUME = new unsigned int[NUMEG];			// Number of elements in each group
	ElementSetList = new Element*[NUMEG];	// Element list in each group

	NUMMAT = new unsigned int[NUMEG];		// Material set number of each group
	MaterialSetList = new Material*[NUMEG];	// Material list in each group

//	Loop over for all element group
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		Input >> ElementTypes[EleGrp] >> NUME[EleGrp] >> NUMMAT[EleGrp];

		switch (ElementTypes[EleGrp])
		{
		case 1:	// Bar element
			ReadBarElementData(EleGrp);
			break;

		default:	// Invalid element type
			return false;
		}
	}
	return true;
}

//	Read bar element data from the input data file
bool CDomain::ReadBarElementData(int EleGrp)
{
//	Read material/section property lines
	MaterialSetList[EleGrp] = new BarMaterial[NUMMAT[EleGrp]];	// Materials for group EleGrp

//	Loop over for all material property sets in group EleGrp
	for (unsigned int mset = 0; mset < NUMMAT[EleGrp]; mset++)
		if (!MaterialSetList[EleGrp][mset].Read(Input, mset))
			return false;

//	Read element data lines
	ElementSetList[EleGrp] = new Bar[NUME[EleGrp]];	// Elements of gorup EleGrp

//	Loop over for all elements in group EleGrp
	for (unsigned int Ele = 0; Ele < NUME[EleGrp]; Ele++)
		if (!ElementSetList[EleGrp][Ele].Read(Input, Ele, MaterialSetList[EleGrp], NodeList))
			return false;

	return true;
}

//	Calculate column heights
void CDomain::CalculateColumnHeights()
{
	clear(ColumnHeights, NEQ);	// Set all elements to zero

	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
		for (unsigned int Ele = 0; Ele < NUME[EleGrp]; Ele++)	//	Loop over for all elements in group EleGrp
			ElementSetList[EleGrp][Ele].CalculateColumnHeight(ColumnHeights);

//	Maximum half bandwidth ( = max(ColumnHeights) + 1 )
	MK = ColumnHeights[0];

	for (unsigned int i=1; i<NEQ; i++)
		if (MK < ColumnHeights[i])
			MK = ColumnHeights[i];

	MK = MK + 1;

#ifdef _DEBUG_
	Outputter* Output = Outputter::Instance();
	Output->PrintColumnHeights();
#endif

}

//	Calculate address of diagonal elements in banded matrix
//	Caution: Address is numbered from 1 !
void CDomain::CalculateDiagnoalAddress()
{
	clear(DiagonalAddress, NEQ + 1);	// Set all elements to zero

//	Calculate the address of diagonal elements
//	M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
	DiagonalAddress[0] = 1;
	for (unsigned int col = 1; col <= NEQ; col++)
		DiagonalAddress[col] = DiagonalAddress[col - 1] + ColumnHeights[col-1] + 1;

//	Number of elements in banded global stiffness matrix
	NWK = DiagonalAddress[NEQ] - DiagonalAddress[0];

#ifdef _DEBUG_
	Outputter* Output = Outputter::Instance();
	Output->PrintDiagonalAddress();
#endif

}

//	Assemble the banded gloabl stiffness matrix
void CDomain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (unsigned int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		unsigned int size = ElementSetList[EleGrp][0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (unsigned int Ele = 0; Ele < NUME[EleGrp]; Ele++)
			ElementSetList[EleGrp][Ele].assembly(Matrix, StiffnessMatrix, DiagonalAddress);

		delete [] Matrix;
	}

#ifdef _DEBUG_
	Outputter* Output = Outputter::Instance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool CDomain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	LoadCaseData* LoadData = &LoadCases[LoadCase - 1];

	clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
	for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
		Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}

//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
void CDomain::AllocateMatrices()
{
//	Allocate for global force/displacement vector
	Force = new double[NEQ];

//	Allocate for column heights 
	ColumnHeights = new unsigned int[NEQ];

//	Allocate for address of diagonal elements
	DiagonalAddress = new unsigned int[NEQ + 1];

//	Calculate column heights
	CalculateColumnHeights();

//	Calculate address of diagonal elements in banded matrix
	CalculateDiagnoalAddress();

//	Allocate for banded global stiffness matrix
	StiffnessMatrix = new double[NWK];
	clear(StiffnessMatrix, NWK);

	Outputter* Output = Outputter::Instance();
	Output->OutputTotalSystemData();
}
