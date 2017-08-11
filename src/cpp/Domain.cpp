/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "Domain.h"
#include "Truss.h"

#include <iomanip>
#include <iostream>

using namespace std;

//	Clear an array
template <class type> void clear( type* a, int N )
{
	for (int i = 0; i < N; i++)
		a[i] = 0;
}

LoadCaseData :: ~LoadCaseData()
{
	delete [] node;
	delete [] dof;
	delete [] load;
}

void LoadCaseData :: Allocate(int num)
{
	nloads = num;
	node = new unsigned int[nloads];
	dof = new unsigned int[nloads];
	load = new double[nloads];
}; 

//	Read load case data from stream Input
bool LoadCaseData :: Read(ifstream& Input, int lcase)
{
//	Load case number (LL) and number of concentrated loads in this load case(NL)
	
	int LL, NL;

	Input >> LL >> NL;	

	if (LL != lcase + 1) 
	{
		cout << "*** Error *** Load case must be inputted in order !" << endl 
			 << "   Expected load case : " << lcase + 1 << endl
			 << "   Provided load case : " << LL << endl;

		return false;
	}

	Allocate(NL);

	for (int i = 0; i < NL; i++)
		Input >> node[i] >> dof[i] >> load[i];

	return true;
}


Domain* Domain::_instance = NULL;

//	Constructor
Domain::Domain()
{
	Title = "";
	MODEX = 0;

	NUMNP = 0;
	NodeList = NULL;
	
	NUMEG = 0;
	ElementSetList = NULL;
	
	NUME = NULL;
	NUMMAT = NULL;
	MaterialSetList = NULL;
	
	NLCASE = 0;
	NLOAD = NULL;
	LoadCases = NULL;
	
	NEQ = 0;
	NWK = 0;
	MK = 0;
	StiffnessMatrix = NULL;
	Force = NULL; 
}

//	Return pointer to the instance of the Domain class
Domain* Domain::Instance()
{
	if (!_instance) 
		_instance = new Domain;
	
	return _instance;
}

//	Read domain data from the input data file
bool Domain::ReadData(string FileName)
{
	Input.open(FileName);

	Outputter* Output = Outputter::Instance();

	if (!Input) 
	{
		cout << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

//	Read the heading line
	Input >> Title;

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
bool Domain::ReadNodalPoints()
{

//	Read nodal point data lines
	NodeList = new Node[NUMNP];

//	Loop over for all nodal points
	for (int np = 0; np < NUMNP; np++)
		if (!NodeList[np].Read(Input, np))
			return false;

	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void Domain::CalculateEquationNumber()
{
	NEQ = 0;
	for (int np = 0; np < NUMNP; np++)	// Loop over for all node
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
bool Domain::ReadLoadCases()
{
//	Read load data lines
	LoadCases = new LoadCaseData[NLCASE];

//	Loop over for all load cases
	for (int lcase = 0; lcase < NLCASE; lcase++)
		if (!LoadCases[lcase].Read(Input, lcase))
			return false;

	return true;
}

// Read element data
bool Domain::ReadElements()
{

//	Read element group control line
	ElementTypes = new unsigned int[NUMEG];
	NUME = new unsigned int[NUMEG];
	ElementSetList = new Element*[NUMEG];

	NUMMAT = new unsigned int[NUMEG];
	MaterialSetList = new Material*[NUMEG];

	for (int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		Input >> ElementTypes[EleGrp] >> NUME[EleGrp] >> NUMMAT[EleGrp];

		switch (ElementTypes[EleGrp])
		{
		case 1:	// Bar element
			ReadBarElementData(EleGrp);
			break;

		default:
			return false;
		}
	}
	return true;
}

//	Read bar element data from the input data file
bool Domain::ReadBarElementData(int EleGrp)
{
//	Read material/section property lines
	BarMaterial* MaterialSets = new BarMaterial[NUMMAT[EleGrp]];
	MaterialSetList[EleGrp] = MaterialSets;

//	Loop over for all material property sets
	for (int mset = 0; mset < NUMMAT[EleGrp]; mset++)
		if (!MaterialSetList[EleGrp][mset].Read(Input, mset))
			return false;

//	Read element data lines
	Bar* ElementGroup = new Bar[NUME[EleGrp]];
	ElementSetList[EleGrp] = ElementGroup;

//	Loop over for all elements in group EleGrp
	for (int Ele = 0; Ele < NUME[EleGrp]; Ele++)
		if (!ElementGroup[Ele].Read(Input, Ele, MaterialSets, NodeList))
			return false;

	return true;
}

//	Calculate column heights
void Domain::CalculateColumnHeights()
{
	clear(ColumnHeights, NEQ);	// Set all elements to zero

	for (int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
		for (int Ele = 0; Ele < NUME[EleGrp]; Ele++)	//	Loop over for all elements in group EleGrp
			ElementSetList[EleGrp][Ele].CalculateColumnHeight(ColumnHeights);

//	Maximum half bandwidth ( = max(ColumnHeights) + 1 )
	MK = ColumnHeights[0];
	for (int i=1; i<NEQ; i++)
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
void Domain::CalculateDiagnoalAddress()
{
	clear(DiagonalAddress, NEQ + 1);	// Set all elements to zero

//	Calculate the address of diagonal elements
//	M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
	DiagonalAddress[0] = 1;
	for (int col = 1; col <= NEQ; col++)
		DiagonalAddress[col] = DiagonalAddress[col - 1] + ColumnHeights[col-1] + 1;

//	Number of elements in banded global stiffness matrix
	NWK = DiagonalAddress[NEQ] - DiagonalAddress[0];

#ifdef _DEBUG_
	Outputter* Output = Outputter::Instance();
	Output->PrintDiagonalAddress();
#endif

}

//	Assemble the banded gloabl stiffness matrix
void Domain::AssembleStiffnessMatrix()
{
//	Loop over for all element groups
	for (int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		unsigned int size = ElementSetList[EleGrp][0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

//		Loop over for all elements in group EleGrp
		for (int Ele = 0; Ele < NUME[EleGrp]; Ele++)
			ElementSetList[EleGrp][Ele].assembly(Matrix, StiffnessMatrix, DiagonalAddress);

		delete [] Matrix;
	}

#ifdef _DEBUG_
	Outputter* Output = Outputter::Instance();
	Output->PrintStiffnessMatrix();
#endif

}

//	Assemble the global nodal force vector for load case LoadCase
bool Domain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	LoadCaseData* LoadData = &LoadCases[LoadCase - 1];

	clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
	for (int lnum = 0; lnum < LoadData->nloads; lnum++)
	{
		int dof = NodeList[LoadData->node[lnum] - 1].bcode[LoadData->dof[lnum] - 1];
		Force[dof - 1] += LoadData->load[lnum];
	}

	return true;
}

//	Calculate stresses 
void Domain::CalculateStress()
{
	for (int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		for (int Ele = 0; Ele < NUME[EleGrp]; Ele++)
		{
		}
	}
}

//	Allocate storage for matrices Force, ColumnHeights, DiagonalAddress and StiffnessMatrix
//	and calculate the column heights and address of diagonal elements
void Domain::AllocateMatrices()
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
