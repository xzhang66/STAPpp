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
	LoadList = NULL;
	
	NEQ = 0;
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

	if (!Input) 
	{
		cout << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}

//	Read the heading line
	Input >> Title;

//	Read the control line
	Input >> NUMNP >> NUMEG >> NLCASE >> MODEX;

//	Read nodal point data lines
	NodeList = new Node[NUMNP];
	
	int N;

//	Loop over for all nodal points
	for (int np = 0; np < NUMNP; np++)
	{
		Input >> N;	// node number
		if (N != np + 1) 
		{
			cout << "*** Error *** Nodes must be inputted in order !" << endl 
				 << "   Expected node number : " << np + 1 << endl
				 << "   Provided node number : " << N << endl;

			return false;
		}

		Input >> NodeList[np].bcode[0] >> NodeList[np].bcode[1] >> NodeList[np].bcode[2]
			  >> NodeList[np].XYZ[0] >> NodeList[np].XYZ[1] >> NodeList[np].XYZ[2];
	}
	
	CalculateEquationNumber();

//	Read load data lines
	LoadList = new LoadData*[NLCASE];
	NLOAD = new unsigned int[NLCASE];
	
	int LL, NL;

//	Loop over for all load cases
	for (int lcase = 0; lcase < NLCASE; lcase++)
	{
//		Load case number (LL) and number of concentrated loads in this load case(NL)
		Input >> LL >> NL;	

		if (LL != lcase + 1) 
		{
			cout << "*** Error *** Load case must be inputted in order !" << endl 
				 << "   Expected load case : " << lcase + 1 << endl
				 << "   Provided load case : " << LL << endl;

			return false;
		}

		NLOAD[lcase] = NL;

		LoadList[lcase] = new LoadData[NLOAD[lcase]];
		for (int load = 0; load < NLOAD[lcase]; load++)
			Input >>LoadList[lcase][load].node >> LoadList[lcase][load].dof >> LoadList[lcase][load].load;
	}

//	Read element group control line
	NUME = new unsigned int[NUMEG];
	ElementSetList = new Element*[NUMEG];

	NUMMAT = new unsigned int[NUMEG];
	MaterialSetList = new Material*[NUMEG];

	unsigned int ElementType;
	for (int EleGrp = 0; EleGrp < NUMEG; EleGrp++)
	{
		Input >> ElementType >> NUME[EleGrp] >> NUMMAT[EleGrp];

//		Will try to move to element class
		switch (ElementType)
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
	int N;

//	Read material/section property lines
	BarMaterial* MaterialGroup = new BarMaterial[NUMMAT[EleGrp]];
	MaterialSetList[EleGrp] = MaterialGroup;

//	Loop over for all property sets
	for (int mset = 0; mset < NUMMAT[EleGrp]; mset++)
	{
		Input >> N;	// Number of property set

		if (N != mset + 1)
		{
			cout << "*** Error *** Material sets must be inputted in order !" << endl 
				 << "   Expected set : " << mset + 1 << endl
				 << "   Provided set : " << N << endl;

			return false;
		}
				
		Input >> MaterialGroup[mset].E >> MaterialGroup[mset].Area;	// Young's modulus and section area
	}

//	Read element data lines
	Bar* ElementGroup = new Bar[NUME[EleGrp]];
	ElementSetList[EleGrp] = ElementGroup;

//	Loop over for all elements in group EleGrp
	for (int Ele = 0; Ele < NUME[EleGrp]; Ele++)
	{
		Input >> N;	// element number

		if (N != Ele + 1)
		{
			cout << "*** Error *** Elements must be inputted in order !" << endl 
				 << "   Expected element : " << Ele + 1 << endl
				 << "   Provided element : " << N << endl;

			return false;
		}

		int MSet;	// Material property set number
		int N1, N2;	// Left node number and right node number

		Input >> N1 >> N2 >> MSet;
		ElementGroup[Ele].ElementMaterial = &MaterialGroup[MSet - 1];
		ElementGroup[Ele].nodes[0] = &NodeList[N1 - 1];
		ElementGroup[Ele].nodes[1] = &NodeList[N2 - 1];
	}
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

//	Calculate column heights
void Domain::CalculateColumnHeights()
{
	clear(ColumnHeights, NEQ);	// Set all elements to zero

	for (int EleGrp = 0; EleGrp < NUMEG; EleGrp++)		//	Loop over for all element groups
		for (int Ele = 0; Ele < NUME[EleGrp]; Ele++)	//	Loop over for all elements in group EleGrp
			ElementSetList[EleGrp][Ele].CalculateColumnHeight(ColumnHeights);
}

//	Calculate address of diagonal elements in banded matrix
void Domain::CalculateDiagnoalAddress()
{
	clear(DiagonalAddress, NEQ + 1);	// Set all elements to zero

//	Calculate the address of diagonal elements
//	M(0) = 1;  M(i+1) = M(i) + H(i) + 1 (i = 0:NEQ)
	DiagonalAddress[0] = 1;
	for (int col = 1; col <= NEQ; col++)
		DiagonalAddress[col] = DiagonalAddress[col - 1] + ColumnHeights[col-1] + 1;
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
}

//	Assemble the global nodal force vector for load case LoadCase
bool Domain::AssembleForce(unsigned int LoadCase)
{
	if (LoadCase > NLCASE) 
		return false;

	LoadData* Load = LoadList[LoadCase - 1];

	clear(Force, NEQ);

//	Loop over for all concentrated loads in load case LoadCase
	for (int lnum = 0; lnum < NLOAD[LoadCase - 1]; lnum++)
	{
		int dof = NodeList[Load[lnum].node - 1].bcode[Load[lnum].dof - 1];
		Force[dof - 1] += Load[lnum].load;
	}

	return true;
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
	StiffnessMatrix = new double[DiagonalAddress[NEQ] - 1];
	clear(StiffnessMatrix, DiagonalAddress[NEQ] - 1);
}

#ifdef _DEBUG_

//	Print debug information
void Domain::Info()
{
	cout << endl << endl;
	cout << "**************** Debug Infomation ******************" << endl << endl;
    
	cout << "StiffnessMatrix : " << endl;
	for (int i = 0; i < DiagonalAddress[NEQ] - 1; i++)
    {
        if (i%6 == 0)
            cout << endl;
        
		cout << setw(15) << StiffnessMatrix[i];
    }
    cout << endl << endl;

	cout << "Address : " << endl;
	for (int i = 0; i < NEQ + 1; i++) 
		cout << setw(8) << DiagonalAddress[i];
	cout << endl << endl;

	cout << "Matrix : " << endl;
	for (int I = 0; I < NEQ; I++)
	{
		for (int J = 0; J < NEQ; J++)
		{
			int i = I;
			int j = J;
			if (i > j)
			{
				int temp = i;
				i = j;
				j = temp;
			}

			cout << setiosflags(ios::scientific);

			int H = DiagonalAddress[j + 1] - DiagonalAddress[j];
			if (j - i - H >= 0) 
				cout << setw(15) << 0.0;
			else 
				cout << setw(15) << StiffnessMatrix[DiagonalAddress[j] + j - i - 1];
		}

		cout << endl;
	}

	cout << endl;

	cout << "U : " << endl;
	for (int I = 0; I < NEQ; I++)
		cout << setw(15) << Force[I];

	cout << endl;
}

#endif