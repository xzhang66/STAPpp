/***************************************************************/
/*  FEM++ £ºA C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "FEM.h"
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

Node::Node(double X = 0, double Y = 0, double Z = 0) 
{
	XYZ[0] = X;		// Coordinates of the node
	XYZ[1] = Y;
	XYZ[2] = Z;

	bcode[0] = 0;	// Boundary codes
	bcode[1] = 0;
	bcode[2] = 0;

	EquationNo[0] = 0;	// Equation numbers 
	EquationNo[1] = 0;
	EquationNo[2] = 0;
};

//	Generate the location matrix of the element
void Element :: GetLocationMatrix()
{
	for (int N = 0; N < NEN; N++)
		for (int D = 0; D < Node::NDF; D++)
			LocationMatrix.push_back(nodes[N]->EquationNo[D]);
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
	ElementList = NULL;
	
	NUME = NULL;
	NUMMAT = NULL;
	MaterialList = NULL;
	
	NLCASE = 0;
	NLOAD = NULL;
	LoadList = NULL;
	
	NEQ = 0;
	StiffnessMatrix = NULL;
	Displacement = NULL;
	Force = NULL; 
}

//	Return pointer to the instance of the Domain class
Domain* Domain::Instance()
{
	if (!_instance) 
		_instance = new Domain;
	
	return _instance;
}

//	Initialize the class data member by reading input data file
bool Domain::Initial(FileReader* Reader)
{
	return Reader->ReadData(this);
}

// Constructor
FileReader::FileReader(string Filename)
{
	Input.open(Filename);

	if (!Input) 
	{
		cout << "*** Error *** File " << Filename << " does not exist !" << endl;
		exit(3);
	}
}

// Read data from input data file and initialize the data member of FEMData
bool FileReader::ReadData(Domain* FEMData)
{
	int N, LL, NL;

//	Read the heading line
	Input >> FEMData->Title;

//	Read the control line
	Input >> FEMData->NUMNP >> FEMData->NUMEG >> FEMData->NLCASE >> FEMData->MODEX;

//	Read nodal point data lines
	FEMData->NodeList = new Node[FEMData->NUMNP];
	Node* NodeList = FEMData->NodeList;
	
	for (int i = 0; i < FEMData->NUMNP; i++)
	{
		Input >> N;
		if (N != i + 1) 
		{
			cout << "*** Error *** Nodes must be inputted in order !" << endl 
				 << "   Expected node number : " << i + 1 << endl
				 << "   Provided node number : " << N << endl;

			return false;
		}

		Input >> NodeList[i].bcode[0] >> NodeList[i].bcode[1] >> NodeList[i].bcode[2]
			  >> NodeList[i].XYZ[0] >> NodeList[i].XYZ[1] >> NodeList[i].XYZ[2];
	}

//	Read load data lines
	FEMData->LoadList = new LoadData*[FEMData->NLCASE];
	FEMData->NLOAD = new unsigned int[FEMData->NLCASE];

	LoadData** LoadList = FEMData->LoadList;
	unsigned int* NLOAD = FEMData->NLOAD;
	
	for (int i = 0; i < FEMData->NLCASE; i++)
	{
		Input >> LL >> NL;

		if (LL != i + 1) 
		{
			cout << "*** Error *** Load case must be inputted in order !" << endl 
				 << "   Expected load case : " << i + 1 << endl
				 << "   Provided load case : " << LL << endl;

			return false;
		}

		NLOAD[i] = NL;

		LoadList[i] = new LoadData[NLOAD[i]];
		for (int j = 0; j < NLOAD[i]; j++)
		{
			Input >>LoadList[i][j].node >> LoadList[i][j].dof >> LoadList[i][j].load;
		}
	}

//	Read element group control line
	FEMData->NUME = new unsigned int[FEMData->NUMEG];
	unsigned int* NUME = FEMData->NUME;
	
	FEMData->ElementList = new Element*[FEMData->NUMEG];

	FEMData->NUMMAT = new unsigned int[FEMData->NUMEG];
	unsigned int* NUMMAT = FEMData->NUMMAT;
	FEMData->MaterialList = new Material*[FEMData->NUMEG];

	unsigned int ElementType;
	for (int i = 0; i < FEMData->NUMEG; i++)
	{
		Input >> ElementType >> NUME[i] >> NUMMAT[i];

//		Will try to move to element class
		switch (ElementType)
		{
		case 1:	// Bar element
		{
//			Read material/section property lines
			BarMaterial* MaterialGroup = new BarMaterial[NUMMAT[i]];
			FEMData->MaterialList[i] = MaterialGroup;

			for (int j = 0; j < NUMMAT[i]; j++)
			{
				Input >> N;

				if (N != j + 1)
				{
					cout << "*** Error *** Material sets must be inputted in order !" << endl 
						 << "   Expected set : " << j + 1 << endl
						 << "   Provided set : " << N << endl;

					return false;
				}
				
				Input >> MaterialGroup[j].E >> MaterialGroup[j].Area;
			}

//			Read element data lines
			Bar* ElementList = new Bar[NUME[i]];
			FEMData->ElementList[i] = ElementList;

			for (int j = 0; j < NUME[i]; j++)
			{
				Input >> N;

				if (N != j + 1)
				{
					cout << "*** Error *** Elements must be inputted in order !" << endl 
						 << "   Expected element : " << j + 1 << endl
						 << "   Provided element : " << N << endl;

					return false;
				}

				int MSet;	// Material property set number
				int N1, N2;	// Node number of the left and right node

				Input >> N1 >> N2 >> MSet;
				ElementList[j].ElementMaterial = &MaterialGroup[MSet - 1];
				ElementList[j].nodes[0] = &NodeList[N1 - 1];
				ElementList[j].nodes[1] = &NodeList[N2 - 1];
			}
			break;
		}

		default:
			return false;
		}
	}
	return true;
}

//	Calculate global equation numbers corresponding to every degree of freedom of each node
void Domain::EquationNumber()
{
	NEQ = 0;
	for (int i = 0; i < NUMNP; i++)
	{
		for (int j = 0; j < Node::NDF; j++)
		{
			if (NodeList[i].bcode[j]) 
				NodeList[i].EquationNo[j] = 0;
			else
			{
				NEQ++;
				NodeList[i].EquationNo[j] = NEQ;
			}
		}
	}
}

//	Allocate storage for the one dimensional array storing the global stiffness matrix,
//	and generate the address of diagonal elements
void Domain::AllocateStiffnessMatrix()
{
	Displacement = new double[NEQ];
	Force = new double[NEQ];
	DiagonalAddress = new unsigned int[NEQ + 1];

	clear(DiagonalAddress, NEQ + 1);

//	Calculate column heights (stored in DiagonalAddress)
	for (int EG = 0; EG < NUMEG; EG++)
	{
		for (int EN = 0; EN < NUME[EG]; EN++)
			ElementList[EG][EN].ColumnHeight(DiagonalAddress);
	}

//	Calculate the address of diagonal elements (Overwrite DiagonalAddress)
//	M(i+1) = M(i) + H(i) + 1
	DiagonalAddress[0] = 1;
	for (int C = 1; C <= NEQ; C++)
		DiagonalAddress[C] = DiagonalAddress[C - 1] + DiagonalAddress[C] + 1;

	StiffnessMatrix = new double[DiagonalAddress[NEQ] - 1];

	clear(StiffnessMatrix, DiagonalAddress[NEQ] - 1);
}

//	Assemble the banded gloabl stiffness matrix
void Domain::AssembleStiffnessMatrix()
{
	for (int EG = 0; EG < NUMEG; EG++)
	{
		unsigned int size = ElementList[EG][0].SizeOfStiffnessMatrix();
		double* Matrix = new double[size];

		for (int E = 0; E < NUME[EG]; E++)
			ElementList[EG][E].assembly(Matrix);

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

	for (int FN = 0; FN < NLOAD[LoadCase - 1]; FN++)
	{
		int dof = NodeList[Load[FN].node - 1].EquationNo[Load[FN].dof - 1];
		Force[dof - 1] += Load[FN].load;
	}

	return true;
}

#ifdef _DEBUG_

//	Print debug information
void Domain::Info()
{
	cout << endl << endl;
	cout << "**************** Debug Infomation ******************" << endl;
	cout << "StiffnessMatrix : " << endl;
	for (int i = 0; i < DiagonalAddress[NEQ] - 1; i++) 
		cout << StiffnessMatrix[i] << " ";
	cout << endl << endl;

	cout << "Address : " << endl;
	for (int i = 0; i < NEQ + 1; i++) 
		cout << DiagonalAddress[i] << " ";
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
				cout << setw(15) << StiffnessMatrix[DiagonalAddress[j] + j - i - 1] << "  ";
		}

		cout << endl;
	}

	cout << endl;

	cout << "U : " << endl;
	for (int I = 0; I < NEQ; I++)
		cout << Displacement[I] << " ";

	cout << endl;
}

#endif