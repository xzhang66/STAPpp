/***************************************************************/
/*  FEM++ : A C++ finite element method code for teaching      */
/*     Computational Dynamics Laboratory                       */
/*     School of Aerospace Engineering, Tsinghua University    */
/*                                                             */
/*     http://www.comdyn.cn/                                   */
/***************************************************************/

#include "Domain.h"
#include "Outputter.h"

#include <iostream>
#include <iomanip>
#include <ctime>

using namespace std;

//	Output current time and date
void Outputter::PrintTime(const struct tm * ptm, ostream& output)
{
	char *weekday[] = {"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};
	char *month[]   = {"January", "February", "March", "April", "May", "June", "July", "August", 
		             "September", "October", "November", "December"};

	output << endl << "        ";
	output << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec << " ";
	output << month[ptm->tm_mon+1] << " " << ptm->tm_mday << ", " << ptm->tm_year + 1900 << " " 
		   << weekday[ptm->tm_wday] << endl << endl;
}

Outputter* Outputter::_instance = NULL;

//	Constructor
Outputter::Outputter(string FileName)
{
	OutputFile.open(FileName);

	if (!OutputFile) 
	{
		cout << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
Outputter* Outputter::Instance(string FileName)
{
	if (!_instance) _instance = new Outputter(FileName);
	return _instance;
}

//	Print program logo
void Outputter::PrintLogo(ostream& output)
{
	output << "***********************************************************" << endl;
	output << "* xxxxxx  xxxxxx  xxx       xx        x            x      *" << endl;
	output << "* xx      xx      xxxx     xxx        x            x      *" << endl;
	output << "* xx      xx      xxxx     x x        x            x      *" << endl;
	output << "* xx      xx      xx x    xx x        x            x      *" << endl;
	output << "* xx      xx      xx xx   xx x        x            x      *" << endl;
	output << "* xxxxxx  xxxxxx  xx xx   x  x   xxxxxxxxxxx  xxxxxxxxxxx *" << endl;
	output << "* xx      xx      xx  xx xx  x        x            x      *" << endl;
	output << "* xx      xx      xx  xx x   x        x            x      *" << endl;
	output << "* xx      xx      xx   xxx   x        x            x      *" << endl;
	output << "* xx      xx      xx   xxx   x        x            x      *" << endl;
	output << "* xx      xxxxxxx xx    x    x        x            x      *" << endl;
	output << "***********************************************************" << endl << endl;
}
//	Print program logo
void Outputter::OutputHeading()
{
	PrintLogo(cout);
	PrintLogo(OutputFile);

	Domain* FEMData = Domain::Instance();

	cout << "TITLE : " << FEMData->GetTitle() << endl;
	OutputFile << "TITLE : " << FEMData->GetTitle() << endl;

	time_t now;
	struct tm *local = new struct tm;
	now = time(NULL);

#if defined(WIN32) 
	localtime_s(local, &now);	//  localtime_s is used only in windows
#else
    localtime_r(&now, local);	//  localtime_r is used only in unix
#endif

	PrintTime(local, cout);
	PrintTime(local, OutputFile);
}

//	Print nodal data
void Outputter::OutputNodeInfo()
{
	Domain* FEMData = Domain::Instance();

	Node* NodeList = FEMData->GetNodeList();

//	Number of lines printed in each page
	int Page = 30;

	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);

	cout << "*********************  N O D E ****************************" << endl;
	OutputFile << "*********************  N O D E ****************************" << endl;

	for (int i = 0; i < FEMData->GetNUMNP(); i++)
	{
		if (i % Page == 0)
		{
			cout << setw(10) << "Node" << setw(15) << "X" << setw(15) << "Y" << setw(15) << "Z" << endl;
			OutputFile << setw(10) << "Node" << setw(15) << "X" << setw(15) << "Y" << setw(15) << "Z" << endl;
		}

		cout << setw(10) << i + 1 << setw(15) << NodeList[i].XYZ[0] << setw(15) << NodeList[i].XYZ[1] 
		     << setw(15) << NodeList[i].XYZ[2] << endl;
		OutputFile << setw(10) << i + 1 << setw(15) << NodeList[i].XYZ[0] << setw(15) << NodeList[i].XYZ[1] 
		           << setw(15) << NodeList[i].XYZ[2] << endl;
	}
	cout << endl;
	OutputFile << endl;
}

//	Print load data for load case LoadCase
void Outputter::OutputLoadInfo(int LoadCase)
{
	Domain* FEMData = Domain::Instance();

	unsigned int* NLOAD = FEMData->GetNLOAD();

	if (LoadCase > FEMData->GetNLCASE()) 
		return;

	LoadData* Load = FEMData->GetLoadList()[LoadCase - 1];

	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);

	cout << "******************** L O A D C A S E " << LoadCase <<  " ********************" << endl;
	OutputFile << "******************** L O A D C A S E " << LoadCase <<  " ********************" << endl;

	cout << setw(10) << "Load case" << setw(14) << "Node" << setw(12) << "DOF" << setw(17) << "Load" << endl;
	OutputFile << setw(10) << "Load case" << setw(14) << "Node" << setw(12) << "DOF" << setw(17) << "Load" << endl;

	for (int i = 0; i < NLOAD[LoadCase - 1]; i++)
	{
		cout << setw(10) << i + 1 << setw(14) << Load[i].node << setw(12) << Load[i].dof 
			 << setw(17) << Load[i].load << endl;
		OutputFile << setw(10) << i + 1 << setw(14) << Load[i].node << setw(12) << Load[i].dof 
			       << setw(17) << Load[i].load << endl;
	}
	cout << endl;
	OutputFile << endl;
}

//	Print nodal displacement
void Outputter::OutputNodalDisplacement()
{
	Domain* FEMData = Domain::Instance();

	int Page = 30;

	Node* NodeList = FEMData->GetNodeList();

	double* Displacement = FEMData->GetDisplacement();

	cout << setiosflags(ios::scientific);
	OutputFile << setiosflags(ios::scientific);

	cout << "************* D I S P L A C E l M E N T *****************" << endl;
	OutputFile << "************* D I S P L A C E l M E N T *****************" << endl;

	for (int i = 0; i < FEMData->GetNUMNP(); i++)
	{
		if (i % Page == 0)
		{
			cout << setw(12) << "Node" << setw(14) << "X" << setw(14) << "Y" << setw(14) << "Z" << endl;
			OutputFile << setw(12) << "Node" << setw(14) << "X" << setw(14) << "Y" << setw(14) << "Z" << endl;
		}

		cout << setw(12) << i + 1;
		OutputFile << setw(12) << i + 1;

		for (int j = 0; j < 3; j++)
		{
			if (NodeList[i].bcode[j] == 0)
			{
				cout << setw(14) << 0.0;
				OutputFile << setw(14) << 0.0;
			}
			else
			{
				cout << setw(14) << Displacement[NodeList[i].bcode[j] - 1];
				OutputFile << setw(14) << Displacement[NodeList[i].bcode[j] - 1];
			}
		}

		cout << endl;
		OutputFile << endl;
	}
}


#ifdef _DEBUG_

//	Print column heights for debuging
void Outputter::PrintColumnHeights()
{
	cout << "************************** Column Heights ****************************" << endl;
	OutputFile << "************************** Column Heights ****************************" << endl;

	Domain* FEMData = Domain::Instance();

	int NEQ = FEMData->GetNEQ();
	unsigned int* ColumnHeights = FEMData->GetColumnHeights();

	for (int col = 0; col < NEQ; col++)
	{
		if (col+1 % 10 == 0)
		{
			cout << endl;
			OutputFile << endl;
		}

		cout << setw(8) << ColumnHeights[col];
		OutputFile << setw(8) << ColumnHeights[col];
	}
	cout << endl << endl;
	OutputFile << endl << endl;
}

//	Print address of diagonal elements for debuging
void Outputter::PrintDiagonalAddress()
{
	cout << "******************** Address of Diagonal Element *********************" << endl;
	OutputFile << "******************** Address of Diagonal Element *********************" << endl;

	Domain* FEMData = Domain::Instance();

	int NEQ = FEMData->GetNEQ();
	unsigned int* DiagonalAddress = FEMData->GetDiagonalAddress();

	for (int col = 0; col <= NEQ; col++)
	{
		if (col+1 % 10 == 0)
		{
			cout << endl;
			OutputFile << endl;
		}

		cout << setw(8) << DiagonalAddress[col];
		OutputFile << setw(8) << DiagonalAddress[col];
	}

	cout << endl << endl;
	OutputFile << endl << endl;
}

//	Print banded and full stiffness matrix for debuging
void Outputter::PrintStiffnessMatrix()
{
	cout << "********************** Banded stiffness matrix ***********************" << endl;
	OutputFile << "********************** Banded stiffness matrix ***********************" << endl;

	Domain* FEMData = Domain::Instance();

	int NEQ = FEMData->GetNEQ();
	unsigned int* DiagonalAddress = FEMData->GetDiagonalAddress();
	double* StiffnessMatrix = FEMData->GetStiffnessMatrix();

	cout << setiosflags(ios::scientific) <<setprecision(5);
	OutputFile << setiosflags(ios::scientific) <<setprecision(5);

	for (int i = 0; i < DiagonalAddress[NEQ]-1; i++)
	{
		if ((i+1) % 6 == 0)
		{
			cout << endl;
			OutputFile << endl;
		}

		cout << setw(14) << StiffnessMatrix[i];
		OutputFile << setw(14) << StiffnessMatrix[i];
	}

	cout << endl << endl;
	OutputFile << endl << endl;

	cout << "*********************** Full stiffness matrix ************************" << endl;
	OutputFile << "*********************** Full stiffness matrix ************************" << endl;

	for (int I = 0; I < NEQ; I++)
	{
		for (int J = 0; J < NEQ; J++)
		{
			int i = I;
			int j = J;
			if (i > j)
				swap(i,j);

			int H = DiagonalAddress[j + 1] - DiagonalAddress[j];
			if (j - i - H >= 0) 
			{
				cout << setw(14) << 0.0;
				OutputFile << setw(14) << 0.0;
			}
			else
			{
				cout << setw(14) << StiffnessMatrix[DiagonalAddress[j] + j - i - 1];
				OutputFile << setw(14) << StiffnessMatrix[DiagonalAddress[j] + j - i - 1];
			}
		}

		cout << endl;
		OutputFile << endl;
	}

	cout << endl;
	OutputFile << endl;
}

//	Print displacement vector for debuging
void Outputter::PrintDisplacement(int loadcase)
{
	cout << "********************** Displacement vector ***********************" << endl;
	OutputFile << "********************** Displacement vector ***********************" << endl;

	Domain* FEMData = Domain::Instance();

	int NEQ = FEMData->GetNEQ();
	double* Force = FEMData->GetForce();

	cout << "  Load case = " << loadcase << endl;
	OutputFile << "  Load case = " << loadcase << endl;

	cout << setiosflags(ios::scientific) <<setprecision(5);
	OutputFile << setiosflags(ios::scientific) <<setprecision(5);

	for (int i = 0; i < NEQ; i++)
	{
		if ((i+1) % 6 == 0)
		{
			cout << endl;
			OutputFile << endl;
		}

		cout << setw(14) << Force[i];
		OutputFile << setw(14) << Force[i];
	}

	cout << endl << endl;
	OutputFile << endl << endl;
}

#endif