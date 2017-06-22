/************************************************/
/*              FEMCPP 教学程序                 */
/*     清华大学航天航空学院计算动力学教研室     */
/*     FEMMain.cpp  FEMCPP程序入口              */
/*     作者：宋言                               */
/*     最近修改：2017/06/12                     */
/************************************************/
#include "FEM.h"
#include "Truss.h"
#include "Outputter.h"
#include <string>
#include <iostream>
using namespace std;
int main()
{
	string InFile = "C:/learn/FEMCPP/src/Input.data";
	string OutFile = "C:/learn/FEMCPP/src/Output.out";
	FileReader* Reader = new FileReader(InFile);
	FEM* FEMData = FEM::Instance();
	bool ifsucessed = FEMData->Initial(Reader);
	if (!ifsucessed)
	{
		cout << "Input Failed!" << endl;
		exit(1);
	}
	Outputter* Output = Outputter::Instance(OutFile);
	Output->OutputLogo();
	Output->OutputNodeInfo();
	FEMData->GenerateFreedom();
	FEMData->AllocateStiffnessMatrix();
	FEMData->AssemblyStiffnessMatrix();
	LDLTSolver* S = new LDLTSolver(FEMData);
	S->Solve();
	return 0;
}