/************************************************/
/*              FEMCPP 教学程序                 */
/*     清华大学航天航空学院计算动力学教研室     */
/*     Outputter.h  输出器头文件                */
/*     作者：宋言                               */
/*     最近修改：2017/06/04                     */
/********************************************；****/
#pragma once
#include <fstream>
#include <iostream>
#include <string>
using namespace std;
// 输出器文件，使用单例模式，可以在各出进行调用
class Outputter
{
private:
	static Outputter* _instance;
	ofstream OutputFile;      //数据文件输出流
	Outputter(string FileName);
public:
	static Outputter* Instance(string FileName = " ");
	void OutputLogo();
	void OutputNodeInfo();
	void OutputLoadInfo(int LoadCase);  //输出第LoadCase个工况的信息
	void OutputDisplacement();          //输出位移
};