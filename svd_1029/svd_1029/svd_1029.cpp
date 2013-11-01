// svd_1028.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include "MatrixSVD.h"
#include "Givens.h"
#include "HouseHold.h"
#include "svd.h"
#include <time.h>

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{

	double start,end;
	start=clock();
	
	svd test1;
	
	test1.BiDiag();

	test1.BiPrint();
	test1.UVPrint();

	end=clock();
	cout<<end-start<<endl;

	start=clock();
	
	test1.CheckConvergence();

	//test1.BiPrint();
	//test1.UVPrint();

	end=clock();
	cout<<end-start<<endl;

	system("pause");
	return 0;

}

