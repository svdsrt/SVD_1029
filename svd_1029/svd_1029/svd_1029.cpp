// svd_1028.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include "MatrixSVD.h"
#include "Givens.h"
#include "HouseHold.h"
#include "svd.h"

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{

	svd test1;
	
	test1.BiDiag();

	test1.BiPrint();
	test1.UVPrint();
	
	test1.CheckConvergence();

	test1.BiPrint();
	test1.UVPrint();

	system("pause");
	return 0;

}

