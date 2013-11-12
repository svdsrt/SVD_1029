// svd_1111.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include "MatrixSVD.h"
#include "svd.h"

using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	Matrix  A(2,4);
	A(1,1) = 1;     A(1,2) = 3;     A(1,3) = 5;     A(1,4) = 7;
	A(2,1) = 2;     A(2,2) = 4;     A(2,3) = 6;     A(2,4) = 8;

	SVD  svd;
	svd.dec(A);

	Matrix  U = svd.getU();
	Matrix  V = svd.getV();
	Matrix  S = svd.getSM();

    cout << setiosflags(ios::fixed) << setprecision(4);
    cout << "Matrix--A: " ;
	printMatrix(A);
	cout<< endl;
	cout << "Matrix--U: "  ;
	printMatrix(U);
	cout<< endl;
	cout << "Vector--S: "  ;
	printMatrix(S);
	cout<< endl;
	cout << "Matrix--V: "  ;
	printMatrix(V);
	cout<< endl;
	cout << "Matrix--A - U * S * V^T:  ";
	printMatrix(minus(A,production(production(U,S),trT(V))));
	cout<< endl;
//         << A- U*multTr(S,V) << endl;

	cout << "The rank of A : " << svd.rank() << endl << endl;
	cout << "The condition number of A : " << svd.cond() << endl << endl;

	system("pause");
	
	return 0;
}

