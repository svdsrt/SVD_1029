//该代码用于测试类

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include "MatrixSVD.h"
#include <string>

using namespace std;

int main()
{
	std::string s1="inputMatrix1.txt";
	std::string s2="inputMatrix2.txt";

	Matrix A(s1);
	Matrix B(s2);

	A.print();
	B.print();

	A.Add(B);
	A.print();

	A.NumProd(2);
	A.print();

	Matrix C(A.N(),A.M());
	A.Trans(C);
	C.print();

	Matrix D(A.M(),C.N());
	D.DotProd(A,C);
	D.print();

	Vector r(D.N());
	for(int i=0;i<D.M();i++)
	{
		r.Row(D,i);
		r.print();
	}

	Vector c(D.M());
	for(int j=0;j<D.N();j++)
	{
		c.Col(D,j);
		c.print();
	}



	system("pause");
	return 0;


}