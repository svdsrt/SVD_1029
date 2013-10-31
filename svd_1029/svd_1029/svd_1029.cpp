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

	/*int p=0,q=0,flag=0;
	double e=0;
	cout<<"";
	cin>>e;
	SVD svd1;
	svd1.B1=B1;
	svd1.B2=B2;


	for(int i=1;i<n;i++)
	{
		if(fabs(B2.vector[i])<=e*(fabs(B1.vector[i])+fabs(B1.vector[i-1])))
		{
			B2.vector[i]=0;
			flag=1;
		}
	}
	while(flag==0)
	{
		for(int j=n-1;j>0;j--)
		{
			if(B2.vector[j]!=0)
			{
				flag=0;
				q=j;
				for(int k=q;k>0;k--)
					if(B2.vector[k]==0)
						p=k;
				double max=B1.max();

				int flag2=0;
				for(int i=p;p<q;p++)
					if(B1.vector[i]<=e*max)
					{
						flag2++;
						B1.vector[i]=0;
						double x=B2.vector[i+1],y=B1.vector[i+1];
						B2.vector[i+1]=0;
						int l=1;

						while(l<q-i)
						{
							Givens g(y,x);
							B1.vector[i+l]=g.r;
							
							g.Update(U,i);
							x=g.s*B2.vector[i+l+1];
							B2.vector[i+l+1]=g.c*B2.vector[i+l+1];
							y=B1.vector[i+l+1];
							l++;
						}
					}
					if(flag2==0)
					{
						svd1.QR(p,q);
						Unit(temp,n);
						for(int ii=p;ii<=q;ii++)
							for(int jj=p;jj<=q;jj++)
								temp.matrix[ii][jj]=svd1.P.matrix[ii-p][jj-p];
						tU.MatrixMultiply(U,temp);

						temp=U;
						U=tU;
						tU=temp;

						Unit(temp,n);
						for(int ii=p;ii<=q;ii++)
							for(int jj=p;jj<=q;jj++)
								temp.matrix[ii][jj]=svd1.P.matrix[ii-p][jj-p];
						tV.MatrixMultiply(V,temp);

						temp=V;
						V=tV;
						tV=temp;
					}
			}
			for(int i=1;i<n;i++)
			{
				if(fabs(B2.vector[i])<=e*(fabs(B1.vector[i])+fabs(B1.vector[i-1])))
				{
					B2.vector[i]=0;
					flag=1;
				}
			}

		}
	}*/


	system("pause");
	return 0;

}

