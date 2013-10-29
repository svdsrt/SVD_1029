// svd_1028.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include "MatrixSVD.h"
#include "Givens.h"
#include "HouseHold.h"

using namespace std;

void BiDiag(Matrix A,Vector B1,Vector B2,Matrix U, Matrix V,int m,int n);//二对角化

int _tmain(int argc, _TCHAR* argv[])
{
	//initalize
	int m,n;
	fstream infile("InputMatrix1.txt");
	infile>>m;
	infile>>n;
	Matrix U(m,m),V(n,n);
	Matrix A(m,n);
	Matrix tU(n,n),tV(n,n),temp(n,n);
	Matrix A("InputMatrix1.txt");

	Vector B1(n);
	Vector B2(n);


	BiDiag(A,B1,B2,U,V,m,n);

	B1.print();
	B2.print();


	int p=0,q=0,flag=0;
	double e=0;
	cout<<"璇疯緭鍏ヨ宸寖鍥?";
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
						temp.~Matrix();

						temp=U;
						U=tU;
						tU=temp;

						Unit(temp,n);
						for(int ii=p;ii<=q;ii++)
							for(int jj=p;jj<=q;jj++)
								temp.matrix[ii][jj]=svd1.P.matrix[ii-p][jj-p];
						tV.MatrixMultiply(V,temp);
						temp.~Matrix();

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
	}


	U.Print();
	cout<<endl;
	B1.Print();
	B2.Print();
	cout<<endl;
	V.Print();


	system("pause");
	return 0;

}

void BiDiag(Matrix A,Vector B1,Vector B2,Matrix U, Matrix V,int m,int n)
{
	Matrix P(m,m),H(n,n),tH(n,n);
	int num=0;
	for(num=0;num<n-2;num++)
	{
		Vector tempU(m-num);
		tempU.HCol(A,num);
		tempU.print();
		Matrix temp1(m-num,m-num);
		HouseHold T1();
		T1.HouseHolder(tempU,m-num);
		B1[num]=T1.delta;
		Matrix tP(m);
		for(int i=0;i<m-num;i++)
			for(int j=0;j<m-num;j++)
				tP.set(num,num,T1.T.matrix[i][j]);
		P.MatrixMultiply(U,tP);
		tP=U;
		U=P;
		P=tP;


		Vector tempV(n-num-1);
		tempV.VectorRowSpecial(A,num);
		tempV.Print();
		Matrix temp2(n-1-num,n-num-1);
		HouseHold T2(temp2);
		T2.HouseHolder(tempV,n-num-1);
		B2.vector[num+1]=T2.delta;
		Matrix tH(n,n);
		Unit(tH,n);
		for(int i=0;i<n-num-1;i++)
			for(int j=0;j<n-num-1;j++)
				tH.matrix[num][num]=T2.T.matrix[i][j];
		H.MatrixMultiply(V,tH);
		tH=V;
		V=H;
		H=tH;
	}

	Vector tempU1(m-n+2),tempU2(m-n+2);
	int num2=num;
	for(;num2<m;num2++)
	{
		tempU1.vector[num2-num]=A.matrix[num2][n-2];
		tempU2.vector[num2-num]=A.matrix[num2][n-1];
	}
	tempU1.Print();
	tempU2.Print();

}