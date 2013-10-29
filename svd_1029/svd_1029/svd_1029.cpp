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
	Matrix *U=new Matrix(m,m);
	Matrix *V=new Matrix(n,n);
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


	U.print();
	cout<<endl;
	B1.print();
	B2.print();
	cout<<endl;
	V.print();


	system("pause");
	return 0;

}

void BiDiag(Matrix A,Vector B1,Vector B2,Matrix *U, Matrix *V,int m,int n)
{
	Matrix *P=new Matrix(m,m);
	Matrix *H=new Matrix(n,n);
	Matrix *swap=NULL;
	int num=0;
	for(num=0;num<n-2;num++)
	{
		Vector tempU(m-num);
		tempU.HCol(A,num);
		tempU.print();
		HouseHold T1(m-num);
		T1.HouseHolder(tempU);
		B1[num]=T1.Delta();
		Matrix tP(m);
		for(int i=0;i<m-num;i++)
			for(int j=0;j<m-num;j++)
				tP.set(num+i,num+j,T1.TMatrix()->a(i,j));
	    U->DotProd(*P,tP);
		swap=U;
		U=P;
		P=swap;


		Vector tempV(n-num-1);
		tempV.HCol(A,num);
		tempV.print();
		HouseHold T2(n-1-num);
		T2.HouseHolder(tempV);
		B2[num+1]=T2.Delta();
		Matrix tH(n);
		for(int i=0;i<n-num-1;i++)
			for(int j=0;j<n-num-1;j++)
				tH.set(num+i,num+j,T2.TMatrix()->a(i,j));
		V->DotProd(*H,tH);
		swap=V;
		V=H;
		H=swap;
	}

	Vector tempU1(m-n+2),tempU2(m-n+2);
	int num2=num;
	for(;num2<m;num2++)
	{
		tempU1.set(num2-num,A.a(num2,n-2));
		tempU2.set(num2-num,A.a(num2,n-1));
	}
	tempU1.print();
	tempU2.print();

	HouseHold T3(m-n+2);
	T3.HouseHolder(tempU1);
		B1[n-2]=T3.Delta();
		Matrix tP3(m);
		for(int i=0;i<n-num-1;i++)
			for(int j=0;j<n-num-1;j++)
				tP3.set(num+i,num+j,T3.TMatrix()->a(i,j));
		U->DotProd(*P,tP3);
		//还要吗？
		swap=U;
		U=P;
		P=swap;

	HouseHold T3(m-n+2);
	T3.HouseHolder(tempU2);
		B1[n-1]=T3.Delta();
		Matrix tP3(m);
		for(int i=0;i<n-num-1;i++)
			for(int j=0;j<n-num-1;j++)
				tP3.set(num+i,num+j,T3.TMatrix()->a(i,j));
		U->DotProd(*P,tP3);
		swap=U;
		U=P;
		P=swap;
}