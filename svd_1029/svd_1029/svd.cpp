#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <string>
#include "MatrixSVD.h"
#include "Givens.h"
#include "HouseHold.h"
#include "svd.h"

using namespace std;

void svd::BiDiag()
{
	static Matrix *P=new Matrix(m);
	static Matrix *H=new Matrix(n);
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
		
		tP.print();
		P->Copy(*U);
		U->DotProd(*P,tP);

		Vector tempV(n-num-1);
		tempV.HRow(A,num);
		tempV.print();
		HouseHold T2(tempV.N());
		T2.HouseHolder(tempV);
		B2[num+1]=T2.Delta();
		Matrix tH(n);
		for(int i=0;i<n-num-1;i++)
			for(int j=0;j<n-num-1;j++)
				tH.set(num+i,num+j,T2.TMatrix()->a(i,j));
		H->Copy(*V);
		V->DotProd(*H,tH);
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
	for(int i=0;i<m-n+2;i++)
		for(int j=0;j<m-n+2;j++)
			tP3.set(num+i,num+j,T3.TMatrix()->a(i,j));
	tP3.print();
	P->Copy(*U);
	U->DotProd(*P,tP3);
	

	HouseHold T4(m-n+2);
	T4.HouseHolder(tempU2);
	B1[n-1]=T4.Delta();
	Matrix tP4(m);
	for(int i=0;i<m-n+2;i++)
		for(int j=0;j<m-n+2;j++)
			tP4.set(num+i,num+j,T4.TMatrix()->a(i,j));
	tP4.print();
	P->Copy(*U);
	U->DotProd(*P,tP4);
}
/*
void svd::CheckConvergence()
{
	int p=0,q=0,flag=0;
	double e=0;
	cout<<"输入收敛系数：";
	cin>>e;
	
	for(int i=1;i<n;i++)
	{
		if(fabs(B2[i])<=e*(fabs(B1[i])+fabs(B1[i-1])))
		{
			B2[i]=0;
			flag=1;
		}
	}
	while(flag==0)
	{
		for(int j=n-1;j>0;j--)
		{
			if(B2[j]!=0)
			{
				flag=0;
				q=j;
				for(int k=q;k>0;k--)
					if(B2[k]==0)
						p=k;
				double max=B1.Max();

				int flag2=0;
				for(int i=p;p<q;p++)
					if(B1[i]<=e*max)
					{
						flag2++;
						B1[i]=0;
						double x=B2[i+1],y=B1[i+1];
						B2[i+1]=0;
						int l=1;

						while(l<q-i)
						{
							Givens g(y,x);
							B1[i+l]=g.getr();
							
							g.update(U,i);
							x=g.gets()*B2[i+l+1];
							B2[i+l+1]=g.getc()*B2[i+l+1];
							y=B1[i+l+1];
							l++;
						}
					}
					if(flag2==0)
					{
						QR or(p,q);
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
				if(fabs(B2[i])<=e*(fabs(B1[i])+fabs(B1[i-1])))
				{
					B2[i]=0;
					flag=1;
				}
			}

		}
	}
}

*/


void svd::BiPrint()
{
	B1.print();
	B2.print();
}

void svd::UVPrint()
{
	U->print();
	V->print();
}

/*
QR::QR(int n)
{
	P=new Matrix(n);
	Q=new Matrix(n);
}

void QR::Iterative(int i1,int i2)
{
	double d=(pow(s.B1[i2-1],2)+pow(B2[i2-1],2)-pow(A[i2],2)-pow(B[i2],2))/2;
	double u;
	if(d>0)
		u=pow(A[i2],2)-pow(B[i2],2)+d-sqrt(pow(d,2)+pow(A[i2-1],2)*pow(B[i2],2));
	else if(d==0)
		u=pow(A[i2],2)-pow(B[i2],2)+d;
	else
		u=pow(A[i2],2)-pow(B[i2],2)+d+sqrt(pow(d,2)+pow(A[i2-1],2)*pow(B[i2],2));
	double x=pow(A[i1],2)-u;
	double y=A[i1]*B[i1+1];
	int k=i1;
	int flag=0;
	while(flag==0)
	{
		Givens a=Givens(x,y);
		double x1=A[k+1];
		double x2=B[k+1];
		x=a.getc()*A[k]-a.gets()*x2;
		y=-a.gets()*x1;
		B.set(k+1,a.gets()*A[k]+a.getc()*x2);
		A.set(k+1,a.getc()*x1);
		a.update(Q,k+1);
		if(k>i1)
			B.set(k,a.getr());
		else
		{
			a=Givens(x,y);
			A.set(k,a.getr());
			a.update(P,k+1);
		}
		if(k<i2-1)
		{
			x1=A[k+1];x2=B[k+2];
			x=a.getc()*B[k+1]-a.gets()*x1;
			y=-a.gets()*x2;
			A.set(k+1,a.gets()*B[k+1]+a.getc()*x1);
			B.set(k+2,a.getc()*x2);
			k++;
		}
		else
		{
			x1=A[i2];x2=B[i2];
			B.set(i2,a.getc()*x2-a.gets()*x1);
			A.set(i2,a.gets()*x2+a.getc()*x1);
			flag=1;
		}
	}
}
*/