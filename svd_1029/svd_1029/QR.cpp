#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "MatrixSVD.h"
#include <string>
#include "Givens.h"
#include "HouseHold.h"
#include <math.h>
#include "QR.h"

QR::QR(int n)
{
	P=new Matrix(n);
	Q=new Matrix(n);
}

void QR::Iterative(Vector &A,Vector &B,int i1,int i2)
{
	double d=(pow(A[i2-1],2)+pow(B[i2-1],2)-pow(A[i2],2)-pow(B[i2],2))/2;
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
