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
	Matrix *tA1=new Matrix(m,n);
	Matrix *tA2=new Matrix(m,n);
	tA1->Copy(A);
	tA2->Copy(A);
	int num=0;
	for(num=0;num<n-1;num++)
	{
		Vector tempU(m-num);
		tempU.HCol(*tA1,num);
		tempU.print();
		HouseHold T1(m-num);
		T1.HouseHolder(tempU);
		B1[num]=T1.Delta();
		Matrix tP(m);
		for(int i=0;i<m-num;i++)
			for(int j=0;j<m-num;j++)
				tP.set(num+i,num+j,T1.TMatrix()->a(i,j));
		P->Copy(*U);
		U->DotProd(*P,tP);
		tP.print();
		tA2->DotProd(tP,*tA1);
		tA2->scanzero();

		Vector tempV(n-num-1);
		tempV.HRow(*tA2,num);
		HouseHold T2(tempV.N());
		T2.HouseHolder(tempV);
		B2[num+1]=T2.Delta();
		Matrix tH(n);
		for(int i=0;i<n-num-1;i++)
			for(int j=0;j<n-num-1;j++)
				tH.set(num+1+i,num+1+j,T2.TMatrix()->a(i,j));
		H->Copy(*V);
		V->DotProd(*H,tH);
		tH.print();
		tA1->DotProd(*tA2,tH);
		tA1->scanzero();

	}

	Vector tempU1(m-n+1);
	tempU1.HCol(*tA1,n-1);
	HouseHold T3(m-n+1);
	T3.HouseHolder(tempU1);
	B1[n-1]=T3.Delta();
	Matrix tP3(m);
	for(int i=0;i<m-n+1;i++)
		for(int j=0;j<m-n+1;j++)
			tP3.set(num+i,num+j,T3.TMatrix()->a(i,j));
	P->Copy(*U);
	U->DotProd(*P,tP3);


	Matrix result1(m,n);
	Matrix result2(m,n);
	Matrix S(m,n);
	for(int mi=0;mi<n;mi++)
		S.set(mi,mi,B1[mi]);
	for(int mj=1;mj<n;mj++)
		S.set(mj-1,mj,B2[mj]);
	S.print();
	Matrix U2(m,m);
	U->Trans(U2);
	result1.DotProd(U2,S);
	result2.DotProd(result1,*V);
	result2.print();
}

void svd::function1()
{
	for(int j=1;j<n;j++)
	{
		if(fabs(B2[j])<=e*(fabs(B1[j])+fabs(B1[j-1])))
		{
			B2[j]=0;
			flag=flag*1;
		}
		else flag=flag*0;
	}
	int flag3=0;
	for(int j=0;j<B2.N();j++)
		if(B2[j]!=0)
		{
			flag3=1;
		}
		if(flag3==0) 
		{
			BiPrint();
			UVPrint();
			cout<<"cao"<<endl;
			system("pause");
			exit(0);
		}
		for(int j=n-1;j>0;j--)
		{
			if(B2[j]!=0)
			{
				flag=0;
				q=j;
				for(int k=q;k>0;k--)
					if(B2[k]==0)
						p=k;
			}
		}
	int flag2=0;
	for(i=p;i<q;i++)
	{
		if(B1[i]<=e*B1.Max())
			flag2++;
	}
	if(flag2>0)
	{
		B1[i]=0;
		double x=B2[i+1],y=B1[i+1];
		B2[i+1]=0;
		int l=1;
		function4();
	}
	else
	{
		QR or(q-p+1);
		or.Iterative(p,q,B1,B2);
		Matrix temp1(n);
		for(int ii=p;ii<=q;ii++)
			for(int jj=p;jj<=q;jj++)
				temp1.set(ii,jj,or.P->a(ii-p,jj-p));
		tU.Copy(*U);
		U->DotProd(tU,temp1);

		Matrix temp2(n);
		for(int ii=p;ii<=q;ii++)
			for(int jj=p;jj<=q;jj++)
				temp2.set(ii,jj,or.Q->a(ii-p,jj-p));
		tV.Copy(*V);
		V->DotProd(tV,temp);
		function1();
	}

}

void svd::function4()
{
	Givens g(y,x);
	B1[i+l]=g.getr();
	g.update(U,i,i+l);
	function5();
}

void svd:: function5()
{
	if(l<q-i)
	{
		x=g.gets()*B2[i+l+1];
		B2[i+l+1]=g.getc()*B2[i+l+1];
		y=B1[i+l+1];
		l++;
		function4();
	}
	else
		function1();
}





void svd::CheckConvergence()
{
	int p=0,q=0,flag=0;
	double e=0;
	cout<<"输入收敛系数：";
	cin>>e;
	function1();	
}



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


QR::QR(int n)
{
	P=new Matrix(n);
	Q=new Matrix(n);
}

void QR::Iterative(int i1,int i2,Vector &B1,Vector &B2)
{
	double d=(pow(B1[i2-1],2)+pow(B2[i2-1],2)-pow(B1[i2],2)-pow(B2[i2],2))/2;
	double u;
	if(d>0)
		u=pow(B1[i2],2)-pow(B2[i2],2)+d-sqrt(pow(d,2)+pow(B1[i2-1],2)*pow(B2[i2],2));
	else if(d==0)
		u=pow(B1[i2],2)-pow(B2[i2],2)+d;
	else
		u=pow(B1[i2],2)-pow(B2[i2],2)+d+sqrt(pow(d,2)+pow(B1[i2-1],2)*pow(B2[i2],2));
	double x=pow(B1[i1],2)-u;
	double y=B1[i1]*B2[i1+1];
	int k=i1;
	int flag=0;
	while(flag==0)
	{
		Givens a;
		a.set(x,y);
		double x1=B1[k+1];
		double x2=B2[k+1];
		x=a.getc()*B1[k]-a.gets()*x2;
		y=-a.gets()*x1;
		B2.set(k+1,a.gets()*B1[k]+a.getc()*x2);
		B1.set(k+1,a.getc()*x1);
		a.update(Q,k+1);
		if(k>i1)
			B2.set(k,a.getr());
		else
		{
			a.set(x,y);
			B1.set(k,a.getr());
			a.update(P,k+1);
		}
		if(k<i2-1)
		{
			x1=B1[k+1];x2=B2[k+2];
			x=a.getc()*B2[k+1]-a.gets()*x1;
			y=-a.gets()*x2;
			B1.set(k+1,a.gets()*B2[k+1]+a.getc()*x1);
			B2.set(k+2,a.getc()*x2);
			k++;
		}
		else
		{
			x1=B1[i2];x2=B2[i2];
			B2.set(i2,a.getc()*x2-a.gets()*x1);
			B1.set(i2,a.gets()*x2+a.getc()*x1);
			flag=1;
		}
	}
}