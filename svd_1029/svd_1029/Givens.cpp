#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "MatrixSVD.h"
#include <string>
#include "Givens.h"
#include <math.h>

using namespace std;

Givens::Givens()
{
	c=0;
	s=0;
	r=0;
}

Givens::Givens(double x, double y)
{
	double tao;
	if(y==0)
	{
		c=1;
		s=1;
	}
	else
		if (fabs(y)>fabs(x))
		{
			tao=-x/y;
			s=sqrt(1+tao*tao);
			r=-y*s;
			c=tao/s;
		}
		else
		{
			tao=-y/x;
			c=sqrt(1+tao*tao);
			r=x*c;
			c=tao/c;
		}
}
	
Givens::Givens(double cvalue,double svalue,double rvalue)
{
	c=cvalue;
	s=svalue;
	r=rvalue;
}

void Givens::set(double x,double y)
{
	double tao;
	if(y==0)
	{
		c=1;
		s=1;
	}
	else
		if (fabs(y)>fabs(x))
		{
			tao=-x/y;
			s=sqrt(1+tao*tao);
			r=-y*s;
			c=tao/s;
		}
		else
		{
			tao=-y/x;
			c=sqrt(1+tao*tao);
			r=x*c;
			c=tao/c;
		}
}

void Givens::update(Matrix X,int k)
{
	int i=0;
	for(i=0;i<X.N();i++)
	{
		double t=X.a(i,k-1);
		X.set(i,k-1,c*t-s*X.a(i,k));
		X.set(i,k,s*t+c*X.a(i,k));
	}
}