#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "MatrixSVD.h"
#include <string>
#include "Givens.h"
#include "HouseHold.h"
#include <math.h>

HouseHold::HouseHold()
{
	T=NULL;
	delta=0;
}

HouseHold::HouseHold(int n)
{
	Matrix M(n,n);
	T=&M;
	delta=0;
}

double HouseHold::Delta()
{
	return delta;
}

void HouseHold::HouseHolder(Vector v)
{
	delta=v.modular();
	if(delta==0)
	{
		Matrix temp(v.N());
		T=&temp;
	}
	else
	{
	    v.set(0,v[0]-delta);
	    if(v.modular()==0)//以后要改
		{
		Matrix temp(v.N());
		T=&temp;
	    }
	    else 
	    v.Normalize();
	    v.Span(*T);
	    Matrix temp(v.N());
		T->NumProd(2);
	    temp.Minus(*T);
	}
}