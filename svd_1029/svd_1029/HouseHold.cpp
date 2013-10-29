#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "MatrixSVD.h"
#include <string>
#include "Givens.h"
#include "HouseHold.h"
#include <math.h>


void HouseHold::HouseHolder(Vector v,int n)
{
	delta=v.modular();
	if(delta==0)
	{
		Matrix temp(n);
		T=&temp;
	}
	else
	{
	    v.set(0,v[0]-delta);
	    if(v.modular()==0)//以后要改
		{
		Matrix temp(n);
		T=&temp;
	    }
	    else 
	    v.Normalize();
	    v.Span(*T);
	    Matrix temp(n);
		T->NumProd(2);
	    temp.Minus(*T);
	}
}