#ifndef GIVENS_H
#define GIVENS_H

class Givens
{
private:
	double c,s,r;

public:
	Givens();//c,s,r全部为0的构造函数
	Givens(double x, double y);//基于x,y计算出c,s,r
	Givens(double cvalue,double svalue,double rvalue);//直接用cvalue,svalue,rvalue赋值

	void set(double x,double y);//重新由x,y计算c,s,r
	void update(Matrix,int k);//Givens变换
};




#endif