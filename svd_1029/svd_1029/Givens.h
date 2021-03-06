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
	void update(Matrix *X,int k);//Givens变换
	void update(Matrix *X,int k,int l);//Givens变换
	double getc();//输出c
	double gets();//输出s
	double getr();//输出r
};




#endif