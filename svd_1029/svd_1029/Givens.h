#ifndef GIVENS_H
#define GIVENS_H

class Givens
{
private:
	double c,s,r;

public:
	Givens();//c,s,rȫ��Ϊ0�Ĺ��캯��
	Givens(double x, double y);//����x,y�����c,s,r
	Givens(double cvalue,double svalue,double rvalue);//ֱ����cvalue,svalue,rvalue��ֵ

	void set(double x,double y);//������x,y����c,s,r
	void update(Matrix,int k);//Givens�任
};




#endif