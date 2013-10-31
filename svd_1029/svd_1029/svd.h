#ifndef SVD_H
#define SVD_H

class svd
{
private:
	int m,n;
	Matrix *U;
	Matrix *V;
	Matrix A,tU,tV,temp;
protected:
	Vector B1,B2;
public:
	svd():A("InputMatrix1.txt"),tU(A.N()),tV(A.N()),temp(A.N()),B1(A.N()),B2(A.N())
	{
		m=A.M();
		n=A.N();
		U=new Matrix(m);
		V=new Matrix(n);
	}

	void BiDiag();//二对角化
	//void CheckConvergence();//收敛检验

	void BiPrint();
	void UVPrint();
};
/*
class QR 
{
private:
	Matrix *P,*Q;
	friend class svd;
public:
	QR(int n);//构造单位阵P、Q的构造函数
	void Iterative(int i1,int i2);//QR迭代的计算函数
};
*/

#endif