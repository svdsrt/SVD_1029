#ifndef SVD_H
#define SVD_H

class svd
{
private:
	int m,n,flag;
	Matrix *U;
	Matrix *V;
	Matrix A,tU,tV,temp;
	double e;
	int l,p,q,i;
	Givens g;
	double x,y;

protected:
	Vector B1,B2;
public:
	svd():A("InputMatrix4.txt"),tU(A.N()),tV(A.N()),temp(A.N()),B1(A.N()),B2(A.N()),g(0,0)
	{
		flag=0;
		m=A.M();
		n=A.N();
		U=new Matrix(m);
		V=new Matrix(n);
		e=0;
		l=0;
		p=0;
		q=0;
		i=0;
		x=0;
		y=0;
	}

	void BiDiag();//���Խǻ�
	void function1();
	void function4();
	void function5();
	void CheckConvergence();//��������

	void BiPrint();
	void UVPrint();

};

class QR
{
public:
	Matrix *P,*Q;
	QR(int n);//���쵥λ��P��Q�Ĺ��캯��
	void Iterative(int i1,int i2,Vector &B1,Vector &B2);//QR�����ļ��㺯��
};


#endif