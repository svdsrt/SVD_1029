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

	void BiDiag();//���Խǻ�
	//void CheckConvergence();//��������

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
	QR(int n);//���쵥λ��P��Q�Ĺ��캯��
	void Iterative(int i1,int i2);//QR�����ļ��㺯��
};
*/

#endif