#ifndef QR_H
#define QR_H
//包含QR迭代所需要的参量与函数

class QR 
{
private:
	
	Matrix *P,*Q;
public:
	QR(int n);//构造单位阵P、Q的构造函数
	void Iterative(Vector &A,Vector &B,int i1,int i2);//QR迭代的计算函数
};




#endif