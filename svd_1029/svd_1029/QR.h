#ifndef QR_H
#define QR_H
//����QR��������Ҫ�Ĳ����뺯��

class QR 
{
private:
	
	Matrix *P,*Q;
public:
	QR(int n);//���쵥λ��P��Q�Ĺ��캯��
	void Iterative(Vector &A,Vector &B,int i1,int i2);//QR�����ļ��㺯��
};




#endif