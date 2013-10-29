#ifndef MATRIXSVD_H
#define MATRIXSVD_H

//��ͷ�ļ�������SVD�п��ܻ��õ��ľ�����������Լ�����


class Matrix
{
private:
	double *matrix;//һά����ָ�룬���Դ������Ԫ��
	int m,n,mn;//m��������n������

public:
	Matrix(int Nnum);//���캯����n*n�ĵ�λ����
	Matrix(int Mnum,int Nnum);//���캯���� Mnum*Nnum�ľ���
	Matrix(std::string s);//���ļ��ж�ȡ�������

	double a(int i,int j);//����a_{ij}��ֵ
	double &operator[](int index);//����matrix[index]
	int M();//���ؾ�������
	int N();//���ؾ�������

	void Add(Matrix A);//��
	void Minus(Matrix A);//��
	void NumProd(double a);//����
	void DotProd(Matrix A,Matrix B);//����A����B��������
	void Trans(Matrix &trans);//ת��

	void set(int i,int j,double value);//�޸�ֵ

	double Max();//������������

	void print();//����Ļ�ϴ�ӡ����
	void fprint();//���ļ��д�ӡ����

	~Matrix();//��������
};

class Vector
{
private:
	double *vector;
	int n;
public:
	Vector(int num); //����nά����

	double &operator[](int index);//���������е�indexά��ֵ����0��ʼ
	int N();//��������ά��

	void Add(Vector A);//��
	void Minus(Vector A);//��
	void NumProd(double a);//����
	double DotProd(Vector A);//���

	void set(int i,double value);//�޸�ֵ

	double modular();//����ģ
	double Max();//������������

	void Row(Matrix &A,int k);//ȡ�����k��������
	void Col(Matrix &A,int k);//ȡ�����k��������

	void HRow(Matrix &A,int k);//Household��ȡ�������ķ���
	void HCol(Matrix &A,int k);//Household��ȡ�������ķ���

	void Normalize();//��һ��
	void Span(Matrix T);//���ɾ���
	void print();//����Ļ�ϴ�ӡ����
	void fprint();//���ļ��д�ӡ����

	~Vector();//��������
};



#endif