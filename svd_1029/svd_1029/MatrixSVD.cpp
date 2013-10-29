#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "MatrixSVD.h"
#include <string>
#include <math.h>

using namespace std;

//Vector类***************************************************************

Vector::Vector(int num)
{
	vector=new double[num];
	n=num;
	int i=0;
	for(i=0;i<n;i++)
		vector[i]=0;
}

double & Vector::operator[](int index)
{
	return vector[index];
}

int Vector::N()
{
	return n;
}

void Vector::Add(Vector A)
{
	if(n!=A.N())
	{
		cout<<"ERRO:向量维数不同，不能相加"<<endl;
		cout<<"维数分别为"<<n<<"和"<<A.N()<<endl;
		cout<<"发生错误在";
		print();
		cout<<" + ";
		A.print();
		cout<<endl;
		system("pause");
		exit(1);
	}
	int i=0;
	for(i=0;i<n;i++)
		vector[i]=vector[i]+A[i];
}

void Vector::Minus(Vector A)
{
	if(n!=A.N())
	{
		cout<<"ERRO:向量维数不同，不能相减"<<endl;
		cout<<"维数分别为"<<n<<"和"<<A.N()<<endl;
		cout<<"发生错误在";
		print();
		cout<<" - ";
		A.print();
		cout<<endl;
		system("pause");
		exit(1);
	}
	int i=0;
	for(i=0;i<n;i++)
		vector[i]=vector[i]-A[i];
}

void Vector::NumProd(double a)
{
	int i=0;
	for(i=0;i<n;i++)
		vector[i]=vector[i]*a;
}

double Vector::DotProd(Vector A)
{
	if(n!=A.N())
	{
		cout<<"ERRO:向量维数不同，不能相乘"<<endl;
		cout<<"维数分别为"<<n<<"和"<<A.N()<<endl;
		cout<<"发生错误在";
		print();
		cout<<" ` ";
		A.print();
		cout<<endl;
		system("pause");
		exit(1);
	}
	int i=0;
	double sum=0;
	for(i=0;i<n;i++)
		sum+=vector[i]*A[i];
	return sum;
}


void Vector::set(int i,double value)
{
	vector[i]=value;
}

double Vector::modular()
{
	double sum=0;
	for(int i=0;i<n;i++)
		sum+=vector[i]*vector[i];
	return sqrt(sum);
}

double Vector::Max()
{
	int i=1;
	double temp=vector[0];
	for(;i<n;i++)
		if(temp<vector[i])
			temp=vector[i];
	return temp;
}



void Vector::print()
{
	int i=0;
	cout<<"(";
	for(i=0;i<n;i++)
		cout<<vector[i]<<" ";
	cout<<")"<<endl;;
}

void Vector::fprint()
{
	ofstream outfile;
	outfile.open("outfile.txt",ios::app,0);
	int i=0;
	cout<<"(";
	for(i=0;i<n;i++)
		outfile<<vector[i]<<" ";
	cout<<")";
	outfile.close();
}

Vector::~Vector()
{
	delete vector;
}
//Matrix类***************************************************************

Matrix::Matrix(int Nnum)
{
	matrix=new double[Nnum*Nnum];
	n=Nnum;
	m=Nnum;
	mn=m*n;
	int i=0;
	for(;i<mn;i++)
		matrix[i]=0;
	for(i=0;i<m;i++)
		matrix[i*n+i]=1;
}

Matrix::Matrix(int Mnum,int Nnum)
{
	matrix=new double[Mnum*Nnum];
	m=Mnum;
	n=Nnum;
	mn=m*n;
	int i=0;
	for(;i<mn;i++)
		matrix[i]=0;
}

Matrix::Matrix(std::string s)
{
	ifstream infile;
	infile.open(s);
	infile>>m;
	infile>>n;

	mn=m*n;
	matrix=new double[mn];
	int i=0,j=0;
	for(;i<m;i++)
		for(j=0;j<n;j++)
			infile>>matrix[i*n+j];
	cout<<"*********数据读取完毕*********"<<endl;
}

double Matrix::a(int i,int j)
{
	return matrix[i*n+j];
}

double & Matrix::operator[](int index)
{
	return matrix[index];
}

int Matrix::M()
{
	return m;
}

int Matrix::N()
{
	return n;
}

void Matrix::Add(Matrix A)
{
	if(n!=A.N())
	{
		cout<<"ERRO:列数不匹配"<<endl;
		cout<<"列数分别为：";
		cout<<n<<"和"<<A.N()<<endl;
		system("pause");
		exit(1);
	}
	else if(m!=A.M())
	{
		cout<<"ERRO:行数不匹配"<<endl;
		cout<<"行数分别为：";
		cout<<n<<"和"<<A.N()<<endl;
		system("pause");
		exit(1);
	}

	int i=0;
	for(;i<mn;i++)
		matrix[i]=matrix[i]+A[i];
}
void Matrix::Minus(Matrix A)
{
	if(n!=A.N())
	{
		cout<<"ERRO:列数不匹配"<<endl;
		cout<<"列数分别为：";
		cout<<n<<"和"<<A.N()<<endl;
		system("pause");
		exit(1);
	}
	else if(m!=A.M())
	{
		cout<<"ERRO:行数不匹配"<<endl;
		cout<<"行数分别为：";
		cout<<n<<"和"<<A.N()<<endl;
		system("pause");
		exit(1);
	}

	int i=0;
	for(;i<mn;i++)
		matrix[i]=matrix[i]-A[i];
}
void Matrix::NumProd(double a)
{
	int i=0;
	for(;i<mn;i++)
		matrix[i]=matrix[i]*a;
}

void Matrix::DotProd(Matrix A,Matrix B)
{
	if(A.N()!=B.M())
	{
		cout<<"ERRO:列数不匹配"<<endl;
		cout<<"列数分别为：";
		cout<<n<<"和"<<A.N()<<endl;
		system("pause");
		exit(1);
	}

	int i=0,j=0,k=0;
	for (;i<m;i++)
		for (j=0;j<A.M();j++)
		{
			matrix[i*B.N()+j]=0;
			for (k=0;k<A.N();k++)
				matrix[i*n+j]+=A[i*A.N()+k]*B[k*B.N()+j];
		}
}

void Matrix::Trans(Matrix &trans)
{
	if(m!=trans.N())
	{
		cout<<"ERRO:列数不匹配"<<endl;
		cout<<"列数分别为：";
		cout<<m<<"和"<<trans.N()<<endl;
		system("pause");
		exit(1);
	}
	else if(n!=trans.M())
	{
		cout<<"ERRO:行数不匹配"<<endl;
		cout<<"行数分别为：";
		cout<<n<<"和"<<trans.M()<<endl;
		system("pause");
		exit(1);
	}
	int i=0,j=0;
	for(;i<n;i++)
		for(j=0;j<m;j++)
			trans[i*m+j]=matrix[j*n+i];
}

void Matrix::set(int i,int j,double value)
{
	matrix[i*n+j]=value;
}

void Vector::Row(Matrix &A,int k)
{
	int i=0;
	int num=A.N();
	for(i=0;i<n;i++)
		vector[i]=A[k*num+i];
}

void Vector::Col(Matrix &A,int k)
{
	int i=0;
	int num=A.N();
	for(i=0;i<n;i++)
		vector[i]=A[i*num+k];
}

void Vector::HRow(Matrix &A,int k)//此处默认m>n，可能需要修改
{
	static Vector v(n-k-1);
	int i=0;
	int num=A.N();
	for(;i<v.N();i++)
		vector[i]=A[k*num+k+i+1];
}

void Vector::HCol(Matrix &A,int k)//此处默认m>n，可能需要修改
{
	static Vector v(n-k);
	int i=0;
	int num=A.N();
	for(;i<v.N();i++)
		vector[i]=A[k*num+i*num+k];
}

void Vector::Normalize()
{
	int i=0,sum=0;
	for(i=0;i<n;i++)
		sum+=vector[i];
	for(i=0;i<n;i++)
		vector[i]=vector[i]/sum;
}
void Vector::Span(Matrix T)
{
	if(T.M()!=T.N())
	{
		cout<<"ERRO：span 发生错误"<<endl;
		cout<<"传入矩阵行列数不等"<<endl;
	}
	else if(T.N()!=n)
	{
		cout<<"ERRO：span 发生错误"<<endl;
		cout<<"传入矩阵行列数维数与向量不等"<<endl;
	}
	int i,j;
	for(int i=0;i<n;i++)
		for(int j=0;j<n;j++)
			T[i*n+j]=vector[i]*vector[j];
}

double Matrix::Max()
{
	int i=1;
	double temp=matrix[0];
	for(;i<mn;i++)
		if(temp<matrix[i])
			temp=matrix[i];
	return temp;
}

void Matrix::print()
{

	int i=0,j=0;
	cout<<endl;
	for(;i<m;i++)
	{
		cout<<"| ";
		for(j=0;j<n;j++)
			cout<<matrix[i*n+j]<<" ";
		cout<<"|"<<endl;
	}
	cout<<endl;
}

void Matrix::fprint()
{
	ofstream outfile;
	outfile.open("outfile.txt",ios::app,0);
	int i=0,j=0;
	for(;i<m;i++)
	{
		outfile<<"| ";
		for(j=0;j<n;j++)
			outfile<<matrix[i*n+j]<<" ";
		outfile<<"|"<<endl;
	}
	outfile.close();
}

Matrix::~Matrix()
{
	delete []matrix;
}