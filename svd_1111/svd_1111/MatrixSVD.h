#ifndef MATRIXSVD_H
#define MATRIXSVD_H

//��ͷ�ļ�������SVD�п��ܻ��õ��ľ���

#include "Vector.h"
#include <iostream>

class Matrix
{
private:

	//����һά����洢ָ��
	// 0-based and 1-based data pointer
	double *pv0, *pv1;

	// 0-based and 1-based row pointer's pointer
	double **prow0, **prow1;

	int	 nRow;//����
	int	 nColumn;//����
	long nTotal;//��Ԫ�ظ���������row*column

	void init( int rows, int columns );//��ʼ��
	void copyFromArray( const double *v );//��һ�������п���
	void setByScalar( const double &x );//�����������Ԫ�ظ�Ϊͬһ��ֵ
	void destroy();//�ͷ��ڴ�

public:

	// ���캯������������
	Matrix();
	Matrix( const Matrix &A );
	Matrix( int rows, int columns, const double &x = double(0) );
	Matrix( int rows, int columns, const double *v );
	~Matrix();

	// assignments
	Matrix& operator=( const Matrix &A );
	Matrix& operator=( const double &x );

	// accessors
	double* operator[]( int i );
	const double* operator[]( int i ) const;
	double& operator()( int row, int column );
	const double& operator()( int row, int column ) const;

	// double conversion
	operator double*();
	operator const double*() const;
	//        operator double**();
	//        operator const double**() const;

	// others
	long size() const;
	int dim( int dimension ) const;
	int rows() const;
	int cols() const;
	Matrix& resize( int rows, int columns );
	Vector getRow( int row ) const;
	Vector getColumn( int column ) const;
	void setRow( const Vector &v, int row );
	void setColumn( const Vector &v, int column );

	// computed assignment
	Matrix& operator+=( const double& );
	Matrix& operator+=( const Matrix& );
	Matrix& operator-=( const double& );
	Matrix& operator-=( const Matrix& );
	Matrix& operator*=( const double& );

	// WARNING: element-by-element
	Matrix& operator*=( const Matrix& );
	Matrix& operator/=( const double& );

	// WARNING: element-by-element
	Matrix& operator/=( const Matrix& );


};

Vector & mult( const Matrix &A, const Vector &b, Vector &c );
Matrix& mult( const Matrix &A, const Matrix &B,Matrix &C );
Matrix trT( const Matrix& );

void printMatrix(const Matrix &A );

// arithmetic operators

Matrix  negtive( const Matrix & );

Matrix  add( const Matrix &, const double& );

Matrix  add( const double&, const Matrix & );

Matrix  add( const Matrix &, const Matrix & );

Matrix  minus( const Matrix &, const double& );

Matrix  minus( const double&, const Matrix & );

Matrix  minus( const Matrix &, const Matrix & );

Matrix  production( const Matrix &, const double& );

Matrix  production( const double&, const Matrix & );

Matrix  production( const Matrix &, const Matrix & );

Vector  production( const Matrix &, const Vector & );

Matrix  quotient( const Matrix &, const double& );

Matrix  quotient( const double&, const Matrix & );

#endif