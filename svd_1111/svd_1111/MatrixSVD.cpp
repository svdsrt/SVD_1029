#include "stdafx.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "MatrixSVD.h"
#include <string>
#include <math.h>
#include "Vector.h"

using namespace std;

void Matrix::init( int rows, int columns )
{
	nRow = rows;
	nColumn = columns;
	nTotal = nRow * nColumn;

	pv0 = new double[nTotal];
	prow0 = new double*[nRow];
	prow1 = new double*[nRow];

	double *p = pv0;
	pv1 = pv0 - 1;
	for( int i=0; i<nRow; ++i )
	{
		prow0[i] = p;
		prow1[i] = p-1;
		p += nColumn;
	}

	prow1--;
}


void Matrix::copyFromArray( const double *v )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = v[i];
}
void Matrix::setByScalar( const double &x )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = x;
}

void Matrix::destroy()
{
	if( pv0 == NULL )
		return ;
	else
		delete []pv0;

	if( prow0 != NULL )
		delete []prow0;

	prow1++;
	if( prow1 != NULL )
		delete []prow1;
}


Matrix::Matrix()
	: pv0(0), pv1(0), prow0(0), prow1(0), nRow(0), nColumn(0), nTotal(0)
{
}


Matrix::Matrix( const Matrix &A )
{
	init( A.nRow, A.nColumn );
	copyFromArray( A.pv0 );
}


Matrix::Matrix( int rows, int columns, const double &x )
{
	init( rows,columns );
	setByScalar(x);
}


Matrix::Matrix( int rows, int columns, const double *arrays )
{
	init( rows,columns );
	copyFromArray( arrays );
}


Matrix::~Matrix()
{
	destroy();
}


Matrix& Matrix::operator=( const Matrix &A )
{
	if( pv0 == A.pv0 )
		return *this;

	if( nRow == A.nRow && nColumn == A.nColumn )
		copyFromArray( A.pv0 );
	else
	{
		destroy();
		init( A.nRow, A.nColumn );
		copyFromArray( A.pv0 );
	}

	return *this;
}


/**
* overload evaluate operator = from scalar to matrix
*/

Matrix& Matrix::operator=( const double &x )
{
	setByScalar( x );

	return *this;
}


/**
* overload operator [] for 0-offset access
*/

double* Matrix::operator[]( int i )
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < nRow );
#endif

	return prow0[i];
}


const double* Matrix::operator[]( int i ) const
{
#ifdef BOUNDS_CHECK
	assert( 0 <= i );
	assert( i < nRow );
#endif

	return prow0[i];
}


/**
* overload operator () for 1-offset access
*/

double& Matrix::operator()( int row, int column )
{
#ifdef BOUNDS_CHECK
	assert( 1 <= row );
	assert( row <= nRow ) ;
	assert( 1 <= column);
	assert( column <= nColumn );
#endif

	return  prow1[row][column];
}


const double& Matrix::operator()( int row, int column ) const
{
#ifdef BOUNDS_CHECK
	assert( 1 <= row );
	assert( row <= nRow ) ;
	assert( 1 <= column);
	assert( column <= nColumn );
#endif

	return  prow1[row][column];
}


/**
* double conversion functions
*/

Matrix::operator double*()
{
	return pv0;
}


Matrix::operator const double*() const
{
	return pv0;
}

//
//  Matrix::operator double**()
//{
//	return prow0;
//}
//
//
//  Matrix::operator const double**() const
//{
//	return prow0;
//}


/**
* get the matrix's size
*/

long Matrix::size() const
{
	return nTotal;
}


/**
* get the matrix's dimension
*/

int Matrix::dim( int dimension ) const
{
#ifdef BOUNDS_CHECK
	assert( dimension >= 1);
	assert( dimension <= 2);
#endif

	if( dimension == 1 )
		return nRow;
	else if( dimension == 2 )
		return nColumn;
	else
		return 0;
}


int Matrix::rows() const
{
	return nRow;
}


int Matrix::cols() const
{
	return nColumn;
}


/**
* reallocate matrix's size
*/

Matrix& Matrix::resize( int rows, int columns )
{
	if(  rows == nRow && columns == nColumn )
		return *this;

	destroy();
	init( rows, columns );

	return *this;
}


/**
* get the matrix's row vector
*/

Vector Matrix::getRow( int row ) const
{
#ifdef BOUNDS_CHECK
	assert( row >= 0 );
	assert( row < nRow );
#endif

	Vector tmp( nColumn );
	for( int j=0; j<nColumn; ++j )
		tmp[j] = prow0[row][j];

	return tmp;
}


/**
* get the matrix's column vector
*/

Vector Matrix::getColumn( int column ) const
{

	Vector tmp( nRow );
	for( int i=0; i<nRow; ++i )
		tmp[i] = prow0[i][column];

	return tmp;
}


/**
* set the matrix's row vector
*/

void Matrix::setRow( const Vector &v, int row )
{

	for( int j=0; j<nColumn; ++j )
		prow0[row][j] = v[j];
}


/**
* set the matrix's column vector
*/

void Matrix::setColumn( const Vector &v, int column )
{
#ifdef BOUNDS_CHECK
	assert( column >= 0 );
	assert( column < nColumn );
	assert( v.dim() == nRow );
#endif

	for( int i=0; i<nRow; ++i )
		prow0[i][column] = v[i];
}


/**
* compound assignment operators +=
*/

Matrix& Matrix::operator+=( const double &x )
{
	double **rowPtr = prow0;
	double *colPtr = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtr = *rowPtr++;
		for( int j=0; j<nColumn; ++j )
			*colPtr++ += x;
	}

	return *this;
}


Matrix& Matrix::operator+=( const Matrix &rhs )
{

	double **rowPtrL = prow0;
	double *colPtrL = 0;
	double **rowPtrR = rhs.prow0;
	const double *colPtrR = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtrL = *rowPtrL++;
		colPtrR = *rowPtrR++;
		for( int j=0; j<nColumn; ++j )
			*colPtrL++ += *colPtrR++;
	}

	return *this;
}


/**
* compound assignment operators -=
*/

Matrix& Matrix::operator-=( const double &x )
{
	double **rowPtr = prow0;
	double *colPtr = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtr = *rowPtr++;
		for( int j=0; j<nColumn; ++j )
			*colPtr++ -= x;
	}

	return *this;
}


Matrix& Matrix::operator-=( const Matrix &rhs )
{

	double **rowPtrL = prow0;
	double *colPtrL = 0;
	double **rowPtrR = rhs.prow0;
	const double *colPtrR = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtrL = *rowPtrL++;
		colPtrR = *rowPtrR++;
		for( int j=0; j<nColumn; ++j )
			*colPtrL++ -= *colPtrR++;
	}

	return *this;
}


/**
* compound assignment operators *=
*/

Matrix& Matrix::operator*=( const double &x )
{
	double **rowPtr = prow0;
	double *colPtr = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtr = *rowPtr++;
		for( int j=0; j<nColumn; ++j )
			*colPtr++ *= x;
	}

	return *this;
}

// WARNING: this is element-by-element multiplication

Matrix& Matrix::operator*=( const Matrix &rhs )
{

	double **rowPtrL = prow0;
	double *colPtrL = 0;
	double **rowPtrR = rhs.prow0;
	const double *colPtrR = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtrL = *rowPtrL++;
		colPtrR = *rowPtrR++;
		for( int j=0; j<nColumn; ++j )
			*colPtrL++ *= *colPtrR++;
	}

	return *this;
}


/**
* compound assignment operators /=
*/

Matrix& Matrix::operator/=( const double &x )
{
	double **rowPtr = prow0;
	double *colPtr = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtr = *rowPtr++;
		for( int j=0; j<nColumn; ++j )
			*colPtr++ /= x;
	}

	return *this;
}

// WARNING: this is element-by-element division

Matrix& Matrix::operator/=( const Matrix &rhs )
{
	double **rowPtrL = prow0;
	double *colPtrL = 0;
	double **rowPtrR = rhs.prow0;
	const double *colPtrR = 0;

	for( int i=0; i<nRow; ++i )
	{
		colPtrL = *rowPtrL++;
		colPtrR = *rowPtrR++;
		for( int j=0; j<nColumn; ++j )
			*colPtrL++ /= *colPtrR++;
	}

	return *this;
}


/**
* Overload the output stream function.
*/

ostream& operator<<( ostream &out, const Matrix &A )
{
	int rows = A.rows();
	int columns = A.cols();

	out << "size: " << rows << " by " << columns << "\n";
	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			out << A[i][j] << "\t";
		out << "\n";
	}

	return out;
}


/**
* Overload the intput stream function.
*/

istream& operator>>( istream &in, Matrix &A )
{
	int rows, columns;
	in >> rows >> columns;

	if( !( rows == A.rows() && columns == A.cols() ) )
		A.resize( rows, columns );

	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			in >> A[i][j];

	return in;
}


/**
* get negative matrix
*/


Matrix negtive( const Matrix &A )
{
	int rows = A.rows();
	int columns = A.cols();

	Matrix tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = -A[i][j];

	return tmp;
}


/**
* matrix-scalar addition
*/

Matrix add( const Matrix &A, const double &x )
{
	Matrix tmp( A );
	return tmp += x;
}


/**
* matrix-matrix addition
*/

Matrix add( const Matrix &A1, const Matrix &A2 )
{
	Matrix tmp( A1 );
	return tmp += A2;
}

/**
* matrix-matrix subtraction
*/

Matrix minus( const Matrix &A1, const Matrix &A2 )
{
	Matrix tmp( A1 );
	return tmp -= A2;
}


Matrix production( const Matrix &A, const double &x )
{
	Matrix tmp( A );
	return tmp *= x;
}

/**
* matrix-matrix multiplication
*/

Matrix production( const Matrix &A1, const Matrix &A2 )
{


	int rows = A1.rows();
	int columns = A2.cols();
	//	int K = A1.cols();

	Matrix tmp( rows, columns );
	//	for( int i=0; i<rows; ++i )
	//		for( int j=0; j<columns; ++j )
	//		{
	//            tmp[i][j] = 0;
	//			for( int k=0; k<K; ++k )
	//			    tmp[i][j] += A1[i][k] * A2[k][j];
	//		}

	mult( A1, A2, tmp );

	return tmp;
}


/**
* matrix-vector multiplication
*/

Vector production( const Matrix &A, const Vector &b )
{


	int rows = A.rows();
	//	int columns = A.cols();

	Vector tmp(rows);
	//	for( int i=0; i<rows; ++i )
	//	{
	//		double sum = 0;
	//		for( int j=0; j<columns; ++j )
	//			sum += A[i][j] * v[j];
	//		tmp[i] = sum;
	//	}

	mult( A, b, tmp );

	return tmp;
}


/**
* matrix-scalar division
*/

Matrix quotient( const Matrix &A, const double &x )
{
	Matrix tmp( A );
	return tmp /= x;
}


Matrix quotient( const double &x, const Matrix &A )
{
	int rows = A.rows();
	int clumns = A.cols();

	Matrix tmp( rows,clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = x / A[i][j];

	return tmp;
}

Matrix& mult( const Matrix &A, const Matrix &B,
			 Matrix &C )
{
	int M = A.rows();
	int N = B.cols();
	int K = A.cols();


	C.resize( M, N );
	double        sum;
	const double  *pRow,
		*pCol;

	for( int i=0; i<M; i++ )
		for( int j=0; j<N; ++j )
		{
			pRow  = &A[i][0];
			pCol  = &B[0][j];
			sum = 0;

			for( int k=0; k<K; ++k )
			{
				sum += (*pRow) * (*pCol);
				pRow++;
				pCol += N;
			}
			C[i][j] = sum;
		}
		return C;
}


Vector & mult( const Matrix &A, const Vector &b,
			  Vector &c )
{
	int M = A.rows();
	int N = A.cols();


	c.resize( M );
	double   sum;
	const double  *pRow,
		*pCol;

	for( int i=0; i<M; i++ )
	{
		pRow  = &A[i][0];
		pCol  = &b[0];
		sum = 0;

		for( int j=0; j<N; ++j )
		{
			sum += (*pRow) * (*pCol);
			pRow++;
			pCol++;
		}
		c[i] = sum;
	}
	return c;
}

/**
* matrix tranpose
*/

Matrix trT( const Matrix &A )
{
	int rows = A.cols();
	int clumns = A.rows();

	Matrix tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = A[j][i];

	return tmp;
}


void printMatrix(const Matrix &A )
{
	int rows = A.rows();
	int columns = A.cols();

	cout << "size: " << rows << " by " << columns << "\n";
	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			cout << A[i][j] << "\t";
		cout << "\n";
	}

}

