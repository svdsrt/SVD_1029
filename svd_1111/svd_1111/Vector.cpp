#include "stdafx.h"
#include <iostream>
#include <cmath>
#include "Vector.h"

/**
* initialize
*/

void Vector::init( int length )
{
	pv0 = new double[length];

	pv1 = pv0 - 1;
	nRow = length;
}


/**
* copy vector from normal array
*/

inline void Vector::copyFromArray( const double *v )
{
	for( int i=0; i<nRow; ++i )
		pv0[i] = v[i];
}


/**
* set vector by a scalar
*/

inline void Vector::setByScalar( const double &x )
{
	for( int i=0; i<nRow; ++i )
		pv0[i] = x;
}


/**
* destroy the vector
*/

void Vector::destroy()
{
	if( pv0 == NULL )
		return;

	delete []pv0;

	pv0 = NULL;
	pv1 = NULL;
}


/**
* constructors and destructor
*/

Vector::Vector()
	: pv0(0), pv1(0), nRow(0)
{
}


Vector::Vector( const Vector &v )
	: pv0(0), pv1(0), nRow(0)
{
	init( v.nRow );
	copyFromArray( v.pv0 );
}


Vector::Vector( int length, const double &x )
	:  pv0(0), pv1(0), nRow(0)
{
	init( length );
	setByScalar( x );
}


Vector::Vector( int length, const double *array )
	:  pv0(0), pv1(0), nRow(0)
{
	init( length );
	copyFromArray( array );
}


Vector::~Vector()
{
	destroy();
}


/**
* overload evaluate operator= from vector to vector
*/

Vector& Vector::operator=( const Vector &v )
{
	if( pv0 == v.pv0 )
		return *this;

	if( nRow == v.nRow )
		copyFromArray( v.pv0 );
	else
	{
		destroy();
		init( v.nRow );
		copyFromArray( v.pv0 );
	}

	return *this;
}


/**
* overload evaluate operator= from scalar to vector
*/

inline Vector& Vector::operator=( const double &x )
{
	setByScalar( x );

	return *this;
}


/**
* overload operator [] for 0-offset access
*/

double& Vector::operator[]( int i )
{
	return pv0[i];
}


const double& Vector::operator[]( int i ) const
{
	return pv0[i];
}


/**
* overload operator () for 1-offset access
*/

double& Vector::operator()( int i )
{
	return pv1[i];
}


const double& Vector::operator()( int i ) const
{
	return pv1[i];
}


/**
* iterators
*/

inline Vector::iterator Vector::begin()
{
	return pv0;
}


inline Vector::const_iterator Vector::begin() const
{
	return pv0;
}


inline Vector::iterator Vector::end()
{
	return pv0 + nRow;
}


inline Vector::const_iterator Vector::end() const
{
	return pv0 + nRow;
}


/**
* double conversion functions
*/

inline Vector::operator double*()
{
	return pv0;
}


inline Vector::operator const double*() const
{
	return pv0;
}


/**
* get the vector's total size
*/

int Vector::size() const
{
	return  nRow;
}


/**
* get the vector's dimension
*/

inline int Vector::dim() const
{
	return  nRow;
}


/**
* reallocate vector's size
*/

Vector& Vector::resize( int length )
{
	if( nRow == length )
		return *this;

	destroy();
	init( length );

	return *this;
}


/**
* compound assignment operators +=
*/

Vector& Vector::operator+=( const double &x )
{
	iterator itr = (*this).begin();
	while( itr != (*this).end() )
		*itr++ += x;

	return *this;
}


Vector& Vector::operator+=( const Vector &rhs )
{

	iterator itrL = (*this).begin();
	const_iterator itrR = rhs.begin();
	while( itrL != (*this).end() )
		*itrL++ += *itrR++;

	return *this;
}


/**
* compound assignment operators -=
*/

Vector& Vector::operator-=( const double &x )
{
	iterator itr = (*this).begin();
	while( itr != (*this).end() )
		*itr++ -= x;

	return *this;
}


Vector& Vector::operator-=( const Vector &rhs )
{

	iterator itrL = (*this).begin();
	const_iterator itrR = rhs.begin();
	while( itrL != (*this).end() )
		*itrL++ -= *itrR++;

	return *this;
}


/**
* compound assignment operators *=
*/

Vector& Vector::operator*=( const double &x )
{
	iterator itr = (*this).begin();
	while( itr != (*this).end() )
		*itr++ *= x;

	return *this;
}


Vector& Vector::operator*=( const Vector &rhs )
{

	iterator itrL = (*this).begin();
	const_iterator itrR = rhs.begin();
	while( itrL != (*this).end() )
		*itrL++ *= *itrR++;

	return *this;
}


/**
* compound assignment operators /=
*/

Vector& Vector::operator/=( const double &x )
{
	iterator itr = (*this).begin();
	while( itr != (*this).end() )
		*itr++ /= x;

	return *this;
}


Vector& Vector::operator/=( const Vector &rhs )
{

	iterator itrL = (*this).begin();
	const_iterator itrR = rhs.begin();
	while( itrL != (*this).end() )
		*itrL++ /= *itrR++;

	return *this;
}