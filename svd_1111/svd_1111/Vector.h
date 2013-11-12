#ifndef VECTOR_H
#define VECTOR_H

class Vector
{
private:

	// data pointer for 0-offset indexing
	double *pv0;

	// data pointer for 1-offset indexing
	double *pv1;

	// the row number of vector
	int	 nRow;

	void init( int length );
	void copyFromArray( const double *v );
	void setByScalar( const double &x );
	void destroy();


public:

	typedef         double*   iterator;
	typedef const   double*   const_iterator;

	// constructors and destructor
	Vector();
	Vector( const Vector &v );
	Vector( int length, const double &x = double(0) );
	Vector( int length, const double *array );
	~Vector();

	// assignments
	Vector& operator=( const Vector &v );
	Vector& operator=( const double &x );

	// accessors
	double& operator[]( int i );
	const double& operator[]( int i ) const;
	double& operator()( int i );
	const double& operator()( int i ) const;

	// iterators
	iterator begin();
	const_iterator begin() const;
	iterator end();
	const_iterator end() const;

	// double conversion
	operator double*();
	operator const double*() const;

	// others
	int size() const;
	int dim() const;
	Vector& resize( int length );

	// computed assignment
	Vector& operator+=( const double& );
	Vector& operator-=( const double& );
	Vector& operator*=( const double& );
	Vector& operator/=( const double& );
	Vector& operator+=( const Vector& );
	Vector& operator-=( const Vector& );
	Vector& operator*=( const Vector& );
	Vector& operator/=( const Vector& );


};

#endif