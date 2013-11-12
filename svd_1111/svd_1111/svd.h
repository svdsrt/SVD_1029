#ifndef SVD_H
#define SVD_H

#include "MatrixSVD.h"

class SVD
{

public:

	SVD();
	~SVD();

	void dec(const Matrix &A );
	Matrix getU() const;
	Matrix getV() const;
	Matrix getSM();
	Vector getSV() const;

	double norm2() const;
	double cond() const;
	int  rank();

private:

	Matrix U;
	Matrix V;
	Vector S;

	void decomposition( Matrix &, Matrix &, Vector &, Matrix & );

};
// class SVD

int min(int x,int y);
int max(int x,int y);
double min(double x,double y);
double max(double x,double y);
void swap(double x,double y);

#endif
// SVD_H