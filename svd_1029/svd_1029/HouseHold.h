#ifndef HOUSEHOLD_H
#define HOUSEHOLD_H

class HouseHold
{
private:
	Matrix * T;
	double delta;
public:
	HouseHold();//构造函数
	HouseHold(int n);//构造n*n的方阵的household
	
	double Delta();//返回delta
	Matrix * TMatrix();//返回函数指针
	
	void HouseHolder(Vector v);
};
#endif // ! HOUSEHOLD_H
