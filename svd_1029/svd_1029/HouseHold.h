#ifndef HOUSEHOLD_H
#define HOUSEHOLD_H

class HouseHold
{
private:
	Matrix * T;
	double delta;
public:
	HouseHold();//���캯��
	HouseHold(int n);//����n*n�ķ����household
	
	double Delta();//����delta
	Matrix * TMatrix();//���غ���ָ��
	
	void HouseHolder(Vector v);
};
#endif // ! HOUSEHOLD_H
