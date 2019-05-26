#ifndef VECTOR_H_FILE
#define VECTOR_H_FILE

template<class myType>
class vector3 {

public:
	vector3(myType x, myType y, myType z);
	vector3();
	myType dot(vector3<myType> other);
	myType norm_square();
	vector3* multiply_by_scalar(myType scalar);
	myType x;
	myType y;
	myType z;
};

#endif
