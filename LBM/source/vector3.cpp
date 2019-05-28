#include "vector3.hpp"

template<class myType>
vector3<myType>::vector3(myType x, myType y, myType z) : x(x), y(y), z(z) {
}

template<class myType>
myType vector3<myType>::dot(vector3<myType> other) {
	return (x * other.get_x() +
		y * other.get_y() +
		z * other.get_z()
	);
}

template<class myType>
vector3<myType>::vector3() : x(0), y(0), z(0) {
}

template<class myType>
myType vector3<myType>::norm_square() {
	return (x * x + y * y + z * z);
}
