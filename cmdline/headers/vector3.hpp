#ifndef VECTOR_H_FILE
#define VECTOR_H_FILE

class vector3 {
private:
	float x;
	float y;
	float z;
public:
	vector3(float x, float y, float z);
	vector3();
	float get_x();
	float get_y();
	float get_z();
	void set_x(float x);
	void set_y(float y);
	void set_z(float z);
	void set_x_y_z(float x, float y, float z);
	void zero_vector();
	float dot(vector3 other);
	float norm_square();
	vector3* multiply_by_scalar(float scalar);
};

#endif
