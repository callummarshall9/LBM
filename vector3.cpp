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

vector3::vector3(float x, float y, float z) {
	this->x = x;
	this->y = y;
	this->z = z;
}

float vector3::dot(vector3 other) {
	return (this->x * other.get_x() + 
		this->y * other.get_y() + 
		this->z * other.get_z()
	);
}

vector3::vector3() {
	this->zero_vector();
}

float vector3::get_x() {
	return this->x;
}

float vector3::get_y() {
	return this->y;
}

float vector3::get_z() {
	return this->z;
}

void vector3::set_x(float x) {
	this->x = x;
}

void vector3::set_y(float y) {
	this->y = y;
}

void vector3::set_x_y_z(float x, float y, float z) {
	this->set_x(x);
	this->set_y(y);
	this->set_z(z);
}

void vector3::zero_vector() {
	this->x = 0;
	this->y = 0;
	this->z = 0;
}

void vector3::set_z(float z) {
	this->z = z;
}

vector3* vector3::multiply_by_scalar(float scalar) {
	return new vector3(this->x * scalar, this->y * scalar, this->z * scalar);
}

float vector3::norm_square() {
	return (this->x * this->x + this->y * this->y + this->z * this->z);
}
