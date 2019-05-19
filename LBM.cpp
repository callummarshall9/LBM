#include <vector>

class LBM {

public:
	void initialise();
	LBM(int grid_size);
	void output_density();
	void output_velocity();
	void output_test();
	float compute_density_moment(int x, int y, int z);
	vector3* compute_momentum_density_moment(int x, int y, int z);
private:
	float* density_field;
	std::vector<vector3> velocity_field;
	float* equilibrium_distribution;
	int grid_size;
	int scalar_index(int x, int y, int z);
	int scalar_index(int x, int y, int z, int w);
	void output_array(float *array);
	//Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
	const int direction_size = 15;
	vector3 directions[15] = {//c_i vectors.
		vector3(0.0,0.0,0.0),
		vector3(1.0,0.0,0.0),
		vector3(-1.0,0.0,0.0),
		vector3(0.0,1.0,0.0),
		vector3(0.0,-1.0,0.0),
		vector3(0.0,0.0,1.0),
		vector3(0.0,0.0,-1.0),
		vector3(1.0,1.0,1.0),
		vector3(-1.0,-1.0,-1.0),
		vector3(1.0,1.0,-1.0),
		vector3(1.0,-1.0,-1.0),
		vector3(1.0,-1.0,1.0),
		vector3(-1.0,1.0,-1.0),
		vector3(-1.0,1.0,1.0),
		vector3(1.0,-1.0,-1.0) 
	};
	float weights[15] = { (2.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0,72.0),(1.0/72.0) };
	//This will result in a change in the equlibrium function which will be reflected below.
};

LBM::LBM(int grid_size) {
	this->grid_size = grid_size;
	this->density_field = (float*)malloc(sizeof(float) * this->grid_size * this->grid_size * this->grid_size);
	velocity_field.resize(this->grid_size * this->grid_size * this->grid_size);
	this->equilibrium_distribution = (float*)malloc(sizeof(float) * this->direction_size * this->grid_size * this->grid_size * this->grid_size);
	this->initialise();
}

void LBM::initialise() {
	for(int i = 0; i < this->grid_size; i++) {
		for(int j = 0; j < this->grid_size; j++) {
			for(int k = 0; k < this->grid_size; k++) {
				this->density_field[this->scalar_index(i,j,k)] = 1;
				this->velocity_field[this->scalar_index(i,j,k)].zero_vector();
				float density = this->density_field[this->scalar_index(i,j,k)];
				vector3 velocity = this->velocity_field[this->scalar_index(i,j,k)];
				for(int w = 0; w < this->direction_size; w++) {
					float dot_product = velocity.dot(this->directions[w]);
					this->equilibrium_distribution[this->scalar_index(i,j,k,w)] = weights[w] * density * (1.0 + 3.0 * dot_product + 4.5 * dot_product * dot_product - 1.5 * velocity.norm_square());
				}
			}
		}
	}
}

int LBM::scalar_index(int x, int y, int z) {
	return (z * this->grid_size * this->grid_size) + (y * this->grid_size) + x;
}

int LBM::scalar_index(int x, int y, int z, int w) {
	return (w * this->direction_size * this->grid_size * this->grid_size) + 
		(z * this->grid_size * this->grid_size) + (y * this->grid_size) + x;
}

void LBM::output_array(float* array) {
	std::cout << "x,y,z value" << std::endl;
	for(int i = 0; i < this->grid_size; i++) {
		for(int j = 0; j < this->grid_size; j++) {
			for(int k = 0; k < this->grid_size; k++) {
				std::cout << i << "," << j << "," << k << ": " << array[this->scalar_index(i,j,k)] << std::endl;
			}
		}
	}
}

void LBM::output_density() {
	this->output_array(this->density_field);
}

void LBM::output_velocity() {
	
}

float LBM::compute_density_moment(int x, int y, int z) {
	float density = 0;
	for(int i = 0; i < this->grid_size; i = i + 1) {
		density = density + this->equilibrium_distribution[this->scalar_index(x,y,z,i)];
	}
	return density;
}

vector3* LBM::compute_momentum_density_moment(int x, int y, int z) {
	vector3 *momentum_density = new vector3();
	for(int i = 0; i < this->grid_size; i = i + 1) {
		momentum_density->set_x(momentum_density->get_x() + this->directions[i].get_x() * this->equilibrium_distribution[this->scalar_index(x,y,z,i)]);
		momentum_density->set_y(momentum_density->get_y() + this->directions[i].get_y() * this->equilibrium_distribution[this->scalar_index(x,y,z,i)]);
		momentum_density->set_z(momentum_density->get_z() + this->directions[i].get_z() * this->equilibrium_distribution[this->scalar_index(x,y,z,i)]);
	}
	return momentum_density;
}
