#include "headers/LBM.hpp"
#include "headers/vector3.hpp"

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
