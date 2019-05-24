#include "LBM/headers/LBM.hpp"
#include "LBM/headers/vector3.hpp"
#include <iomanip>

LBM::LBM(int grid_size) {
	this->grid_size = grid_size;
	int box_flatten_length = grid_size * grid_size * grid_size;
	int equilibrium_flatten_length = box_flatten_length * this->direction_size;
	this->density_field = new float[box_flatten_length];
	velocity_field.resize(box_flatten_length);
	this->equilibrium_distribution = new float[equilibrium_flatten_length];
	this->initialise();
}

void LBM::free_memory() {
	free(density_field);
	free(equilibrium_distribution);
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
					this->equilibrium_distribution[this->scalar_index(i,j,k,w)] =
						weights[w] * density * (1.0 + 3.0 * dot_product + 4.5 * dot_product * dot_product - 1.5 * velocity.norm_square());
				}
			}
		}
	}
}

int LBM::scalar_index(int x, int y, int z) {
	return (z * this->grid_size * this->grid_size) + (y * this->grid_size) + x;
}

int LBM::scalar_index(int x, int y, int z, int w) {
	int index = x + y * this->grid_size + z * this->grid_size * grid_size + w * grid_size * grid_size * grid_size;
	return index;
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
	for(int i = 0; i < this->direction_size; i = i + 1) {
		density = density + this->equilibrium_distribution[this->scalar_index(x,y,z,i)];
		this->density_field[this->scalar_index(x,y,z)] = density;
	}
	return density;
}

void LBM::compute_density_moment() {
	for(int x = 0; x < this->grid_size; x = x + 1) {
		for(int y = 0; y < this->grid_size; y = y + 1) {
			for(int z = 0; z < this->grid_size; z = z + 1) {
				this->compute_density_moment(x,y,z);
			}
		}
	}
}

void LBM::compute_momentum_density_moment(int x, int y, int z) {
	vector3 momentum_density;
	for(int i = 0; i < this->direction_size; i = i + 1) {
		momentum_density.set_x(momentum_density.get_x() + this->directions[i].get_x() * this->equilibrium_distribution[this->scalar_index(x,y,z,i)]);
		momentum_density.set_y(momentum_density.get_y() + this->directions[i].get_y() * this->equilibrium_distribution[this->scalar_index(x,y,z,i)]);
		momentum_density.set_z(momentum_density.get_z() + this->directions[i].get_z() * this->equilibrium_distribution[this->scalar_index(x,y,z,i)]);
	}
	this->velocity_field[this->scalar_index(x,y,z)] = vector3(
		momentum_density.get_x() / this->density_field[this->scalar_index(x,y,z)],
		momentum_density.get_y() / this->density_field[this->scalar_index(x,y,z)],
		momentum_density.get_z() / this->density_field[this->scalar_index(x,y,z)]
	);
}

void LBM::compute_momentum_density_moment() {
	for(int x = 0; x < this->grid_size; x = x + 1) {
		for(int y = 0; y < this->grid_size; y = y + 1) {
			for(int z = 0; z < this->grid_size; z = z + 1) {
				this->compute_momentum_density_moment(x,y,z);
			}
		}
	}
}

void LBM::stream() {
	float* new_distribution = (float*)malloc(sizeof(float) * this->grid_size * this->grid_size * this->grid_size * this->direction_size);
	for(int x = 0; x < this->grid_size; x = x + 1) {
		for(int y = 0; y < this->grid_size; y = y + 1) {
			for(int z = 0; z < this->grid_size; z = z + 1) {
				for(int i = 0; i < this->direction_size; i = i + 1) {
					int xmd = (this->grid_size + x - (int)this->directions[i].get_x()) % this->grid_size;
					int ymd = (this->grid_size + y - (int)this->directions[i].get_y()) % this->grid_size;
					int zmd = (this->grid_size + z - (int)this->directions[i].get_z()) % this->grid_size;
					new_distribution[this->scalar_index(x,y,z,i)] = this->equilibrium_distribution[this->scalar_index(xmd,ymd,zmd,i)];
				}
			}
		}
	}
	free(this->equilibrium_distribution);
	this->equilibrium_distribution = new_distribution;
}

void LBM::collision() {//Performs the collision step.
	for(int x = 0; x < this->grid_size; x = x + 1) {
		for(int y = 0; y < this->grid_size; y = y + 1) {
			for(int z = 0; z < this->grid_size; z = z + 1) {
				vector3 velocity = this->velocity_field[this->scalar_index(x,y,z)];
				float density = this->density_field[this->scalar_index(x,y,z)];
				for(int i = 0; i < this->direction_size; i = i + 1) {
					float dot_product = velocity.dot(this->directions[i]);
					this->equilibrium_distribution[this->scalar_index(x,y,z,i)] =
						weights[i] * density * (1.0 + 3.0 * dot_product + 4.5 * dot_product * dot_product - 1.5 * velocity.norm_square());
				}
			}
		}
	}
}

void LBM::perform_timestep() {
	this->stream();
	this->compute_density_moment();
	this->compute_momentum_density_moment();
	this->collision();
}

void LBM::output_lbm_data(std::string filename) {
	std::ofstream output_stream;
  output_stream.open (filename, std::ofstream::out | std::ofstream::app);
	if(!output_stream.is_open()) {
		std::cout << "Cannot open filename: " + filename + " for opening." << std::endl;
		return;
	}
	output_stream << "x,y,z,p,u_x,u_y,u_z" << std::endl;
	output_stream << std::fixed;
	output_stream << std::setprecision(6);
	for(int x = 0; x < grid_size; x++) {
		for(int y = 0; y < grid_size; y++) {
			for(int z = 0; z < grid_size; z++) {
				float density = density_field[this->scalar_index(x,y,z)];
				vector3 fluid_velocity = velocity_field[this->scalar_index(x,y,z)];
				output_stream << x << "," << y << "," << z << "," << density << "," <<
					(float)fluid_velocity.get_x() << "," << (float)fluid_velocity.get_y() << "," << (float)fluid_velocity.get_z() << std::endl;
			}
		}
	}
	output_stream.close();
}
