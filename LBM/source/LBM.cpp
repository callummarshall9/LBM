#include "LBM.hpp"
#include "vector3.cpp"
#include <iomanip>

LBM::LBM(int grid_size) : NX(grid_size),NY(grid_size),NZ(grid_size) {
	initialise();
}

LBM::LBM(int nx, int ny, int nz) : NX(nx), NY(ny), NZ(nz) {
	 initialise();
}

void LBM::initialise() {
	int box_flatten_length = NX * NY * NZ;
	int equilibrium_flatten_length = box_flatten_length * direction_size;
	density_field = new double[box_flatten_length];
	velocity_field = new vector3<double>[box_flatten_length];
	equilibrium_distribution = new double[equilibrium_flatten_length];
	for(int i = 0; i < NX * NY * NZ; i++) {
		density_field[i] = 1;
	}
	for(int i = 0; i < NX; i++) {
		for(int j = 0; j < NY; j++) {
			for(int k = 0; k < NZ; k++) {
				for(int w = 0; w < direction_size; w++) {
					double dot_product = (double)velocity_field[scalar_index(i,j,k)].x * (double)directions[w].x + (double)velocity_field[scalar_index(i,j,k)].y * (double)directions[w].y +
						(double)velocity_field[scalar_index(i,j,k)].z * (double)directions[w].z;
					double feq = weights[w] * density_field[scalar_index(i,j,k)] * (1.0 + 3.0 * dot_product + 4.5 * dot_product * dot_product - 1.5 * velocity_field[scalar_index(i,j,k)].norm_square());
					equilibrium_distribution[scalar_index(i,j,k,w)] = feq;
				}
			}
		}
	}
}

void LBM::free_memory() {
	delete[] density_field;
	delete[] equilibrium_distribution;
}

int LBM::scalar_index(int x, int y, int z) {
	return (z * NX * NY) + (y * NX) + x;
}

int LBM::scalar_index(int x, int y, int z, int w) {
	int index = x + y * NX + z * NX * NY + w * NX * NY * NZ;
	return index;
}

void LBM::output_array(double* array) {
	std::cout << "x,y,z value" << std::endl;
	for(int i = 0; i < NX; i++) {
		for(int j = 0; j < NY; j++) {
			for(int k = 0; k < NZ; k++) {
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

void LBM::compute_density_momentum_moment(double* f) {
	for(int i = 0; i < NX; i++) {
		for(int j = 0; j < NY; j++) {
			for(int k = 0; k < NZ; k++) {
				double new_density = 0;
				vector3<double> u;
				for(int w = 0; w < direction_size; w++) {
					new_density += f[scalar_index(i,j,k,w)];
					u.x += f[scalar_index(i,j,k,w)] * directions[w].x;
					u.y += f[scalar_index(i,j,k,w)] * directions[w].y;
					u.z += f[scalar_index(i,j,k,w)] * directions[w].z;
				}
				density_field[scalar_index(i,j,k)] = new_density;
				velocity_field[scalar_index(i,j,k)].x = u.x / new_density;
				velocity_field[scalar_index(i,j,k)].y = u.y / new_density;
				velocity_field[scalar_index(i,j,k)].z = u.z / new_density;
			}
		}
	}
}

double* LBM::stream() {
	double* f2 = new double[NX * NY * NZ * direction_size];
	for(int x = 0; x < NX; x++) {
		for(int y = 0; y < NY; y++) {
			for(int z = 0; z < NZ; z++) {
				for(int i = 0; i < direction_size; i++) {
					int xmd = (NX + x - (int)directions[i].x) % NX;
					int ymd = (NY + y - (int)directions[i].y) % NY;
					int zmd = (NZ + z - (int)directions[i].z) % NZ;
					f2[this->scalar_index(x,y,z,i)] = equilibrium_distribution[this->scalar_index(xmd,ymd,zmd,i)];
				}
			}
		}
	}
	return f2;
}

void LBM::collision(double* f) {//Performs the collision step.
	const double tauinv = 1;
	const double omtauinv = 1.0-tauinv;     // 1 - 1/tau
	for(int x = 0; x < NX; x++) {
		for(int y = 0; y < NY; y++) {
			for(int z = 0; z < NZ; z++) {
				for(int i = 0; i < direction_size; i++) {
					double dot_product = (double)velocity_field[scalar_index(x,y,z)].x * (double)directions[i].x + (double)velocity_field[scalar_index(x,y,z)].y * (double)directions[i].y +
						(double)velocity_field[scalar_index(x,y,z)].z * (double)directions[i].z;
					double norm_square = (double)velocity_field[scalar_index(x,y,z)].norm_square();
					double feq = (double)weights[i] * (double)density_field[scalar_index(x,y,z)] * (1.0 + 3.0 * dot_product + 4.5 * dot_product * dot_product - 1.5 * norm_square);
					f[scalar_index(x,y,z,i)] = omtauinv * equilibrium_distribution[scalar_index(x,y,z,i)] + tauinv * feq;
				}
			}
		}
	}
}

void LBM::perform_timestep() {
	double* f2 = this->stream();
	compute_density_momentum_moment(f2);
	collision(f2);
	delete[] equilibrium_distribution;
	equilibrium_distribution = f2;
}

void LBM::output_lbm_data(std::string filename) {
	std::ofstream output_stream;
  output_stream.open (filename, std::ofstream::out | std::ofstream::app);
	output_stream << "p,u_x,u_y,u_z" << '\n';
	for(int x = 0; x < NX; x++) {
		for(int y = 0; y < NY; y++) {
			for(int z = 0; z < NZ; z++) {
				output_stream << density_field[this->scalar_index(x,y,z)] << "," <<
					velocity_field[this->scalar_index(x,y,z)].x << "," <<
					velocity_field[this->scalar_index(x,y,z)].y << "," <<
					velocity_field[this->scalar_index(x,y,z)].z << '\n';
			}
		}
	}
	output_stream.close();
}
