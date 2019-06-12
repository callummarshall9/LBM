#include "LBM.hpp"
#include "vector3.cpp"
#include <iomanip>

LBM::LBM(int grid_size, std::string velocity_set, double m_c_s, double m_tau) : NX(grid_size),NY(grid_size),NZ(grid_size), c_s(m_c_s), tau(m_tau) {
	set_velocity_set(velocity_set);
	initialise();
}

LBM::LBM(int nx, int ny, int nz, std::string velocity_set, double m_c_s, double m_tau) : NX(nx), NY(ny), NZ(nz), c_s(m_c_s), tau(m_tau) {
	set_velocity_set(velocity_set);
	initialise();
}

void LBM::initialise() {
	int box_flatten_length = NX * NY * NZ;
	int equilibrium_flatten_length = box_flatten_length * direction_size;
	density_field = new double[box_flatten_length];
	velocity_field = new vector3<double>[box_flatten_length];
	previous_equilibrium_distribution = new double[equilibrium_flatten_length];
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
						//Equation 3.4 with c_s^2 = 1/3
					double feq = weights[w] * density_field[scalar_index(i,j,k)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - velocity_field[scalar_index(i,j,k)].norm_square() / (2 * c_s * c_s));
					previous_equilibrium_distribution[scalar_index(i,j,k,w)] = feq;
				}
			}
		}
	}
}

void LBM::set_velocity_set(std::string velocity_set) {
	if(velocity_set == "D3Q15") {
		direction_size = 15;
	 	directions = new vector3<int>[15]{ {0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},{1,-1,1},{-1,1,-1},{-1,1,1},{1,-1,-1} };
		weights = new double[15] { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0  };
	} else if(velocity_set == "D3Q27") {
		direction_size = 27;
		directions = new vector3<int>[27]{{0,0,0}, {1,0,0},{-1,0,0},
		{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,0,1},
		{-1,0,-1},{0,1,1},{0,-1,-1},{1,-1,0},{-1,1,0},{1,0,-1},{-1,0,1},
		{0,1,-1},{0,-1,1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},{1,-1,1},
		{-1,1,-1},{-1,1,1},{1,-1,-1}};
		weights = new double[27] {8.0/27.0,2.0/27.0,2.0/27.0,2.0/27.0,
			2.0/27.0,2.0/27.0,2.0/27.0, 1.0/54.0,1.0/54.0,1.0/54.0
		,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,1.0/54.0,
	    1.0/54.0,1.0/54.0,1.0/54.0, 1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0
        ,1.0/216.0, 1.0/216.0, 1.0/216.0, 1.0/216.0};
	}
}

LBM::~LBM() {
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

void LBM::set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z) {
	velocity_field[scalar_index(x_field,y_field,z_field)].x = u_x;
	velocity_field[scalar_index(x_field,y_field,z_field)].y = u_y;
	velocity_field[scalar_index(x_field,y_field,z_field)].z = u_z;
}

void LBM::set_density(int x_field, int y_field, int z_field, double density) {
	density_field[scalar_index(x_field,y_field,z_field)] = density;
}


void LBM::compute_density_momentum_moment() {
	for(int i = 0; i < NX; i++) {
		for(int j = 0; j < NY; j++) {
			for(int k = 0; k < NZ; k++) {
				//Equation 3.1
				double new_density = 0;
				vector3<double> u;
				for(int w = 0; w < direction_size; w++) {
					new_density += equilibrium_distribution[scalar_index(i,j,k,w)];
					u.x += equilibrium_distribution[scalar_index(i,j,k,w)] * directions[w].x;
					u.y += equilibrium_distribution[scalar_index(i,j,k,w)] * directions[w].y;
					u.z += equilibrium_distribution[scalar_index(i,j,k,w)] * directions[w].z;
				}
				density_field[scalar_index(i,j,k)] = new_density;
				velocity_field[scalar_index(i,j,k)].x = u.x / new_density;
				velocity_field[scalar_index(i,j,k)].y = u.y / new_density;
				velocity_field[scalar_index(i,j,k)].z = u.z / new_density;
			}
		}
	}
}

void LBM::stream() {
	for(int x = 0; x < NX; x++) {
		for(int y = 0; y < NY; y++) {
			for(int z = 0; z < NZ; z++) {
				for(int i = 0; i < direction_size; i++) {
					int xmd = (NX + x - (int)directions[i].x) % NX;
					int ymd = (NY + y - (int)directions[i].y) % NY;
					int zmd = (NZ + z - (int)directions[i].z) % NZ;
					//Equation 3.10 with periodic boundary conditions.
					equilibrium_distribution[this->scalar_index(x,y,z,i)] = previous_equilibrium_distribution[this->scalar_index(xmd,ymd,zmd,i)];
				}
			}
		}
	}
}

void LBM::collision() {//Performs the collision step.
	const double tauinv = 1.0 / tau;
	const double omtauinv = 1.0-tauinv;     // 1 - 1/tau
	for(int x = 0; x < NX; x++) {
		for(int y = 0; y < NY; y++) {
			for(int z = 0; z < NZ; z++) {
				for(int i = 0; i < direction_size; i++) {
					double dot_product = (double)velocity_field[scalar_index(x,y,z)].x * (double)directions[i].x + (double)velocity_field[scalar_index(x,y,z)].y * (double)directions[i].y +
						(double)velocity_field[scalar_index(x,y,z)].z * (double)directions[i].z;
					double norm_square = (double)velocity_field[scalar_index(x,y,z)].norm_square();
					double feq = (double)weights[i] * (double)density_field[scalar_index(x,y,z)] * (1.0 + dot_product / (c_s * c_s) + dot_product * dot_product / (2 * c_s * c_s * c_s * c_s) - norm_square / (2 * c_s * c_s));
					//Equation 3.9
					equilibrium_distribution[scalar_index(x,y,z,i)] = omtauinv * previous_equilibrium_distribution[scalar_index(x,y,z,i)] + tauinv * feq;
				}
			}
		}
	}
}

void LBM::perform_timestep() {
	this->stream();
	compute_density_momentum_moment();
	collision();
	previous_equilibrium_distribution = equilibrium_distribution;
}


void LBM::output_lbm_data(std::string filename, bool header) {
	std::ofstream output_stream;
  output_stream.open (filename, std::ofstream::out | std::ofstream::app);
	if(header) {
		output_stream << "p,u_x,u_y,u_z" << '\n';
	}
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
