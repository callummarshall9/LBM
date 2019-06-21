#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "vector3.hpp"
#ifndef LBM_H_FILE
#define LBM_H_FILE

class LBM {

public:
	LBM(int grid_size, std::string velocity_set, double c_s, double tau, std::string boundary_conditions, double gamma_dot);
	LBM(int nx, int ny, int nz, std::string velocity_set, double c_s, double tau, std::string boundary_condition, double gamma_dot);
	~LBM();
	void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z);//Set velocity at position in velocity field.
	void set_density(int x_field, int y_field, int z_field, double density);//Set density at position in density field.
	double calculate_feq(int i, int j, int k, int w);
    double calculate_feq(int i, int j, int k, int w, double u_le_x);
	void output_density();
	void output_velocity();
	void output_test();
	void output_indices_file();
	void compute_density_momentum_moment();
	void stream();//Stream the current equilibrium distribution to the next distribution.
	void collision();//Perform the collision step. Assumes delta t / tau = 1.
	void perform_timestep();//Delta t = 1 lattice unit.
	void output_lbm_data(std::string filename, bool header=true);
	void lookup_reverse();
    void output_f_array(double* f_array, int z_index);
	int get_time();
private:
    int time = 0;
	int NX;
	int NY;
	int NZ;
	double c_s;
	double nu;
	double tau;
	double* density_field;
	vector3<double>* velocity_field;
	double* previous_particle_distributions;
	double* particle_distributions;
	inline int scalar_index(int x, int y, int z) const {
        return (z * NX * NY) + (y * NX) + x;
    }
	inline int scalar_index(int x, int y, int z, int w) const {
        return (x + y * NX + z * NX * NY + w * NX * NY * NZ);
	}
	void output_array(double *array);
	void set_velocity_set(std::string velocity_set);//Used to internally generate velocity_set.
	void initialise();
	//Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
	int direction_size = 15;
	vector3<int>* directions;
	double* weights;
	int* reverse_indexes;
	std::string boundary_condition;
	double gamma_dot;
	std::string velocity_set;
	//This will result in a change in the equlibrium function which will be reflected below.
};

#endif
