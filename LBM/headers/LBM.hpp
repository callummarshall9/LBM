#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "vector3.hpp"
#ifndef LBM_H_FILE
#define LBM_H_FILE

class LBM {

public:
	LBM(int grid_size, std::string velocity_set, double c_s, double nu, double tau);
	LBM(int nx, int ny, int nz, std::string velocity_set, double c_s, double nu, double tau);
	~LBM();
	void set_velocity(int x_field, int y_field, int z_field, double u_x, double u_y, double u_z);//Set velocity at position in velocity field.
	void set_density(int x_field, int y_field, int z_field, double density);//Set density at position in density field.
	void output_density();
	void output_velocity();
	void output_test();
	void compute_density_momentum_moment();
	void stream();//Stream the current equilibrium distribution to the next distribution.
	void collision();//Perform the collision step. Assumes delta t / tau = 1.
	void perform_timestep();//Delta t = 1 lattice unit.
	void output_lbm_data(std::string filename, bool header=true);
private:
	int NX;
	int NY;
	int NZ;
	double c_s;
	double nu;
	double tau;
	double* density_field;
	vector3<double>* velocity_field;
	double* previous_equilibrium_distribution;
	double* equilibrium_distribution;
	int scalar_index(int x, int y, int z);
	int scalar_index(int x, int y, int z, int w);
	void output_array(double *array);
	void set_velocity_set(std::string velocity_set);//Used to internally generate velocity_set.
	void initialise();
	//Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
	int direction_size = 15;
	vector3<int>* directions;
	double* weights;
	//This will result in a change in the equlibrium function which will be reflected below.
};

#endif
