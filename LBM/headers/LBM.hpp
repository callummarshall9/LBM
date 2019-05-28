#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include "vector3.hpp"
#ifndef LBM_H_FILE
#define LBM_H_FILE

class LBM {

public:
	LBM(int grid_size);
	LBM(int nx, int ny, int nz);
	void free_memory();
	void output_density();
	void output_velocity();
	void output_test();
	void compute_density_momentum_moment(double *f);
	double* stream();//Stream the current equilibrium distribution to the next distribution.
	void collision(double *f);//Perform the collision step. Assumes delta t / tau = 1.
	void perform_timestep();//Delta t = 1 lattice unit.
	void output_lbm_data(std::string filename);
private:
	int NX;
	int NY;
	int NZ;
	const double nu = 1.0/6.0;
	double* density_field;
	vector3<double>* velocity_field;
	double* equilibrium_distribution;
	int scalar_index(int x, int y, int z);
	int scalar_index(int x, int y, int z, int w);
	void output_array(double *array);
	void initialise();
	//Lattice directions using D3DQ15. assumed speed of sound c_s = 1/sqrt(3).
	static const int direction_size = 15;
	vector3<int> directions[direction_size] = {//c_i vectors.
		{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,1},{-1,-1,-1},{1,1,-1},{-1,-1,1},{1,-1,1},
		{-1,1,-1},{-1,1,1},{1,-1,-1}
	};//////////////////////1			2		3		4			5			6		7		8			9		10		11		12		13		14		15
	double weights[15] = { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0  };
	//This will result in a change in the equlibrium function which will be reflected below.
};

#endif
