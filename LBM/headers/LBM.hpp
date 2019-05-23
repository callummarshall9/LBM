#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include "LBM/headers/vector3.hpp"
#ifndef LBM_H_FILE
#define LBM_H_FILE

class LBM {

public:
	void initialise();
	LBM(int grid_size);
	void output_density();
	void output_velocity();
	void output_test();
	float compute_density_moment(int x, int y, int z);
	void compute_density_moment();
	void compute_momentum_density_moment(int x, int y, int z);
	void compute_momentum_density_moment();
	void stream();//Stream the current equilibrium distribution to the next distribution.
	void collision();//Perform the collision step. Assumes delta t / tau = 1.
	void perform_timestep();//Delta t = 1 lattice unit.
	void output_lbm_data(std::string filename);
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
	float weights[15] = { (2.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/9.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0),(1.0/72.0) };
	//This will result in a change in the equlibrium function which will be reflected below.
};

#endif
