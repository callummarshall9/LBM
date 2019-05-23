#include <stdio.h>
#include <iostream>
#include "LBM/headers/LBM.hpp"
#include "LBM/headers/vector3.hpp"

extern char   __BUILD_DATE;
extern char   __BUILD_NUMBER;

int main(int argc, char** argv) {

	LBM solver(64);
	int scale = 1;
	int runs = 200 * scale * scale * scale;
	for(int i = 0; i < runs; i = i + 1) {
		float percentage = (float)(i+1) / (float)(runs) * 100.0;
		std::cout << "Saving data - " << (i+1)  << "/" << runs << " (" << percentage << "%)" << std::endl;
		solver.output_lbm_data("output/" + std::to_string(i) + ".csv");
		std::cout << "Performing timestep - " << i << "/" << runs << std::endl;
		solver.perform_timestep();
	}
}
