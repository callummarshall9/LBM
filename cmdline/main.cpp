#include <stdio.h>
#include <iostream>
#include "LBM/headers/LBM.hpp"
#include "LBM/headers/vector3.hpp"

int main(int argc, char** argv) {
	std::cout << "Do you want to clean the previous run? (1 - Yes, 0 - No): ";
	int choice;
	std::cin >> choice;
	if(choice == 1) {
		system("rm -rf /home/callummarshall/output");
		system("mkdir /home/callummarshall/output");
	}
	std::cout << "Enter grid size: ";
	int grid_size;
	std::cin >> grid_size;
	LBM solver(grid_size);
	solver.output_lbm_data("output/0.csv");
	int scale = 1;
	int runs = 1000 * scale * scale * scale;
	for(int i = 0; i < runs; i = i + 1) {
		solver.perform_timestep();
		if((i+1) % 200 == 0) {
			float percentage = (float)(i+1) / (float)(runs) * 100.0;
			std::cout << "Saving data - " << (i+1)  << "/" << runs << " (" << percentage << "%)" << std::endl;
			solver.output_lbm_data("output/" + std::to_string(i+1) + ".csv");
		}


	}
	solver.free_memory();
	return 0;
}
