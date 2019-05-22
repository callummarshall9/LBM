#include <stdio.h>
#include <iostream>
#include "headers/LBM.hpp"
#include "headers/vector3.hpp"

extern char   __BUILD_DATE;
extern char   __BUILD_NUMBER;

int main(int argc, char** argv) {
	LBM solver(64);
	int scale = 1;
	for(int i = 0; i < 200 * scale * scale * scale; i = i + 1) {
		std::cout << "Performing timestep" << std::endl;
		solver.perform_timestep();
	}
}
