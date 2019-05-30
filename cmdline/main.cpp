#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
#include "LBM.hpp"
#include "vector3.hpp"
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"

using namespace rapidjson;

int main(int argc, char** argv) {
	std::cout << "Do you want to clean the previous run? (1 - Yes, 0 - No): ";
	int choice;
	std::cin >> choice;
	if(choice == 1) {
		system("rm -rf output");
		system("mkdir output");
	}
	std::ifstream t("options.json");
	std::string str((std::istreambuf_iterator<char>(t)),
	                 std::istreambuf_iterator<char>());
	rapidjson::Document d;
	d.Parse(str.c_str());
	Value& save_every_value = d["save_every"];
	int save_every = save_every_value.GetInt();
	std::cout << save_every << std::endl;
	auto grid_size = d["grid_size"].GetArray();
	int NX = grid_size[0].GetInt();
	int NY = grid_size[1].GetInt();
	int NZ = grid_size[2].GetInt();
	LBM solver(NX,NY,NZ);
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
