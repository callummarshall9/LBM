#include <stdio.h>
#include <iostream>
#include <fstream>
#include <streambuf>
//Now Linux only.
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "LBM.hpp"
#include "vector3.hpp"
//RapidJSON files.
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"
#include "csv.h"

using namespace rapidjson;

inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char** argv) {
  if(!file_exists("options.json")) {
    std::cout << "Please ensure that options.json exists. If not, it can be obtained from root directory of GitHub repo." << std::endl;
    return -1;
  }
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
	std::cout << "Save every: " << save_every << std::endl;
	auto grid_size = d["grid_size"].GetArray();
	int NX = grid_size[0].GetInt();
	int NY = grid_size[1].GetInt();
	int NZ = grid_size[2].GetInt();
	std::cout << "Grid size: " << NX << "x" << NY << "x" << NZ << std::endl;
	Value& m_c_s = d["c_s"];
	double c_s = m_c_s.GetDouble();
	std::cout << "c_s (Speed of sound): " << c_s << std::endl;
	Value& m_tau = d["tau"];
	double tau = m_tau.GetDouble();
	std::cout << "tau: " << tau << std::endl;
	Value& m_velocity_set = d["velocity_set"];
	std::string velocity_set = m_velocity_set.GetString();
	Value& m_boundary_conditions = d["boundary_conditions"];
	std::string boundary_conditions = m_boundary_conditions.GetString();
	if(velocity_set != "D3Q15" && velocity_set != "D3Q27" && velocity_set != "D2Q9") {
		std::cout << "Error: Please specify a valid velocity set such as D3Q15,D3Q27 or D2Q9." << std::endl;
		return -1;
	}
	std::cout << "Velocity set: " << velocity_set << std::endl;
	if(boundary_conditions != "lees_edwards" && boundary_conditions != "periodic" && boundary_conditions != "couette") {
	    std::cout << "Errors: boundary_conditions in options.json can either be: periodic, Couette (D2Q9 only) or lees_edwards (Lees-Edwards Shear, Please see research paper by Alexander Wagner)";
	    return -1;
	}
	std::cout << "Boundary conditions: " << boundary_conditions << std::endl;
	    Value& m_gamma_dot = d["gamma_dot"];
	    double gamma_dot = m_gamma_dot.GetDouble();
	std::cout << "Shear rate (gamma_dot): " << gamma_dot << std::endl;
	if(NZ != 1 && velocity_set == "D2Q9") {
	    std::cout << "Warning: NZ=1 for D2Q9.";
	    return -1;
	}
    Value& m_n_steps = d["n_steps"];
    double nsteps = m_n_steps.GetInt();
	LBM *solver = new LBM(NX,NY,NZ, velocity_set, c_s, tau, boundary_conditions, gamma_dot);
	for(int i = 0; i < argc; i++) {
		if(std::string(argv[i]) == "generate_ic") {
			solver->output_lbm_data("ic.csv", true);
			std::cout << "Generated ic.csv" << std::endl;
			return 0;
		}
	}
	if(file_exists("ic.csv")) {
		std::cout << "Loading initial conditions" << std::endl;
		io::CSVReader<4> in("ic.csv");
		in.read_header(io::ignore_extra_column, "p","u_x","u_y","u_z");
		double density,u_x,u_y,u_z;
		for(int i = 0; i < NX; i++) {
			for(int j = 0; j < NY; j++) {
				for(int k = 0; k < NZ; k++) {
					in.read_row(density,u_x,u_y,u_z);
					solver->set_density(i,j,k,density);
					solver->set_velocity(i,j,k,u_x,u_y,u_z);
				}
			}
		}
		std::cout << "Loaded initial conditions" << std::endl;
	} else {
		std::cout << "Using default of p=1 for all x,y,z and u(x,t=0)=0 for all x,y,z. (Steady state)" << std::endl;
		std::cout << "If you wish to use your own initial conditions, please run the program but with command: generate_ic as a argument which will output ic.csv in format of p,u_x,u_y,u_z, assume indexes are incrementing i,j,k for i<NX,j<NY and k<NZ" << std::endl;
	}

	solver->output_lbm_data("output/0.csv");
	int scale = 1;
	int runs = 1000 * scale * scale * scale;
	for(int i = 0; i < runs; i = i + 1) {
		solver->perform_timestep();
		if((i+1) % save_every == 0) {
            float percentage = (float) (i + 1) / (float) (runs) * 100.0;
            std::cout << "Saving data - " << (i + 1) << "/" << runs << " (" << percentage << "%)" << std::endl;
            solver->output_lbm_data("output/" + std::to_string(i + 1) + ".csv");
            solver->output_velocity();
        }
	}
	delete solver;
	return 0;
}
