/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/
#include "Input.h"
#include "Constant.h"
#include <fstream>
#include <sstream>
#include <map>

Input::Input(char *filename, vector<RadialWF> &Orbitals, Grid &Lattice, ofstream & log)
{
	log << "Input file: " << filename << endl;

	name = filename;
	size_t lastdot = name.find_last_of(".");
	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
	size_t lastslash = name.find_last_of("/");
	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

	ifstream infile(filename);

	map<string, vector<string>> FileContent;
	string comment = "//";
	string curr_key;
	
	while (!infile.eof() && infile.is_open()) {
		string line;
		getline(infile, line);

		if (!line.compare(0, 2, comment)) continue;
		if (!line.compare(0, 2, "")) continue;
		if (!line.compare(0, 1, "#")) {
			if ( FileContent.find(line) == FileContent.end() ) {
				FileContent[line] = vector<string>(0);
			}
			curr_key = line;
		} else {
			FileContent[curr_key].push_back(line);
		}
	}

	for (int n = 0; n < FileContent["#PULSE"].size(); n++) {
		stringstream stream(FileContent["#PULSE"][n]);
		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> fluence;
	}

	for (int n = 0; n < FileContent["#OUTPUT"].size(); n++) {
		stringstream stream(FileContent["#OUTPUT"][n]);
		if (n == 0) stream >> out_time_steps;
		if (n == 1) {
			char tmp;
			stream >> tmp;
			if (tmp == 'Y') write_charges = true;
		}
		if (n == 2) {
			char tmp;
			stream >> tmp;
			if (tmp == 'Y') write_intensity = true;
		}
	}

	int line_indentifier = 0, num_grid_pts, num_orbitals, N, L, current_orbital = 0, occupancy;
	double r_min, r_box;
	map<char, int> angular = {{'s', 0}, {'p', 1}, {'d', 2}, {'f', 3}};

	for (int n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
		stringstream stream(FileContent["#NUMERICAL"][n]);
		if (n == 0) stream >> num_grid_pts;
		if (n == 1) stream >> r_min;
		if (n == 2) stream >> r_box;
		if (n == 3) stream >> num_time_steps;
		if (n == 4) stream >> omp_threads;
		if (n == 5) {
			int tmp = 0;
			stream >> tmp;
			HF_tollerance = pow(10, tmp);
		}
		if (n == 6) {
			int tmp = 0;
			stream >> tmp;
			Master_tollerance = pow(10, tmp);
		}
		if (n == 7) {
			int tmp = 0;
			stream >> tmp;
			No_exchange_tollerance = pow(10, tmp);
		}
		if (n == 8) stream >> max_HF_iterations;
	}

	for (int n = 0; n < FileContent["#ATOM"].size(); n++) {
		stringstream stream(FileContent["#ATOM"][n]);
		if (n == 0) stream >> Z;
		if (n == 1) stream >> model;
		if (n == 2) stream >> hamiltonian;
		if (n == 3) stream >> num_orbitals;
		if (n > 3 && n < num_orbitals + 4) {
			char tmp;
			stream >> N >> tmp >> occupancy;
			if (occupancy == 0) log << "Orbital with N=" << N << ", L=" << angular[tmp] << " has occupancy=0!" << endl;
				Orbitals.push_back(RadialWF(num_grid_pts));
				Orbitals.back().set_N(N);
				Orbitals.back().set_L(angular[tmp]);//setting L overwrites occupancy with 4L+2
				Orbitals.back().set_occupancy(occupancy);
		}
		if (n == num_orbitals + 4) stream >> potential;
		if (n == num_orbitals + 5) stream >> me_gauge;
	}

	Grid lattice(num_grid_pts, r_min/Z, r_box, 4);
	Lattice = lattice;
	fluence /= omega/Constant::eV_in_au;
}

Input::Input(const Input & Other)
{
	name = Other.name;
	model = Other.model;
	potential = Other.potential;
	me_gauge = Other.me_gauge;
	hamiltonian = Other.hamiltonian;
	omega = Other.omega;
	width = Other.width;
	fluence = Other.fluence;
	num_time_steps = Other.num_time_steps;
	omp_threads = Other.omp_threads;
	Z = Other.Z;
	write_charges = Other.write_charges;
	write_intensity = Other.write_intensity;
	out_time_steps = Other.out_time_steps;
	Master_tollerance = Other.Master_tollerance;
	No_exchange_tollerance = Other.No_exchange_tollerance;
	HF_tollerance = Other.HF_tollerance;
	max_HF_iterations = Other.max_HF_iterations;
}

int Input::Hamiltonian()
{
	if (hamiltonian == "HF") return 0;
	else return 1;
}

Input::~Input()
{
}

MolInp::MolInp(char* filename, ofstream & log)
{
	// Input file for molecular ionization calculation.
	map<string, vector<string>> FileContent;

	name = filename;
	size_t lastdot = name.find_last_of(".");
	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
	size_t lastslash = name.find_last_of("/");
	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

	ifstream infile(filename);
	string comment = "//";
	string curr_key = "";

	while (!infile.eof())
	{
		string line;
		getline(infile, line);
		if (!line.compare(0, 2, comment)) continue;
		if (!line.compare(0, 1, "")) continue;
		if (!line.compare(0, 1, "#")) {
			if ( FileContent.find(line) == FileContent.end() ) {
				FileContent[line] = vector<string>(0);
			}
			curr_key = line;
		} else {
			FileContent[curr_key].push_back(line);
		}
	}

	int num_atoms = FileContent["#ATOMS"].size();

	Orbits.clear();
	Orbits.resize(num_atoms);
	Latts.clear();
	Latts.resize(num_atoms, Grid(0));
	Pots.clear();
	Pots.resize(num_atoms);
	Atomic.clear();
	Store.clear();
	Store.resize(num_atoms);
	Index.clear();
	Index.resize(num_atoms);

	for (int n = 0; n < FileContent["#VOLUME"].size(); n++) {
		stringstream stream(FileContent["#VOLUME"][n]);

		if (n == 0) stream >> unit_V;
		if (n == 1) stream >> radius;
	}

	for (int n = 0; n < FileContent["#OUTPUT"].size(); n++) {
		stringstream stream(FileContent["#OUTPUT"][n]);
		char tmp;

		if (n == 0) stream >> out_T_size;
		if (n == 1) {
			stream >> tmp;
			if (tmp == 'Y') write_charges = true;
		}
		if (n == 2) {
			stream >> tmp;
			if (tmp == 'Y') write_intensity = true;
		}
		if (n == 3) {
			stream >> tmp;
			if (tmp != 'Y') write_md_data = false;
		}
	}

	for (int n = 0; n < FileContent["#PULSE"].size(); n++) {
		stringstream stream(FileContent["#PULSE"][n]);

		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> fluence;
		if (n == 3) stream >> num_time_steps;
	}

	for (int n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
		stringstream stream(FileContent["#NUMERICAL"][n]);

		if (n == 0) stream >> num_time_steps;
		if (n == 1) stream >> omp_threads;
	}

  // Convert to number of photon flux.
  fluence /= omega/Constant::eV_in_au;
	radius /= Constant::au_in_Angs;
	unit_V /= Constant::au_in_Angs*Constant::au_in_Angs*Constant::au_in_Angs; 
	
	for (int i = 0; i < num_atoms; i++) {
		string at_name;
		double at_num;

		stringstream stream(FileContent["#ATOMS"][i]);
		stream >> at_name >> at_num;

		Store[i].nAtoms = at_num/unit_V;
		Store[i].name = at_name;
		Store[i].R = radius;

		at_name = "input/" + at_name + ".inp";

		Atomic.push_back(Input((char*)at_name.c_str(), Orbits[i], Latts[i], log));
		Atomic.back().Set_Pulse(omega, fluence, width);
		Atomic.back().Set_Num_Threads(omp_threads);

		Potential U(&Latts[i], Atomic[i].Nuclear_Z(), Atomic[i].Pot_Model());

		Pots[i] = U;
	}


}
