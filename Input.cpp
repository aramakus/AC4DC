#include "Input.h"
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
	string comment = "//";
	int line_indentifier = 0, num_grid_pts, num_orbitals, n, current_orbital = 0, occupancy;
	double r_min, r_box;
	char angular;

	Orbitals.clear();
	int line_num = 0;

	while (!infile.eof())
	{
		string line;
		getline(infile, line);
		line_num++;
		if (line.compare(0, 2, comment))//skip if comment which starts with "//"
		{
			stringstream stream(line);
			switch (line_indentifier)
			{
			case 0:
				stream >> hamiltonian;
				line_indentifier++;
				break;
			case 1:
				stream >> Z;
				line_indentifier++;
				break;
			case 2:
				stream >> model;
				if (model.compare("sphere")) stream >> r_box;
				line_indentifier++;
				break;
			case 3:
				stream >> num_grid_pts;
				line_indentifier++;
				break;
			case 4:
				stream >> r_min;
				line_indentifier++;
				break;
			case 5:
				stream >> r_box;
				line_indentifier++;
				break;
			case 6:
				stream >> num_orbitals;
				line_indentifier++;
				break;
			case 7:
				if (current_orbital < num_orbitals)
				{
					current_orbital++;
					Orbitals.push_back(RadialWF(0));
					stream >> n >> angular >> occupancy;
					Orbitals.back().set_N(n);
					if (angular == 's') n = 0;
					if (angular == 'p') n = 1;
					if (angular == 'd') n = 2;
					if (angular == 'f') n = 3;
					Orbitals.back().set_L(n);//setting L overwrites occupancy with 4L+2
					Orbitals.back().set_occupancy(occupancy);
				}
				if (current_orbital == num_orbitals) line_indentifier++;
				break;
			case 8:
				stream >> potential;
				line_indentifier++;
				break;
			case 9:
				stream >> me_gauge;
				line_indentifier++;
				break;
			case 10:
				stream >> omega;
				line_indentifier++;
				break;
			case 11:
				stream >> width;
				line_indentifier++;
				break;
			case 12:
				stream >> fluence;
				line_indentifier++;
				break;
			case 13:
				stream >> num_time_steps;
				line_indentifier++;
				break;
			case 14:
				break;
			default:
				log << "In Configuration: extra line or comment sign missing in line" << line_num << endl;
			}

		}

	}
	infile.close();

	Grid lattice(num_grid_pts, r_min/Z, r_box, "exponential");
	Lattice = lattice;

	for (int i = 0; i < Orbitals.size(); i++) Orbitals[i].resize(lattice.size());
}

Input::Input(char *filename, vector<RadialWF> &Orbitals, vector<RadialWF> &Virtual, Grid &Lattice, ofstream & log)
{
	log << "Input file: " << filename << endl;

	name = filename;
	size_t lastdot = name.find_last_of(".");
	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
	size_t lastslash = name.find_last_of("/");
	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

	ifstream infile(filename);
	string comment = "//";
	int line_indentifier = 0, num_grid_pts, num_orbitals, n, l, current_orbital = 0, occupancy;
	double r_min, r_box;
	char angular;

	Orbitals.clear();
	int line_num = 0;

	log << "====================================================" << endl;
	log << "Input file copy: " << endl;
	log << "====================================================" << endl;
	while (!infile.eof())
	{
		string line;
		getline(infile, line);
		log << line << endl;
		line_num++;
		if (line.compare(0, 2, comment))//skip if comment which starts with "//"
		{
			stringstream stream(line);
			switch (line_indentifier)
			{
			case 0:
				stream >> hamiltonian;
				line_indentifier++;
				break;
			case 1:
				stream >> Z;
				line_indentifier++;
				break;
			case 2:
				stream >> model;
				if (model.compare("sphere")) stream >> r_box;
				line_indentifier++;
				break;
			case 3:
				stream >> num_grid_pts;
				line_indentifier++;
				break;
			case 4:
				stream >> r_min;
				line_indentifier++;
				break;
			case 5:
				stream >> r_box;
				line_indentifier++;
				break;
			case 6:
				stream >> num_orbitals;
				line_indentifier++;
				break;
			case 7:
				if (current_orbital < num_orbitals)
				{
					current_orbital++;
					stream >> n >> angular >> occupancy;
					if (angular == 's') l = 0;
					if (angular == 'p') l = 1;
					if (angular == 'd') l = 2;
					if (angular == 'f') l = 3;
					if (occupancy != 0) {
						Orbitals.push_back(RadialWF(num_grid_pts));
						Orbitals.back().set_N(n);
						Orbitals.back().set_L(l);//setting L overwrites occupancy with 4L+2
						Orbitals.back().set_occupancy(occupancy);
					}
					else {
						Virtual.push_back(RadialWF(num_grid_pts));
						Virtual.back().set_N(n);
						Virtual.back().set_L(l);//setting L overwrites occupancy with 4L+2
						Virtual.back().set_occupancy(occupancy);
					}
				}
				if (current_orbital == num_orbitals) line_indentifier++;
				break;
			case 8:
				stream >> potential;
				line_indentifier++;
				break;
			case 9:
				stream >> me_gauge;
				line_indentifier++;
				break;
			case 10:
				stream >> omega;
				line_indentifier++;
				break;
			case 11:
				stream >> width;
				line_indentifier++;
				break;
			case 12:
				stream >> fluence;
				line_indentifier++;
				break;
			case 13:
				stream >> num_time_steps;
				line_indentifier++;
				break;
			case 14:
				break;
			default:
				log << "In Configuration: extra line or comment sign missing in line" << line_num << endl;
			}

		}

	}
	infile.close();
	
	log << "====================================================" << endl;

	//Grid lattice(num_grid_pts, r_min/Z, r_box, "exponential");
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
	Z = Other.Z;
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
		if (!line.compare(0, 1, "#")) {
			if ( FileContent.find(line) == FileContent.end() ) {
				FileContent[line] = vector<string>(0);
			}
			curr_key = line;
		} else {
			FileContent[curr_key].push_back(line);
		}
	}

	int num_atoms = FileContent["#ATOM"].size();

	Orbits.clear();
	Orbits.resize(num_atoms);
	Latts.clear();
	Latts.resize(num_atoms, Grid(0));
	Pots.clear();
	Pots.resize(num_atoms);
	Atomic.clear();
	Store.clear();
	Store.resize(num_atoms);
  AuxStore.clear();
	AuxStore.resize(num_atoms);
	Index.clear();
	Index.resize(num_atoms);

	for (int n = 0; n < FileContent["#MOLECULE"].size(); n++) {
		stringstream stream(FileContent["#MOLECULE"][n]);
    string line_key;
    stream >> line_key;

		if (n == 0 && line_key == "Y") calculate_r_ion = true;
		if (n == 1 && line_key == "Y") calculate_pol_ion = true;
    if (n == 2 && line_key == "Y") calculate_form_fact_ion = true;
	}

	for (int n = 0; n < FileContent["#VOLUME"].size(); n++) {
		stringstream stream(FileContent["#VOLUME"][n]);

		if (n == 0) stream >> unit_V;
		if (n == 1) stream >> radius;
	}

	for (int n = 0; n < FileContent["#PULSE"].size(); n++) {
		stringstream stream(FileContent["#PULSE"][n]);

		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> fluence;
		if (n == 3) stream >> num_time_steps;
	}
  // Convert to number of photon flux.
  fluence /= omega/Constant::eV_in_au;

	for (int i = 0; i < num_atoms; i++) {
		string at_name;
		double at_num;

		stringstream stream(FileContent["#ATOM"][i]);
		stream >> at_name >> at_num;

		Store[i].nAtoms = at_num/unit_V;
		Store[i].name = at_name;
		Store[i].R = radius;

		at_name = "input/" + at_name + ".inp";

		vector<RadialWF> Virtual;
		Atomic.push_back(Input((char*)at_name.c_str(), Orbits[i], Virtual, Latts[i], log) );
		Atomic.back().Set_Pulse(omega, fluence, width);

		Potential U(&Latts[i], Atomic[i].Nuclear_Z(), Atomic[i].Pot_Model());

		Pots[i] = U;
	}


}
