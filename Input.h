#pragma once

#include "RadialWF.h"
#include "Potential.h"
#include "Grid.h"
#include <vector>
#include <string>
#include "IntegrateRateEquation.h"
#include <cmath>

using namespace std;

class Input
{
public:
	//Read or set the configuration. Assign all the quantum numbers and trial energies, lattice and potential
	Input(char* filename, vector<RadialWF> &Orbitals, Grid &Lattice, ofstream & log);
	Input(char* filename, vector<RadialWF> &Orbitals,  vector<RadialWF> &Virtual, Grid &Lattice, ofstream & log);
	Input(const Input & Other);

	string Pot_Model() { return model; }
	string Exited_Pot_Model() { return potential; }
	string Gauge() { return me_gauge; }
	string Name() { return name; }
	double Omega() { return omega; }
	double Width() { return width; }
	double Fluence() { return 10000*fluence; }
  void Set_Width(double ext_width) {width = ext_width;}
  void Set_Fluence(double ext_fluence) {fluence = ext_fluence;}
	
  void Set_Pulse(double ext_omega, double ext_fluence, double ext_width, int ext_T_size = 0) {
		omega = ext_omega;
		fluence = ext_fluence;
		width = ext_width;
		num_time_steps = ext_T_size;
	}
	int TimePts() { return num_time_steps; }
	int Hamiltonian();
	int Nuclear_Z() { return Z; }

	~Input();
protected:
	string name = "";// A name of argv[1] suppilied as filename without extension. Is added to output files.
	string model;
	string potential = "V_N";
	string me_gauge = "length";
	string hamiltonian = "LDA";
	double omega = 5000;// XFEL field frequency.
	double width = 5; // XFEL pulse width. Gaussian profile hardcoded.
	double fluence = 0; // XFEL pulse fluence.
	int num_time_steps = 0; // Guess number of time steps for time dynamics.
	int omp_threads = 1;
	int Z;

	bool write_charges = false;
	bool write_intensity = false;
	int out_time_steps = 500; // Guess number of time steps for time dynamics.

	double Master_tollerance = pow(10, -10);
	double No_exchange_tollerance = pow(10, -3);
	double HF_tollerance = pow(10, -6);
	int max_HF_iterations = 500;
	int max_Virt_iterations = 70;
};


struct AuxAtomData
{
  // Container for static atomic data later aggregated with time-dependent occupancies.
  vector<double> r_ion = vector<double>(0); // Ionic radius = \int dr r \rho(r, I) / \int r \rho(r, I), I - configuration.
  vector<vector<double>> form_fact_ion = vector<vector<double>>(0); // Time-dependent ionic form-factor. 
  vector<double> pol_ion = vector<double>(0); // Scalar polarizability.
};


class MolInp
{
	// Molecular input fpr coupled atom/electron plasma calcualtions.
public:
	MolInp(char* filename, ofstream & log);
	~MolInp() {}

	vector<Input> Atomic;
	vector<AtomRateData> Store; 
  vector<AuxAtomData> AuxStore;

	vector<Potential> Pots;
	vector<vector<RadialWF>> Orbits;
	vector<Grid> Latts;
	vector<vector<vector<int>>> Index;

	double Omega() {return omega;}
	double Width() {return width;}
	double Fluence() {return 10000*fluence;}
	int ini_T_size() {return num_time_steps;}
	double dropl_R() {return radius;}

  void Set_Fluence(double new_fluence) {fluence = new_fluence;}

  // Auxilary data calculation options.
  bool Calc_R_ion() {return calculate_r_ion;}
  bool Calc_Pol_ion() {return calculate_pol_ion;}
  bool Calc_FF_ion() {return calculate_form_fact_ion;}

	string name = "";
private:
	double omega = 5000;// XFEL photon energy, eV.
	double width = 5; // XFEL pulse width in femtoseconds. Gaussian profile hardcoded.
	double fluence = 0; // XFEL pulse fluence, 10^4 J/cm^2.
	int num_time_steps = 1000; // Guess number of time steps for time dynamics.
	double radius = 1000.;

	double unit_V = 1.;

  bool calculate_r_ion = false;
  bool calculate_pol_ion = false;
  bool calculate_form_fact_ion = false;
};

