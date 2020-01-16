#include "stdafx.h"
#include "Input.h"
#include "Grid.h"
#include "Potential.h"
#include "RadialWF.h"
#include "DecayRates.h"
#include <fstream>
#include "HartreeFock.h"
#include "RateEquationSolver.h"
#include "Constant.h"
#include <sys/stat.h>
#include <ctime>
#include <eigen3/Eigen/Dense>
#include "DecayRates.h"
#include "Plasma.h"


using namespace std;

int main(int argc, char *argv[])
{
	string logname = "./output/log_";
	string filename;
	string tail;
	if (argc > 1) {
		filename = argv[1];
		size_t lastdot = filename.find_last_of(".");
		if (lastdot != std::string::npos) {
			tail = filename.substr(lastdot, filename.size());
			filename = filename.substr(0, lastdot);
		}
		size_t lastslash = filename.find_last_of("/");
		if (lastslash != std::string::npos) filename = filename.substr(lastslash+1);

		logname = logname + filename + ".txt";
	}
	else logname = logname + ".txt";
	ofstream log(logname);

	if (argc <= 1) {
		log << "no input file found. Exiting..." << endl;
		log.flush();
		log.close();
		return 1;
	}

	if (argc != 2)
	{
		log << "no input file found" << endl;

		return 1;
	}

	int c_start = clock();

	if (tail == ".txt") {

		// Molecular input. Default case.

		MolInp Init(argv[1], log);
		
		for (int a = 0; a < Init.Atomic.size(); a++) {
			printf("Nuclear charge: %d\n", Init.Atomic[a].Nuclear_Z());
			HartreeFock HF(Init.Latts[a], Init.Orbits[a], Init.Pots[a], Init.Atomic[a], log);

			RateEquationSolver Dynamics(Init.Latts[a], Init.Orbits[a], Init.Pots[a], Init.Atomic[a]);
			vector<int> final_occ(Init.Orbits[a].size(), 0);
			vector<int> max_occ(Init.Orbits[a].size(), 0);
			for (int i = 0; i < max_occ.size(); i++) {
				if (fabs(Init.Orbits[a][i].Energy) > Init.Omega()) final_occ[i] = Init.Orbits[a][i].occupancy();
				max_occ[i] = Init.Orbits[a][i].occupancy();
			}

			string name = Init.Store[a].name;
			double nAtoms = Init.Store[a].nAtoms;

			Init.Store[a] = Dynamics.SolvePlasmaBEB(max_occ, final_occ, log);
			Init.Store[a].name = name;
			Init.Store[a].nAtoms = nAtoms;
      Init.Store[a].R = Init.dropl_R();
			Init.Index[a] = Dynamics.Get_Indexes();
		}

		RateEquationSolver Dynamics(Init.Latts[0], Init.Orbits[0], Init.Pots[0], Init.Atomic[0]);

		Dynamics.SetupAndSolve(Init, log);
    
	} else {

		// Atomic input. 

		Grid Lattice(0);//dummy grid. Will be modified by configuration class
		vector<RadialWF> Orbitals;
		Input Init(argv[1], Orbitals, Lattice, log);
		
		std::cout << "Nuclear charge: " << Init.Nuclear_Z() << endl;
		std::cout << "Nuclear potential: pointlike Coulomb" << endl << endl;

		Potential U(&Lattice, Init.Nuclear_Z(), Init.Pot_Model());
		HartreeFock HF(Lattice, Orbitals, U, Init, log);
	
		if (Init.TimePts() != 0) {
			RateEquationSolver Dynamics(Lattice, Orbitals, U, Init);
			vector<int> final_occ(Orbitals.size(), 0);
			vector<int> max_occ(Orbitals.size(), 0);

			for (int i = 0; i < max_occ.size(); i++) {
				if (fabs(Orbitals[i].Energy) > Init.Omega()) final_occ[i] = Orbitals[i].occupancy();
				max_occ[i] = Orbitals[i].occupancy();
			}

			Dynamics.SolveFrozen(max_occ, final_occ, log);
			Dynamics.SetupAndSolve(log);
		}
	}


	int c_stop = clock();
	log << "====================================================" << endl
		<<"Total execution time " << float(c_stop - c_start)/1000000 << " sec" << endl;
	log.close();

	return 0;
}
