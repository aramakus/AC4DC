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
#include <string>
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

		// Molecular input.

		MolInp Init(argv[1], log);
		
		for (int a = 0; a < Init.Atomic.size(); a++) {
			printf("Nuclear charge: %d\n", Init.Atomic[a].Nuclear_Z());
			HartreeFock HF(Init.Latts[a], Init.Orbits[a], Init.Pots[a], Init.Atomic[a].Hamiltonian(), log);

			RateEquationSolver Dynamics(Init.Latts[a], Init.Orbits[a], Init.Pots[a], Init.Atomic[a]);
			vector<int> final_occ(Init.Orbits[a].size(), 0);
			vector<int> max_occ(Init.Orbits[a].size(), 0);
			for (int i = 0; i < max_occ.size(); i++) {
				if (fabs(Init.Orbits[a][i].Energy) > Init.Omega()) final_occ[i] = Init.Orbits[a][i].occupancy();
				max_occ[i] = Init.Orbits[a][i].occupancy();
			}

			string name = Init.Store[a].name;
			double nAtoms = Init.Store[a].nAtoms;
      vector<bool> args = {Init.Calc_R_ion(), Init.Calc_Pol_ion(), Init.Calc_FF_ion()}; // Save auxiliary data.

			Init.Store[a] = Dynamics.SolvePlasmaBEB(max_occ, final_occ, log, args);
			Init.Store[a].name = name;
			Init.Store[a].nAtoms = nAtoms;
      Init.Store[a].R = Init.dropl_R();
			Init.Index[a] = Dynamics.Get_Indexes();
      Init.AuxStore[a] = Dynamics.AtomAuxStore;
		}

		RateEquationSolver Dynamics(Init.Latts[0], Init.Orbits[0], Init.Pots[0], Init.Atomic[0]);

		Dynamics.SetupAndSolve(Init, log);
    
		} else {

		// Atomic input.

		Grid Lattice(0);//dummy grid. Will be modified by configuration class
		vector<RadialWF> Orbitals;
		vector<RadialWF> Virtual;
		Input Init(argv[1], Orbitals, Virtual, Lattice, log);
		
		std::cout << "Nuclear charge: " << Init.Nuclear_Z() << endl;
		std::cout << "Nuclear potential: pointlike Coulomb" << endl << endl;

		Potential U(&Lattice, Init.Nuclear_Z(), Init.Pot_Model());
		HartreeFock HF(Lattice, Orbitals, U, Init.Hamiltonian(), log);
	
		if (Init.TimePts() != 0) {
			RateEquationSolver Dynamics(Lattice, Orbitals, U, Init);
			vector<int> final_occ(Orbitals.size(), 0);
			vector<int> max_occ(Orbitals.size(), 0);

			for (int i = 0; i < max_occ.size(); i++) {
				if (fabs(Orbitals[i].Energy) > Init.Omega()) final_occ[i] = Orbitals[i].occupancy();
				max_occ[i] = Orbitals[i].occupancy();
			}

			if (Virtual.size() != 0) Dynamics.SolveFrozen(Virtual, max_occ, final_occ, log);
			else {
				Dynamics.SolveFrozen(max_occ, final_occ, log);//Dynamics.SolvePlasmaBEB(max_occ, final_occ, log);//
				//Dynamics.SetupAndSolve(log);
			}
			//vector<RateEquationSolver> Solutions({ Dynamics });
			// The problem is with time integration. Check the mesh and the nucmer of points.
		//	DamageMatrix Modes(Solutions, Lattice);
			/*
			RateEquationSolver Dynamics2(Lattice, Orbitals2, U2, Init2);
			vector<int> final_occ2(Orbitals2.size(), 0);
			vector<int> max_occ2(Orbitals2.size(), 0);
			for (int i = 0; i < max_occ2.size(); i++) {
				//if (fabs(Orbitals[i].Energy) > Init.Omega()) final_occ[i] =  Orbitals[i].occupancy();
				max_occ2[i] = Orbitals2[i].occupancy();
			}
			if (Virtual.size() != 0) Dynamics2.SolveFrozen(Virtual, max_occ2, final_occ2, log);
			else Dynamics2.SolveFrozen(max_occ2, final_occ2, log);*/
		} else {
			HF.LDA_Get_Virtual(Virtual, Orbitals, U, log);
			// Polarizability calculation. Average over configuration estimate.	
			polarize MixMe = HF.Hybrid(Orbitals, Virtual, U);
			MixMe.index = 0;
			vector<RadialWF*> OrbList(0);
			for (auto& Orb: Orbitals) OrbList.push_back(&Orb);
			for (auto& Vrt: Virtual) OrbList.push_back(&Vrt);
			
			double Pol = 0;
			int g = 0, e = 0;
			for (int i = 0; i < MixMe.excited.size(); i++) {
				g = MixMe.excited[i][0];
				e = MixMe.excited[i][1];
				double wght = 2*OrbList[g]->occupancy()*(1 - 1.*OrbList[e]->occupancy()/(4*OrbList[e]->L() + 2))/3/(2*OrbList[g]->L() + 1);	
				Pol += wght*MixMe.Dipoles[i]*MixMe.Dipoles[i]/(MixMe.extEnergy[i] - MixMe.refEnergy);
			}
			cout << "Polarizability estimate: " << Pol << endl;

			for (vector<RadialWF>::iterator Orb = Orbitals.begin(); Orb != Orbitals.end(); ++Orb) {
				printf("n = %d | l = %d | Energy = %2.8f\n", Orb->N(), Orb->L(), Orb->Energy);
			}

			DecayRates test(Lattice, Orbitals, U, Init);
			
		}

	}


	int c_stop = clock();
	log << "====================================================" << endl
		<<"Total execution time " << float(c_stop - c_start)/1000000 << " sec" << endl;
	log.close();

	system("pause");

	return 0;
}
