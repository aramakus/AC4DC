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
#include "RateEquationSolver.h"
#include "HartreeFock.h"
#include "DecayRates.h"
#include "Numerics.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include <sstream>
#include <string>
#include <omp.h>
#include <algorithm>
#include "EigenSolver.h"
#include "Plasma.h"
#include <utility>

inline bool exists_test(const std::string&);
vector<double> generate_dT(int);
vector<double> generate_T(vector<double>&);
vector<double> generate_I(vector<double>&, double, double);
void SmoothOrigin(vector<double>&, vector<double>&);

using namespace CustomDataType;

RateEquationSolver::RateEquationSolver(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, Input & Inp) :
lattice(Lattice), orbitals(Orbitals), u(U), input(Inp)
{
  omp_set_num_threads(input.Num_Threads());
//Orbitals are HF wavefunctions. This configuration is an initial state.
//Assuming there are no unoccupied states in initial configuration!!!
}


int RateEquationSolver::SolveFrozen(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	// Solves system of rate equations exactly.
	// Final_occ defines the lowest possible occupancies for the initiall orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return 1;

	cout << "Check if there are pre-calculated rates..." << endl;
    // Make an output folder, if it doesn't exist yet
	if (!exists_test("./output")) {
		mkdir("output", ACCESSPERMS);
	}

	// Make a data folder inside output
	if (!exists_test("./output/" + input.Name())) {
		string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}
    // make a folder for the Rates
    string RateLocation = "./output/" + input.Name() + "/Rates/";
	if (!exists_test("./output/" + input.Name())) {
		mkdir(RateLocation.c_str(), ACCESSPERMS);
	}

	bool existPht = !ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool existFlr = !ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool existAug = !ReadRates(RateLocation + "Auger.txt", Store.Auger);
	bool existFF = !exists_test(RateLocation + "Form_Factor.txt");

	if (existPht) printf("Photoionization rates found. Reading...\n");
	if (existFlr) printf("Fluorescence rates found. Reading...\n");
	if (existAug) printf("Auger rates found. Reading...\n");
    if (existFF) printf("Auger rates found. Reading...\n");

	string PolarFileName = "./output/Polar_" + input.Name() + ".txt";

	{
		cout << "No rates found. Calculating..." << endl;
		cout << "Total number of configurations: " << dimension << endl;
		Rate Tmp;
		vector<Rate> LocalPhoto(0);
		vector<Rate> LocalFluor(0);
		vector<Rate> LocalAuger(0);
		vector<ffactor> LocalFF(0);

		#pragma omp parallel default(none) \
		shared(cout, runlog, existAug, existFlr, existPht, existFF) \
		private(Tmp, Max_occ, LocalPhoto, LocalAuger, LocalFluor, LocalFF)
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0; i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				vector<RadialWF> Orbitals = orbitals;
				cout << "configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (int j = 0; j < Orbitals.size(); j++)
				{
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
				}
				Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				Potential U(&Lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(Lattice, Orbitals, U, input, runlog);

				DecayRates Transit(Lattice, Orbitals, u, input);

				// ======= Experimental =========
				if (existFF) LocalFF.push_back({i, Transit.FT_density()});

				Tmp.from = i;

				if (existPht) {
					vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega()/Constant::eV_in_au, runlog);
					for (int k = 0; k < PhotoIon.size(); k++)
					{
						if (PhotoIon[k].val <= 0) continue;
						Tmp.val = PhotoIon[k].val;
						Tmp.to = i + hole_posit[PhotoIon[k].hole];
						Tmp.energy = input.Omega()/Constant::eV_in_au - Orbitals[PhotoIon[k].hole].Energy;
						LocalPhoto.push_back(Tmp);
					}
				}

				if (i != 0)
				{
					if (existFlr) {
						vector<fluor> Fluor = Transit.Fluor();
						for (int k = 0; k < Fluor.size(); k++)
						{
							if (Fluor[k].val <= 0) continue;
							Tmp.val = Fluor[k].val;
							Tmp.to = i - hole_posit[Fluor[k].hole] + hole_posit[Fluor[k].fill];
							Tmp.energy = Orbitals[Fluor[k].fill].Energy - Orbitals[Fluor[k].hole].Energy;
							LocalFluor.push_back(Tmp);
						}
					}

					if (existAug) {
						vector<auger> Auger = Transit.Auger(Max_occ, runlog);
						for (int k = 0; k < Auger.size(); k++)
						{
							if (Auger[k].val <= 0) continue;
							Tmp.val = Auger[k].val;
							Tmp.to = i - hole_posit[Auger[k].hole] + hole_posit[Auger[k].fill] + hole_posit[Auger[k].eject];
							Tmp.energy = Auger[k].energy;
							LocalAuger.push_back(Tmp);
						}
					}
				}
			}

			#pragma omp critical
			{
				Store.Photo.insert(Store.Photo.end(), LocalPhoto.begin(), LocalPhoto.end());
				Store.Fluor.insert(Store.Fluor.end(), LocalFluor.begin(), LocalFluor.end());
				Store.Auger.insert(Store.Auger.end(), LocalAuger.begin(), LocalAuger.end());
				FF.insert(FF.end(), LocalFF.begin(), LocalFF.end());
			}
		}

        // Add the last configuration
        if (existFF) FF.push_back({dimension-1, vector<double> (FF.back().val.size(), 0)});

		sort(Store.Photo.begin(), Store.Photo.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(FF.begin(), FF.end(), [](ffactor A, ffactor B) { return (A.index < B.index); });
		GenerateRateKeys(Store.Auger);

		if (existPht) {
			string dummy = RateLocation + "Photo.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Photo) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existFlr) {
			string dummy = RateLocation + "Fluor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Fluor) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existPht) {
			string dummy = RateLocation + "Auger.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Auger) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existFF) {
		string dummy = RateLocation + "Form_Factor.txt";
		FILE * fl = fopen(dummy.c_str(), "w");
		for (auto& ff : FF) {
			for (int i = 0; i < ff.val.size(); i++) fprintf(fl, "%3.5f ", ff.val[i]);
			fprintf(fl, "\n");
		}
		fclose(fl);
	}
	}

	string IndexTrslt = "./output/" + input.Name() + "/index.txt";
	ofstream config_out(IndexTrslt);
	for (int i = 0; i < Index.size(); i++) {
		for (int j = 0; j < Max_occ.size(); j++) {
			config_out << Max_occ[j] - Index[i][j] << " ";
		}
		config_out << endl;
	}

 	return dimension;
}


AtomRateData RateEquationSolver::SolvePlasmaBEB(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	// Solves system of Atomic rate equations + Temperature and electron concentration in Plasma exactly.
	// Final_occ defines the lowest possible occupancies for the initiall orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

  // Example args:    vector<bool> args = {Init.Calc_R_ion(), Init.Calc_Pol_ion(), Init.Calc_FF_ion()};

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return Store;

	Store.num_conf = dimension;

	cout << "Check if there are pre-calculated rates..." << endl;
	string RateLocation = "./output/" + input.Name() + "/Rates/";
	if (!exists_test("./output/" + input.Name())) {
		string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}
	if (!exists_test(RateLocation)) {
		string dirstring = "output/" + input.Name() + "/Rates";
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	bool existPht = !ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool existFlr = !ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool existAug = !ReadRates(RateLocation + "Auger.txt", Store.Auger);
	bool existFF = !ReadFFactors(RateLocation + "Form_Factor.txt", Store.FF);

	if (existPht) printf("Photoionization rates found. Reading...\n");
	if (existFlr) printf("Fluorescence rates found. Reading...\n");
	if (existAug) printf("Auger rates found. Reading...\n");
    if (existFF) printf("Form Factors found. Reading...\n");

	if (true)// EII parameters are not currently stored.
	{
		cout << "No rates found. Calculating..." << endl;
		cout << "Total number of configurations: " << dimension << endl;
		Rate Tmp;
		vector<Rate> LocalPhoto(0);
		vector<Rate> LocalFluor(0);
		vector<Rate> LocalAuger(0);
		vector<ffactor> LocalFF(0);
		// Electron impact ionization orbital enerrgy storage.
		EIIdata tmpEIIparams;
		int MaxBindInd = 0;
		// Slippery assumption - electron impact cannot ionize more than the XFEL photon.
		while(Final_occ[MaxBindInd] == orbitals[MaxBindInd].occupancy()) MaxBindInd++;
		tmpEIIparams.kin.clear();
		tmpEIIparams.kin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.ionB.clear();
		tmpEIIparams.ionB.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.fin.clear();
		tmpEIIparams.fin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.occ.clear();
		tmpEIIparams.occ.resize(orbitals.size() - MaxBindInd, 0);
		vector<EIIdata> LocalEIIparams(0);

	    #pragma omp parallel default(none) \
		shared(cout, runlog, MaxBindInd, existAug, existFlr, existPht, existFF) \
		private(Tmp, Max_occ, LocalPhoto, LocalAuger, LocalFluor, LocalFF, LocalEIIparams, tmpEIIparams)
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0; i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				vector<RadialWF> Orbitals = orbitals;
				cout << "configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (int j = 0; j < Orbitals.size(); j++) {
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
				}
				//Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				// Change Lattice to lattice for electron density evaluation.
				Potential U(&lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(lattice, Orbitals, U, input, runlog);

				// EII parameters to store for Later BEB model calculation.
				tmpEIIparams.init = i;
				int size = 0;
				for (int n = MaxBindInd; n < Orbitals.size(); n++) if (Orbitals[n].occupancy() != 0) size++;
				tmpEIIparams.kin = U.Get_Kinetic(Orbitals, MaxBindInd);
				tmpEIIparams.ionB = vector<float>(size, 0);
				tmpEIIparams.fin = vector<int>(size, 0);
				tmpEIIparams.occ = vector<int>(size, 0);
				size = 0;
				//tmpEIIparams.inds.resize(tmpEIIparams.vec2.size(), 0);
				for (int j = MaxBindInd; j < Orbitals.size(); j++) {
					if (Orbitals[j].occupancy() == 0) continue;
					int old_occ = Orbitals[j].occupancy();
					Orbitals[j].set_occupancy(old_occ - 1);
					tmpEIIparams.fin[size] = mapOccInd(Orbitals);
					tmpEIIparams.occ[size] = old_occ;
					Orbitals[j].set_occupancy(old_occ);
					tmpEIIparams.ionB[size] = float(-1*Orbitals[j].Energy);
					tmpEIIparams.kin[size] /= tmpEIIparams.ionB[size];
					size++;
				}
				LocalEIIparams.push_back(tmpEIIparams);

				DecayRates Transit(lattice, Orbitals, u, input);

				// ======= Experimental =========
				if (existFF) LocalFF.push_back({i, Transit.FT_density()});

				Tmp.from = i;

				if (existPht) {
					vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega()/Constant::eV_in_au, runlog);
					for (int k = 0; k < PhotoIon.size(); k++)
					{
						if (PhotoIon[k].val <= 0) continue;
						Tmp.val = PhotoIon[k].val;
						Tmp.to = i + hole_posit[PhotoIon[k].hole];
						Tmp.energy = input.Omega()/Constant::eV_in_au + Orbitals[PhotoIon[k].hole].Energy;
						LocalPhoto.push_back(Tmp);
					}
				}

				if (i != 0)
				{
					if (existFlr) {
						vector<fluor> Fluor = Transit.Fluor();
						for (int k = 0; k < Fluor.size(); k++)
						{
							if (Fluor[k].val <= 0) continue;
							Tmp.val = Fluor[k].val;
							Tmp.to = i - hole_posit[Fluor[k].hole] + hole_posit[Fluor[k].fill];
							Tmp.energy = Orbitals[Fluor[k].fill].Energy - Orbitals[Fluor[k].hole].Energy;
							LocalFluor.push_back(Tmp);
						}
					}

					if (existAug) {
						vector<auger> Auger = Transit.Auger(Max_occ, runlog);
						for (int k = 0; k < Auger.size(); k++)
						{
							if (Auger[k].val <= 0) continue;
							Tmp.val = Auger[k].val;
							Tmp.to = i - hole_posit[Auger[k].hole] + hole_posit[Auger[k].fill] + hole_posit[Auger[k].eject];
							Tmp.energy = Auger[k].energy;
							LocalAuger.push_back(Tmp);
						}
					}
				}
			}

			#pragma omp critical
			{
				Store.Photo.insert(Store.Photo.end(), LocalPhoto.begin(), LocalPhoto.end());
				Store.Fluor.insert(Store.Fluor.end(), LocalFluor.begin(), LocalFluor.end());
				Store.Auger.insert(Store.Auger.end(), LocalAuger.begin(), LocalAuger.end());
				Store.EIIparams.insert(Store.EIIparams.end(), LocalEIIparams.begin(), LocalEIIparams.end());
				Store.FF.insert(Store.FF.end(), LocalFF.begin(), LocalFF.end());
			}
		}       

        // Add the last configuration
        if (existFF) Store.FF.push_back({dimension-1, vector<double> (Store.FF.back().val.size(), 0)});

		sort(Store.Photo.begin(), Store.Photo.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.FF.begin(), Store.FF.end(), [](ffactor A, ffactor B) { return (A.index < B.index); });
		sort(Store.EIIparams.begin(), Store.EIIparams.end(), [](EIIdata A, EIIdata B) {return (A.init < B.init);});
		GenerateRateKeys(Store.Auger);

		if (existPht) {
			string dummy = RateLocation + "Photo.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Photo) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existFlr) {
			string dummy = RateLocation + "Fluor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Fluor) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existPht) {
			string dummy = RateLocation + "Auger.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Auger) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existFF) {
			string dummy = RateLocation + "Form_Factor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& ff : Store.FF) {
				for (int i = 0; i < ff.val.size(); i++) fprintf(fl, "%3.5f ", ff.val[i]);
				fprintf(fl, "\n");
			}
			fclose(fl);
		}
	}

	string IndexTrslt = "./output/" + input.Name() + "/index.txt";
	ofstream config_out(IndexTrslt);
	for (int i = 0; i < Index.size(); i++) {
		config_out << i << " | ";
		for (int j = 0; j < Max_occ.size(); j++) {
			config_out << Max_occ[j] - Index[i][j] << " ";
		}
		config_out << endl;
	}

 	return Store;
}


int RateEquationSolver::SetupAndSolve(ofstream & runlog)
{
	// initial conditions for rate equation
	// first index represents the configuration
	// second index represents the time. The routine will expand it itself unless
	// you provide first adams_n points yourself.
	double fluence = input.Fluence() / Constant::Fluence_in_au;
	double Sigma = input.Width()/Constant::fs_in_au/ (2*sqrt(2*log(2.)));
	int T_size = input.TimePts();
	vector<double> InitCond(dimension, 0);
	InitCond[0] = 1;
	P.clear();

	vector<double> Intensity;
	double scaling_T = 1;

	int converged = 1;
	Plasma Mxwll(T.size());
	while (converged != 0)
	{
		dT.clear();
		dT = generate_dT(T_size);
		T.clear();
		T = generate_T(dT);
		scaling_T = 4 * input.Width()/Constant::fs_in_au / T.back();
		for (int i = 0; i < T.size(); i++) {
			T[i] *= scaling_T;
			dT[i] *= scaling_T;
		}
		Intensity = generate_I(T, fluence, Sigma);
		//SmoothOrigin(T, Intensity);
		IntegrateRateEquation Calc(dT, T, Store, InitCond, Intensity);
		converged = Calc.Solve(0, 1, input.Out_T_size());
		if (converged == 0)	{
			cout << "Final number of time steps: " << T_size << endl;
			P = Calc.GetP();
			T.clear();
			T = Calc.GetT();
			dT.clear();
			dT = vector<double>(T.size(), 0);
            for (int m = 1; m < T.size(); m++) dT[m-1] = T[m] - T[m-1];
            dT[T.size()-1] = dT[T.size()-2];
		} else {
			cout << "Diverged at step: " << converged << " of " << T_size << endl;
			T_size *= 2;
		}
	}


	Intensity.clear();
	Intensity = generate_I(T, fluence, Sigma);
	SmoothOrigin(T, Intensity);

	int tmp = 0;
	for (int i = 0; i < Index.back().size(); i++) tmp += Index.back()[i] - Index.begin()->at(i);
	charge.clear();
	charge.resize(tmp + 1);
	for (auto& ch: charge) ch = vector<double>(T.size(), 0);
	for (int i = 0; i < P.size(); i++)
	{
		tmp = Charge(i);
		for (int m = 0; m < P[i].size(); m++) {
			charge[tmp][m] += P[i][m];
		}
	}

	double t_tmp = 0;
	for (int m = 0; m < T.size(); m++) {
		T[m] = (T[m]-0.5*T.back())*Constant::fs_in_au;
		dT[m] *= Constant::fs_in_au;
	}

	if (input.Write_Charges()) {
		string ChargeName = "./output/Charge_" + input.Name() + ".txt";
		ofstream charge_out(ChargeName);
		double chrg_tmp = 0;
		for (int m = 0; m < T.size(); m++) {
			for (int i = 0; i < charge.size(); i++) {
				chrg_tmp = charge[i][m];
				if (chrg_tmp <= 0.00000001) charge_out << 0 << " ";
				else charge_out << chrg_tmp << " ";
			}
			charge_out << T[m];
			if (m != T.size() - 1) charge_out << endl;
		}

		charge_out.close();
	}

	if (input.Write_Intensity()) {
		string IntensityName = "./output/Intensity_" + input.Name() + ".txt";
		ofstream intensity_out(IntensityName);

		double I_max = *max_element(begin(Intensity), end(Intensity));
		for (int m = 0; m < T.size(); m++) {
			intensity_out << Intensity[m] / I_max << " " << T[m];
			if (m != T.size() - 1) intensity_out << endl;
		}

		intensity_out.close();
	}



	return 0;
}


int RateEquationSolver::SetupAndSolve(MolInp & Input, ofstream & runlog)
{
	// initial conditions for rate equation
	// first index represents the configuration
	// second index represents the time. The routine will expand it itself unless
	// you provide first adams_n points yourself.
	double fluence = Input.Fluence()/ Constant::Fluence_in_au;
	double Sigma = Input.Width()/Constant::fs_in_au / (2*sqrt(2*log(2.)));
	int T_size = Input.ini_T_size();

	P.clear();

	vector<double> Intensity;
    vector<double> old_Intensity;
    vector<double> old_T;
	double scaling_T = 1;

	int converged = 1;
	Plasma Mxwll(T.size());
	while (converged != 0)
	{
		dT.clear();
		dT = generate_dT(T_size);
		T.clear();
		T = generate_T(dT);
		scaling_T = 4*input.Width()/Constant::fs_in_au/T.back();
		for (int i = 0; i < T.size(); i++) {
			T[i] *= scaling_T;
			dT[i] *= scaling_T;
		}

        T_size = T.size();
		Intensity = generate_I(T, fluence, Sigma);
        // Extend mesh to 100 fs (extra 20fs added as pulse starts at -20fs).
        extend_I(Intensity, (2*input.Width() + 100)/Constant::fs_in_au, dT.back());
		
        //SmoothOrigin(T, Intensity);
 		Mxwll.resize(T.size());
		IntegrateRateEquation Calc(dT, T, Input.Store, Mxwll, Intensity);

		converged = Calc.Solve(Mxwll, Input.Store, Input.Out_T_size());

		if (converged == 0) {
			cout << "Final number of time steps: " << T_size << endl;
			P = Calc.GetP();
            old_T = T;
			T = Calc.GetT();
            old_Intensity = Intensity;
		} else {
			cout << "Diverged at step: " << converged << " of " << T_size << endl;
			T_size *= 2;
		}
	}

	Intensity.clear();
    Intensity.resize(T.size(), 0);
    // Intepolate intensity on a new, smaller time grid.
    Interpolation Intp(5);
    for (int m = 0; m < T.size(); m++) Intensity[m] = Intp.get_value(old_Intensity, old_T, T[m])[0];


	double t_tmp = 0;
	double I_max = *max_element(begin(Intensity), end(Intensity));

    // Set the Gaussian pulse peak at T = 0 fs.
	for (int m = 0; m < T.size(); m++) T[m] = T[m]*Constant::fs_in_au - 2*Input.Width();

	int shift = 0;
	vector<int> P_to_charge(0);

    // Aggregate and output charges, plasma parameters, and other parameters into an output.
    // Charges.
    vector<vector<double>> AllAtomCharge(Input.Atomic.size(), vector<double>(T.size(), 0));
    T_size = T.size();
    
	for (int a = 0; a < Input.Atomic.size(); a++) {

		// Occupancies associated with the atom "a".
		vector<int> map_p(Input.Store[a].num_conf);
		for (int i = 0; i < Input.Store[a].num_conf; i++) {
			map_p[i] = i + shift;
		}

		// Work out charges of all configurations.
		double chrg_tmp = 0;
		int tmp = Input.Index[a][0].size();
		P_to_charge.clear();
		P_to_charge.resize(Input.Store[a].num_conf, 0);
		for (int i = 0; i < Input.Index[a].size(); i++) {
			for (int j = 0; j < tmp; j++) P_to_charge[i] += Input.Index[a][i][j];
			chrg_tmp = 0;
		}

		charge.clear();
		charge.resize(*max_element(begin(P_to_charge), end(P_to_charge)) + 1);
		for (auto& ch: charge) ch = vector<double>(T.size(), 0);
		for (int i = 0; i < map_p.size(); i++)
		{
            int conf_ind = map_p[i];
			tmp = P_to_charge[i];
			for (int m = 0; m < P[i].size(); m++) charge[tmp][m] += P[conf_ind][m];
		}

        for (int m = 0; m < T.size(); m++) {
            for (int i = 1; i < charge.size(); i++) {
                AllAtomCharge[a][m] += i*charge[i][m];
            }
        }

		shift += Input.Store[a].num_conf;

		if (Input.Write_Charges()) {
			string ChargeName = "./output/Charge_" + Input.Store[a].name + ".txt";
			ofstream charge_out(ChargeName);
			double chrg_tmp = 0;
			for (int m = 0; m < T.size(); m++) {
				for (int i = 0; i < charge.size(); i++) {
					chrg_tmp = charge[i][m];
					if (chrg_tmp <= 0.00000001) charge_out << 0 << " ";
					else charge_out << chrg_tmp << " ";
				}
				charge_out << T[m];
				if (m != T.size() - 1) charge_out << endl;
			}

			charge_out.close();
		}

	}

	if (Input.Write_MD_data()) {
		string MD_Data = "./output/MD_Data.txt";
		ofstream OutFile(MD_Data);
		OutFile << T.size() << endl;
		OutFile << "Time ";
		for (int a = 0; a < Input.Atomic.size(); a++) OutFile << Input.Atomic[a].Nuclear_Z() << " ";
		OutFile << "N(elec) E(elec)";
		for (int m = 0; m < T.size(); m++) {
			OutFile << endl << T[m] << " ";
			for (int a = 0; a < AllAtomCharge.size(); a++) {
				OutFile << AllAtomCharge[a][m] << " ";
			}
			if (m != 0) OutFile << Mxwll.N[m] << " " << Mxwll.E[m];
			else OutFile << 0 << " " << 0;
		}

		OutFile.close();
	}


	if (Input.Write_Intensity()) {
		string IntensityName = "./output/Intensity_" + Input.name + ".txt";
		ofstream intensity_out(IntensityName);

		double I_max = *max_element(begin(Intensity), end(Intensity));
		for (int m = 0; m < T.size(); m++) {
			intensity_out << Intensity[m] / I_max << " " << T[m];
			if (m != T.size() - 1) intensity_out << endl;
		}

		intensity_out.close();
	}

    
    if (Input.Write_FormFactor()) {
        // Aggregate form factor for each atom, for each time point
        // and output to a file.
        string FFName = "./output/FormFactor_" + Input.name + ".txt";
		ofstream FF_out(FFName);
        shift = 0;

        // Read Form-Factors for each configuration.
        double agg_FF[20];

		for (int a = 0; a < Input.Atomic.size(); a++) {
            // Occupancies associated with the atom "a".
            int num_conf = Input.Store[a].num_conf;
            
            vector<double> p_t(num_conf);
            for (int m = 0; m < T.size(); m++) {
                for (int q = 0; q < 20; q++) agg_FF[q] = 0;
                for (int i = 0; i < num_conf; i++) p_t[i] = P[i + shift][m];            

                FF_out << T[m] << " " << Input.Store[a].name;
                for (auto& conf_ff: Input.Store[a].FF) {
                    int conf_ind = conf_ff.index;
                    double conf_prob = p_t[conf_ind];

                    for (int q = 0; q < 20; q++) agg_FF[q] += conf_prob * conf_ff.val[q];
                }
                for (int q = 0; q < 20; q++) FF_out << " " << agg_FF[q];
                FF_out << "\n";
            }

            shift += num_conf;
        }

		FF_out.close();
    }

	return 0;
}




bool RateEquationSolver::ReadRates(const string & input, vector<Rate> & PutHere)
{
	if (exists_test(input))
	{
		Rate Tmp;
		ifstream infile(input);

		char type;
		while (!infile.eof())
		{
			string line;
			getline(infile, line);

			stringstream stream(line);
			stream >> Tmp.val >> Tmp.from >> Tmp.to >> Tmp.energy;

			PutHere.push_back(Tmp);
		}

		infile.close();
		return true;
	}
	else return false;
}


bool RateEquationSolver::ReadFFactors(const string & input, vector<ffactor> & PutHere)
{
	if (exists_test(input))
	{
		ffactor Tmp;
        Tmp.val.clear();
        Tmp.val.resize(20, 0);
		ifstream ifs(input);
        string line;

        int i = 0;
		while (!ifs.eof())
		{
			getline(ifs, line);
            if (line == "") continue;
			stringstream ss(line);
            Tmp.index = i;
            
            // TODO: Q_mesh hardcoded to have 20 points.
            for (int q = 0; q < 20; q++) ss >> Tmp.val[q];

			PutHere.push_back(Tmp);
            i++;
		}

		ifs.close();
		return true;
	}
	else return false;
}


int RateEquationSolver::Symbolic(const string & input, const string & output)
{
	if (Store.Photo.size() == 0)
	{
		if (exists_test(input)) {
			ifstream Rates_in(input);
			ofstream Rates_out(output);

			Rate Tmp;
			char type;

			while (!Rates_in.eof())
			{
				string line;
				getline(Rates_in, line);

				stringstream stream(line);
				stream >> Tmp.val >> Tmp.from >> Tmp.to >> Tmp.energy;

				Rates_out << Tmp.val << " " << InterpretIndex(Tmp.from) << " " << InterpretIndex(Tmp.to) << " " << Tmp.energy << endl;
			}

			Rates_in.close();
			Rates_out.close();

			return 0;
		}
		else return 1;
	} else {
		ofstream Rates_out(output);

		for (auto& v : Store.Photo) {
			Rates_out << v.val << " " << InterpretIndex(v.from) << " " << InterpretIndex(v.to) << " " << v.energy << endl;
		}

		Rates_out.close();

		return 0;
	}

}


string RateEquationSolver::InterpretIndex(int i)
{
	// LaTeX output format for direct insertion into table.
	string Result;
	if (!Index.empty()) {
		if (Index[i].size() == orbitals.size())	{
			ostringstream tmpstr;
			tmpstr << '$';
			for (int j = 0; j < orbitals.size(); j++) {
				tmpstr << orbitals[j].N();
				switch (orbitals[j].L()) {
				case 0:
					tmpstr << 's';
					break;
				case 1:
					tmpstr << 'p';
					break;
				case 2:
					tmpstr << 'd';
					break;
				case 3:
					tmpstr << 'f';
					break;
				default:
					tmpstr << 'x';
					break;
				}
				tmpstr << '^' << orbitals[j].occupancy() - Index[i][j];
			}
			tmpstr << '$';
			Result = tmpstr.str();
		}
	}
	return Result;
}


int RateEquationSolver::Charge(int Iconf)
{
	int Result = 0;
	for (int j = 0; j < Index[Iconf].size(); j++) Result += Index[Iconf][j] - Index[0][j];
	return Result;
}


bool RateEquationSolver::SetupIndex(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	if (orbitals.size() != Final_occ.size())
	{
		runlog << "Final occupancies should be provided for all orbitals." << endl;
		return false;
	}
	if (Max_occ.size() != orbitals.size()) {
		Max_occ.resize(orbitals.size());
		for (int i = 0; i < orbitals.size(); i++) {
			Max_occ[i] = 4*orbitals[i].L() + 2;
		}
	}
	// Work out the number of allowed configurations
	dimension = 1;
	int orbitals_size = orbitals.size();
	int l_hole = 0;//lowest energy orbital with allowed holes
	hole_posit.clear();
	hole_posit.resize(orbitals.size());

	vector<int> max_holes(orbitals.size(), 0);
	for (int i = orbitals_size - 1; i >= 0; i--)
	{
		max_holes[i] = Max_occ[i] - Final_occ[i] + 1;
		if (Max_occ[i] == Final_occ[i]) l_hole++;
		hole_posit[i] = dimension;
		dimension *= max_holes[i];
	}

	//configuration information. Index[i][j] denotes number of holes in oorbital [j], that together form configuration [i].
	Index.resize(dimension);
	for (auto& v : Index) {
		v.resize(orbitals.size());
	}

	for (int i = 0; i < dimension; i++)
	{
		int tmp = i;
		for (int j = 0; j < orbitals_size; j++)
		{
			Index[i][j] = tmp / hole_posit[j];
			if (Index[i][j] > Max_occ[j]) { Index[i][j] = Max_occ[j]; }
			tmp -= Index[i][j] * hole_posit[j];
		}
	}
	return true;
}


RateEquationSolver::~RateEquationSolver() {}


inline bool exists_test(const std::string& name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}


vector<double> RateEquationSolver::generate_dT(int num_elem)//default time interval
{
	vector<double> Result(num_elem, 0);
	double tmp = 1;
	for (int i = 0; i < num_elem; i++) {
	tmp = fabs(1.*i / (num_elem-1) - 0.5) + 0.01;
		Result[i] = tmp;
	}
	return Result;
}


vector<double> RateEquationSolver::generate_T(vector<double>& dT)//default time
{
	vector<double> Result(dT.size(), 0);
	vector<double> Bashforth_4{ 55. / 24., -59. / 24., 37. / 24., -9. / 24. }; //Adams�Bashforth method
	for (int i = 1; i < Bashforth_4.size(); i++)//initial few points
	{
		Result[i] = Result[i - 1] + dT[i - 1];
	}
	for (int i = Bashforth_4.size(); i < Result.size(); i++)//subsequent points
	{
		Result[i] = Result[i - 1];
		for (int j = 0; j < Bashforth_4.size(); j++)
		{
			Result[i] += Bashforth_4[j] * dT[i - j - 1];
		}
	}

	return Result;
}


vector<double> RateEquationSolver::generate_I(vector<double>& Time, double Fluence, double Sigma)//intensity of the Gaussian X-ray pulse
{

	vector<double> Result(Time.size(), 0);

	double midpoint = 0.5*(T.back() + T[0]);
	double denom = 2*Sigma*Sigma;
	double norm = 1./sqrt(denom*Constant::Pi);
	//include the window function to make it exactly 0 at the beginning and smoothly increase toward Gaussian


  int smooth = T.size()/10;
  double tmp = 0;
	for (int i = 0; i < Time.size(); i++)
	{
    Result[i] = Fluence * norm * exp(-(Time[i] - midpoint)*(Time[i] - midpoint) / denom);
    if (i < smooth) {
      tmp = fabs(T[i] - T[0]) / fabs(T[smooth] - T[0]);
      Result[i] *= tmp*tmp*(3 - 2 * tmp);
    }
  }

  /*
  for (int i = 0; i < Time.size(); i++)
	{
    Result[i] = Fluence / T.back();
	}
  */



	return Result;
}


int RateEquationSolver::extend_I(vector<double>& Intensity, double new_max_T, double step_T)
{
  // Uniform mesh is added.
  double last_T = T.back(), dT_last = dT.back();
  //double MidT = 0.5*(T.back() - T[0]);
  //double denom = 2*Sigma*Sigma;
  //double norm = 1./sqrt(denom*Constant::Pi);

  while (last_T < new_max_T) {
    dT.push_back(dT_last);
    T.push_back(last_T + dT_last);
    last_T = T.back();
    Intensity.push_back(0);
  }
  return 0;
}


void SmoothOrigin(vector<double> & T, vector<double> & F)
{
	int smooth = T.size() / 10;
	for (int i = 0; i < smooth; i++)
	{
		F[i] *= (T[i] / T[smooth])*(T[i] / T[smooth])*(3 - 2 * (T[i] / T[smooth]));
	}
}


vector<double> RateEquationSolver::generate_G()
{
	// Intensity profile normalized to 1.
	// Time is assumbed to be in FEM
	double Sigma = input.Width()/(2*sqrt(2*log(2.)));//Constant::fs_in_au

	return generate_I(T, 1, Sigma);
}


void RateEquationSolver::GenerateRateKeys(vector<Rate> & ToSort)
{
	int CurrentFrom = 0;
	int start = 0;
	RatesFromKeys.push_back(0);
	for (int i = 1; i < ToSort.size(); i++) {
		if (ToSort[i].from != CurrentFrom) {
			CurrentFrom = ToSort[i].from;
			start = RatesFromKeys.back();
			RatesFromKeys.push_back(i);
			sort(ToSort.begin() + start, ToSort.begin() + RatesFromKeys.back(), sortRatesTo);
		}
	}
}


int RateEquationSolver::mapOccInd(vector<RadialWF> & Orbitals)
{
	int Result = 0;
	for (int j = 0; j < hole_posit.size(); j++)	{
		Result += (orbitals[j].occupancy() - Orbitals[j].occupancy())*hole_posit[j];
	}

	return Result;
}


double RateEquationSolver::T_avg_RMS(vector<pair<double, int>> conf_RMS)
{
  // Calculate pulse-averaged root mean square radius of an atom.
  double tmp = 0;

  vector<double> intensity = generate_G();
  if (P.size()-1 != density.size()) return -1;
  for (int m = 0; m < T.size(); m++) {
    tmp = 0;
    for (int i = 0; i < conf_RMS.size(); i++) tmp += P[i][m]*conf_RMS[i].first;
    intensity[m] *= tmp;
  }

  Grid Time(T, dT);
  Adams I(Time, 10);

  return I.Integrate(&intensity, 0, T.size()-1);
}


double RateEquationSolver::T_avg_Charge()
{
  // Calculate pulse-averaged charge of atom.
  double tmp = 0;

  vector<double> intensity = generate_G();
  for (int m = 0; m < T.size(); m++) {
    tmp = 0;
    for (int i = 0; i < charge.size(); i++) tmp += (input.Nuclear_Z() - i)*charge[i][m];
    intensity[m] *= tmp;
  }

  Grid Time(T, dT);
  Adams I(Time, 10);

  return I.Integrate(&intensity, 0, T.size()-1);
}