#pragma once

#include "RadialWF.h"
#include "Grid.h"
#include "Potential.h"
#include <vector>
#include "Constant.h"
#include "IntegrateRateEquation.h"
#include "Input.h"

using namespace std;

//Class for solving system of linear coupled diffecrential equations of the form
//
// dP(i)/dt = Sum{j!=i} (R(j->i)P(j) - R(i->j)P(i))
//
//where the rates R(j->i) represent electronic transitions from configuration j to
//configuration i (Photo-ionization, fluorescense, and Auger). Photoionization im
//external beam with a given flux I(t) is assumed. Assumption - electron in the
//continuum orbital is out of the system.


class RateEquationSolver
{
public:
	RateEquationSolver(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, Input & Inp);
	~RateEquationSolver();

	// Halfwidth = 5/Constant::Time -> 5 fs half width.
	int SolveFrozen(vector<int> Max_occ, vector<int> Final_occ, ofstream & log);
	int SolveFrozen(vector<RadialWF> & Virtual, vector<int> Max_occ, vector<int> Final_occ, ofstream & log);
	AtomRateData SolvePlasmaBEB(vector<int> Max_occ, vector<int> Final_occ, ofstream & log, vector<bool> args = vector<bool>(3, false));
	// Atomic.
	int SetupAndSolve(ofstream & log, int out_T_size = 500);
	// Molecular.
	int SetupAndSolve(MolInp & Input, ofstream & log, int out_T_size = 500);

	//string CompareRates(string RateFile1, string RateFile2, ofstream & log);// Find the difference in rate equation using two different rates.

	bool ReadRates(const string & input, vector<Rate> & PutHere);
	int Symbolic(const string & input, const string & output);//convertes configuration indexes in human readable format
	int Charge(int Iconf);
	vector<double> PerturbMe(vector<RadialWF> & Virtual, double Dist, double Einit);
	vector<double> Secular(vector<RadialWF> & Virtual, double Dist, double Einit);
	
	int NumPath() { return dimension; }
	vector<double> generate_G();
	vector<double> Times() { return T; }
	vector<double> dTimes() { return dT; }
	vector<double> Probs(int i) { return P[i]; }
	vector<vector<double>> AllProbs() {return P;}


	bool SetupIndex(vector<int> Max_occ, vector<int> Final_occ, ofstream & log);
	vector<vector<int>> Get_Indexes() { return Index; }

  // Atomic data containers.
	vector<vector<double>> density = vector<vector<double>>(0);
  AuxAtomData AtomAuxStore;

  Grid & Atom_Mesh() { return lattice; }
protected:
	Grid & lattice;
	Input & input;
	vector<RadialWF> & orbitals;
	Potential& u;
	vector<CustomDataType::polarize> MixMe; 
	int dimension;//number of configurations
	vector<vector<double>> charge;
	vector<double> T;// Time grid points.
	vector<double> dT;// Accurate differentials.
	vector<vector<double>> P;// P[i][m] is the probabilities of having configurations "i" at time T[m].
	vector<vector<int> > Index;
	int mapOccInd(vector<RadialWF> & Orbitals);// Inverse of what Index returns.

	//int SetupAndSolve(vector<Rate> rates, double I_max, double HalfWidth, ofstream & log, int & start_T_size);

private:
	string InterpretIndex(int i);

	AtomRateData Store;
	
	vector<CustomDataType::ffactor> FF;
	vector<int> hole_posit;
  
	int extend_I(vector<double>& Intensity, double new_max_T, double step_T);
  vector<double> generate_I(vector<double>& T, double I_max, double HalfWidth);
	vector<double> generate_T(vector<double>& dT);
	vector<double> generate_dT(int num_elem);
  double T_avg_RMS(vector<pair<double, int>> conf_RMS);

	static bool sortEIIbyInd(CustomDataType::EIIdata A, CustomDataType::EIIdata B) { return (A.init < B.init); } 
	static bool sortRatesFrom(Rate A, Rate B) { return (A.from < B.from); }
	static bool sortRatesTo(Rate A, Rate B) { return (A.to < B.to); }
	// Keys allow to quickly find the required element. See the GenerateFromKeys().
	vector<int> RatesFromKeys;
	void GenerateRateKeys(vector<Rate> & ToSort); 
};

