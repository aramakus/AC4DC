#pragma once
#include <vector>
#include <list>
#include <fstream>
#include "Constant.h"
#include "Plasma.h"

using namespace std;

struct Rate
{
	double val = 0;
	long int from = 0;
	long int to = 0;
	double energy = 0;
};

struct AtomRateData
{
	string name = "";
	double nAtoms = 1.;// atomic number density
	double R = 189.; // 100nm focal spot radius.
	int num_conf = 1;
	vector<Rate> Photo = vector<Rate>(0);
	vector<Rate> Fluor = vector<Rate>(0);
	vector<Rate> Auger = vector<Rate>(0);
	vector<CustomDataType::EIIdata> EIIparams = vector<CustomDataType::EIIdata>(0);
};

class IntegrateRateEquation
{
	vector<vector<double>> dpdt;
	vector<vector<double>> p;
	vector<vector<double>> p_storage;

	vector<double> time_storage;
	vector<double>& t;
	vector<double>& dt;
	vector<double> f;
	AtomRateData& store;

	int adams_n = 5;
public:
	// Rate equations for single atom. 
	IntegrateRateEquation(vector<double>& dT, vector<double>& T, AtomRateData& Store, vector<double> InitCond, const vector<double>& Intensity = vector<double>());
	int Solve(double P_min = 0, double P_max = 1, int storage_time_pts = 500);
	// Rate equations for single chemical element + electron plasma.
	IntegrateRateEquation(vector<double>& dT, vector<double>& T, AtomRateData& Store, Plasma & Electrons, vector<double> InitCond, const vector<double>& Intensity = vector<double>());
	int Solve(Plasma & Electrons, double P_min = 0, double P_max = 1, int storage_time_pts = 500);
	// Rate equations for molecule + electron plasma.
	IntegrateRateEquation(vector<double>& dT, vector<double>& T, vector<AtomRateData> & Store, Plasma & Electrons, const vector<double>& Intensity = vector<double>());
	int Solve(Plasma & Electrons, vector<AtomRateData> & Store, int storage_time_pts = 500);
	int WriteCharge(vector<AtomRateData> & Store);
	// P_min and P_max are upper and lower bound for P.
	// storage_time_pts is the number of time points to store. Sometimes calculations might require
	// too many time points, so it is cheaper to store some and interpolate later if needed.
	
	int Write(ofstream & charge);

	vector<vector<double>> GetP() { return p_storage; }
	vector<double> GetT() { return time_storage; }

	~IntegrateRateEquation();
};


class MonteCarlo
{
	// Monte Carlo implementation of rate equation solver. 
};

