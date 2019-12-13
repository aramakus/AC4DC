#pragma once

#include "RadialWF.h"
#include "Grid.h"
#include "Potential.h"
#include <vector>
#include "Constant.h"
#include <fstream>
#include "Input.h"

using namespace std;
using namespace CustomDataType;
//Class calculates total photoionization crossection for a given orbital

class DecayRates
{
public:
	DecayRates(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, Input & Input);

	vector<photo> Photo_Ion(double omega, ofstream & log); // All photoinonization crossections. Position in vector indicates orbital
	vector<fluor> Fluor(); // Fluorescence rates for all channels.
	vector<auger> Auger(vector<int> Max_occ, ofstream & log); // Auger decay rates for all channels.

	~DecayRates();

private:
	int IntegrateContinuum(Grid &Lattice, Potential &U, vector<RadialWF> &Core, RadialWF* Current, int c = 0);
	Grid& lattice;
	vector<RadialWF>& orbitals;
	Potential& u;
	Input & input;

	// EII internal intepolated data.
	Grid CntLattice = Grid(0);// Grid for continuum wave calculations.
	Potential CntU;
	vector<RadialWF> CntOrbitals;// Interpolated onto continuum grid "orbitals".
};

class FormFactor
{
//======================================================================
//	Least square fitting of electronic density
//	with analytic functions r^(l+1)*exp(-alpha_i*r). Used to calculate
//	scattering formfactors, which are Fourier transforms of the
//	density.
//======================================================================
public:
	FormFactor(Grid & Lattice, vector<double> & density, int infinity);//l sets asymptotic for small r (as r^(l_asymp)), num_func is number of basis functions
	~FormFactor();

	//density[infinity] < 10^-15
	double getFF(double Q); //returns least square expansion coefficients
	vector<double> getAllFF();
	vector<double> getAllQ() { return Q_mesh; }
protected:
	Grid lattice;//dense lattice
	vector<double> iDensity;//interpolated density
	double Q0 = 0;
	vector<double> Q_mesh;
};
