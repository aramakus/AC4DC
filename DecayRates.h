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