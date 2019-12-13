#ifndef SPLINEINTEGR_H
#define SPLINEINTEGR_H

#include "Spline.h"
#include "GaussQadrature.h"
#include <vector>

class BsplineIntegrals
{
	/* this class calculates Hamiltonian and Overlap matrix and
	   passes if to LinearAlgebra to calculate Hydrogen atom states.
	   Update for Hartree-fock later.
	*/
public:
	BsplineIntegrals(Bspline &Bas);
	~BsplineIntegrals() {}

	double Overlap(int i, int j) { return overlap[i][j]; }
	double Kinetic(int i, int j) { return kinetic[i][j]; }
	double HydroPotential(int i, int j) { return hydro_potential[i][j]; }

protected:
	std::vector<std::vector<double>> overlap;
	std::vector<std::vector<double>> kinetic;
	std::vector<std::vector<double>> hydro_potential;
};


#endif