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
#pragma once

#include <vector>
#include "Grid.h"
#include "RadialWF.h"
#include <cmath>

using namespace std;
/*this class deploys the Adams routine for solution of the linear coupled system of ODE dy/dx = f(x)*y
y is a vector with two components (F G), and f is a 2x2 matrix, so the equation take the following form

 dF/dx = A * F + B * G + X
 dG/dx = C * F + D * G + Y

Also ordinary integration can be performed using overloaded Integrate function
*/
class Adams
{
public:
	Adams(Grid &Lattice, int AdamsOrder);
	~Adams(void);

	vector<double> A, B, C, D, X, Y;

	// Start from 0 and calculare 9 subsequent points if forward=true, 
  // or start from the end and calculate 9 previous points if forward=false.
  // Gets initial points for Adams integrator.
	void StartAdams(RadialWF* Psi, int start_pt, bool forward);

	//finds maximum of Psi.F before classical turning point
	int FindMaximum(int Turn)
	{
		if (FirstMaxima == 0) { return Turn; }
		else { return FirstMaxima; }
	}

	//finds maximum closest to R_box

	void Integrate(RadialWF* Psi, int start_pt, int end_pt);
	void Integrate(vector<double> &Func, vector<double> &Result, int start_pt, int end_pt);
	double Integrate(vector<double> * Func, int start_pt, int end_pt);
	int Nodes() { return NumNodes;  }
	vector<double> GreenOrigin(RadialWF* Psi);
	vector<double> GreenInfinity(RadialWF* Psi);

	void Integrate_ODE(vector<double> &f, int start_pt, int end_pt);

	int GetAdamsOrder() { return Adams_N; }

protected:
	int Adams_N;

	int FirstMaxima = 0;
	int NumNodes = 0;

	Grid& Lattice;

	vector<double> Adams_Coeff;
};



//Interpolation routine from Numerical Receipts. Polynomial interpolation with Neville's algorithm
class Interpolation
{
public:
	Interpolation(int Order = 6);
	~Interpolation() {};

	vector<double> get_value(const vector<double> &f, const vector<double> &x_ini, double X);//first argument is the value, second is the derivative at X
	vector<double> get_value(PairFunction &S_old, Grid &Lattice_old, double X);//first argument is the value, second is the derivative at X
	int RecalcWF(RadialWF &S_old, Grid &Lattice_old, RadialWF &S_new, Grid &Lattice_new);
  
  // Interpolate the function Vals, defined on a grid R (with intervals dR) with a sum of "Order" Gaussians. Weights and exponenets are optimized.
  void gaussian_sum(vector<double> & Vals, Grid & Lattice, int max_iter = 1000, double conv_toll = 0.001);

protected:
	int order;
};



class GaussQuad
{
public:
	GaussQuad(int Order = 6);
	~GaussQuad() {};

	void set_order(int new_Order);
	vector<double> get_Gauss_X(double a, double b);
	double Integrate(vector<double> &F, vector<double> &x, double a = -1., double b = 1.);
private:
	int order;
	vector<double> GaussW;
	vector<double> GaussX;
};

