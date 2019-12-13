#ifndef SPLINE_H
#define SPLINE_H

#include <vector>

class Bspline
{
	/* This class defines and calculates B-splines at a given knot sequence. */
public:
	Bspline(std::vector<double>* Knots);
	Bspline(std::vector<double>* Knots, int Order);
	~Bspline() {}

	double Value(unsigned int k, unsigned int i, double x);
	double dValue(unsigned int k, unsigned int i, unsigned int d_Order, double x); //d_Order derivative of B_{i, k} at point "x"
	std::vector<double> Values(unsigned int k, unsigned int d_Order, double x); // calculates non-trivial B-spline values at "x" if d_Order=0. Otherwise returns derivatives of order d_Order. 
	std::vector<int> Values_i() { return spln_nums; } // addition to Values(...). Returnes positions of non-trivial Bsplines at "x" in B-spline set. Can be run only after Values.
	
	int Size() { return abs(size - order); }
	double Knot(int i);
	int Order() { return order; }

protected:
	std::vector<double>* knots;
	std::vector<double> values;  // Non-trivial B-spline values at "x"
	std::vector<int> spln_nums; // Indexes of above B-splines,  SplnNumbers[i] = j & values[i] = B_j(x) !=0
	int order;
	int size;
};
	

#endif